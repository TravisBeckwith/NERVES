#!/usr/bin/env python3
"""
nerves_qc.py
============
Tensor-fit health and scalar-reliability QC for diffusion derivatives.

Reads tensor scalar maps that already exist (QSIRecon workflows such as
DSI Studio / DIPY / TORTOISE tensor fits, or NERVES/FSL ``dti_*`` outputs),
and for every subject/workflow reports

* **fit health**: fraction of non-physical voxels, FA saturation, MD
  plausibility, internal consistency (MD vs (AD+2RD)/3, FA vs AD/RD);
* **expected error** of FA / MD / mode / NDE at the subject's SNR
  (Monte-Carlo, see ``nerves_reliability.py``) when an SNR is known;
* **group outliers** (robust z-scores) when enough subjects are present.

    python nerves_qc.py DERIVATIVES_DIR [--out qc_out] [--snr 25]
                        [--dirs 30] [--b 1000]

Outputs in ``--out``: ``nerves_qc_summary.tsv``, ``nerves_qc_summary.json``,
``nerves_qc_report.html``. Exit code: 0 all pass, 3 if any subject fails a
hard check, 1 if nothing was found.

The brain mask is inferred from voxels where MD is finite and > 0 (QSIRecon
maps are zero outside the brain mask). No registration or atlas is needed;
all checks are voxelwise on the maps as written.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import nibabel as nib

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nerves_dti as nd  # noqa: E402
import nerves_reliability as rel  # noqa: E402

__version__ = "2.0.0-dev"

# --- thresholds (documented in docs/reliability.md) -----------------------
FAIL_NONPHYSICAL_FRAC = 0.05     # > 5 % non-physical voxels -> FAIL
WARN_NONPHYSICAL_FRAC = 0.01
WARN_FA_SATURATED_FRAC = 0.02    # FA >= 0.99 in > 2 % of voxels
WARN_CONSISTENCY = 0.05          # median |MD-(AD+2RD)/3|/MD
WARN_EXPECTED_FA_BIAS = 0.05     # expected GM-like FA bias at this SNR
WARN_EXPECTED_MODE_SD = 0.5      # GM-like mode SD: mode not interpretable
OUTLIER_Z = 3.5
MIN_GROUP_FOR_OUTLIERS = 5
MD_MAX_PLAUSIBLE = 4.0           # 1e-3 mm^2/s; free water is ~3

QSIRECON_RE = re.compile(
    r"^(?P<prefix>.+)_param-(?P<param>[A-Za-z0-9]+)_dwimap\.nii(\.gz)?$")
NII_SUFFIX = re.compile(r"\.nii(\.gz)?$")
FSL_NAMES = {"dti_FA": "fa", "dti_MD": "md", "dti_L1": "l1",
             "dti_L2": "l2", "dti_L3": "l3"}
KNOWN = {"fa", "md", "ad", "rd", "l1", "l2", "l3"}


# ---------------------------------------------------------------- discovery
def discover(root: Path) -> Dict[str, Dict[str, Path]]:
    """Map ``label -> {param: path}`` for every tensor scalar set under
    ``root``. QSIRecon files group by (workflow folder, filename prefix);
    FSL/NERVES ``dti_*`` files group by directory."""
    root = Path(root)
    found: Dict[str, Dict[str, Path]] = {}
    for p in sorted(root.rglob("*.nii*")):
        name = p.name
        m = QSIRECON_RE.match(name)
        if m and m.group("param").lower() in KNOWN:
            rp = p.parent.relative_to(root).parts
            wf = next((x for x in rp if x.startswith("qsirecon")), "")
            label = f"{wf}/{m.group('prefix')}" if wf else m.group("prefix")
            found.setdefault(label, {})[m.group("param").lower()] = p
            continue
        stem = NII_SUFFIX.sub("", name)
        if stem in FSL_NAMES:
            label = str(p.parent.relative_to(root)) or "."
            found.setdefault(label, {})[FSL_NAMES[stem]] = p
    return {k: v for k, v in found.items() if "fa" in v and "md" in v}


def _load(path: Path) -> np.ndarray:
    return np.asarray(nib.load(str(path)).dataobj, dtype=np.float64)


def _find_report(paths: Dict[str, Path], root: Path) -> Optional[dict]:
    """NERVES report next to the dti_* files, if any (for SNR / NDE stats)."""
    d = next(iter(paths.values())).parent
    for cand in (d, d.parent, d.parent.parent):
        if cand == root.parent:
            break
        r = cand / "nde_pipeline_report.json"
        if r.exists():
            try:
                return json.loads(r.read_text())
            except (OSError, ValueError):
                return None
    return None


# ------------------------------------------------------------ per-subject QC
def analyse_subject(paths: Dict[str, Path], snr: Optional[float] = None,
                    n_dirs: int = 30, b: float = 1000.0,
                    method: str = "ols") -> dict:
    """Fit-health, consistency and expected-error metrics for one set of
    tensor scalar maps. Returns a flat dict (``flags`` is a list)."""
    fa = _load(paths["fa"])
    md = _load(paths["md"])
    if fa.shape != md.shape:
        raise ValueError(f"FA {fa.shape} and MD {md.shape} grids differ")

    inb = np.isfinite(md) & (md != 0)                  # inferred brain mask
    n = int(inb.sum())
    if n == 0:
        return {"n_voxels": 0, "status": "FAIL",
                "flags": ["empty: no non-zero MD voxels"]}

    # units: mm^2/s (~7e-4) -> 1e-3 mm^2/s
    scale = 1e3 if np.nanmedian(np.abs(md[inb])) < 0.05 else 1.0
    mdv = md[inb] * scale
    fav = fa[inb]
    flags: List[str] = []
    row: dict = {"n_voxels": n, "md_units_scaled": scale != 1.0}

    bad_fa = ~np.isfinite(fav) | (fav < -1e-6) | (fav > 1 + 1e-6)
    nonphys = bad_fa | ~np.isfinite(mdv) | (mdv <= 0)

    have_axial = "ad" in paths and "rd" in paths
    if have_axial:
        adv = _load(paths["ad"])[inb] * scale
        rdv = _load(paths["rd"])[inb] * scale
        nonphys |= ~np.isfinite(adv) | ~np.isfinite(rdv) | (rdv < 0) | \
            (adv < rdv - 1e-6 * np.maximum(np.abs(adv), 1e-12))
        good = ~nonphys & (mdv > 0)
        if good.any():
            resid = np.abs(mdv[good] - (adv[good] + 2 * rdv[good]) / 3) \
                / mdv[good]
            row["consistency_md_ad_rd"] = float(np.median(resid))
            # FA implied by (AD, RD) under axial symmetry
            ad_, rd_ = adv[good], rdv[good]
            den = np.sqrt(ad_ ** 2 + 2 * rd_ ** 2)
            fa_imp = np.sqrt(1.5) * np.sqrt(
                (ad_ - mdv[good]) ** 2 + 2 * (rd_ - mdv[good]) ** 2) \
                / np.where(den > 0, den, np.nan)
            row["fa_vs_axial_fa_median_abs_diff"] = float(
                np.nanmedian(np.abs(fav[good] - fa_imp)))
        if row.get("consistency_md_ad_rd", 0) > WARN_CONSISTENCY:
            flags.append("MD inconsistent with (AD+2RD)/3: median error "
                         f"{row['consistency_md_ad_rd']:.1%}")

    if all(k in paths for k in ("l1", "l2", "l3")):
        e = np.stack([_load(paths[k])[inb] * scale
                      for k in ("l1", "l2", "l3")], axis=1)
        neg = (e <= 0).any(axis=1)
        order = ~((e[:, 0] >= e[:, 1]) & (e[:, 1] >= e[:, 2]))
        nonphys |= neg | ~np.isfinite(e).all(axis=1)
        row["frac_negative_eigenvalue"] = float(neg.mean())
        row["frac_misordered_eigenvalues"] = float(order.mean())
        ok = ~nonphys
        if ok.any():
            nde = nd.nde_from_eigenvalues(e[ok])
            mode = nd.tensor_mode(e[ok])
            row["nde_median"] = float(np.median(nde))
            row["mode_median"] = float(np.nanmedian(mode))
            row["frac_mode_undefined"] = float(np.isnan(mode).mean())

    row["frac_nonphysical"] = float(nonphys.mean())
    valid = ~nonphys
    row["frac_fa_saturated"] = float((fav[valid] >= 0.99).mean()) \
        if valid.any() else float("nan")
    row["frac_md_implausible"] = float(
        (mdv[valid] > MD_MAX_PLAUSIBLE).mean()) if valid.any() else float("nan")
    if valid.any():
        row["md_median"] = float(np.median(mdv[valid]))
        row["fa_mean"] = float(fav[valid].mean())
        row["fa_p95"] = float(np.percentile(fav[valid], 95))
        row["frac_fa_gt_0p4"] = float((fav[valid] > 0.4).mean())

    # ---- hard / soft flags
    status = "PASS"
    if row["frac_nonphysical"] > FAIL_NONPHYSICAL_FRAC:
        status = "FAIL"
        flags.append(f"{row['frac_nonphysical']:.1%} non-physical voxels "
                     f"(> {FAIL_NONPHYSICAL_FRAC:.0%})")
    elif row["frac_nonphysical"] > WARN_NONPHYSICAL_FRAC:
        flags.append(f"{row['frac_nonphysical']:.1%} non-physical voxels")
    if row.get("frac_fa_saturated", 0) > WARN_FA_SATURATED_FRAC:
        flags.append(f"FA >= 0.99 in {row['frac_fa_saturated']:.1%} of "
                     "voxels (noise floor / fit failure)")
    if row.get("frac_md_implausible", 0) > 0.05:
        flags.append(f"MD > {MD_MAX_PLAUSIBLE:g}e-3 in "
                     f"{row['frac_md_implausible']:.1%} of voxels")
    if row.get("frac_misordered_eigenvalues", 0) > 0.01:
        flags.append("eigenvalues not sorted L1>=L2>=L3 "
                     f"({row['frac_misordered_eigenvalues']:.1%})")

    # ---- expected error at this SNR
    if snr:
        row["snr"] = float(snr)
        for tag in ("WM", "GM"):
            r = rel.expected_error(snr, n_dirs, b, tag, method)
            for k in rel.METRICS:
                row[f"exp_{tag}_{k}_bias"] = r[k]["bias"]
                row[f"exp_{tag}_{k}_sd"] = r[k]["sd"]
            row[f"exp_{tag}_invalid_pct"] = r["invalid_pct"]
        if row["exp_GM_FA_bias"] > WARN_EXPECTED_FA_BIAS:
            flags.append(f"SNR {snr:g}: expected GM-like FA bias "
                         f"{row['exp_GM_FA_bias']:+.3f} (> "
                         f"{WARN_EXPECTED_FA_BIAS})")
        if row["exp_GM_MODE_sd"] > WARN_EXPECTED_MODE_SD:
            flags.append(f"SNR {snr:g}: tensor mode in low-FA tissue is "
                         f"not interpretable (voxel SD "
                         f"{row['exp_GM_MODE_sd']:.2f})")
    row["status"] = status if not flags or status == "FAIL" else "WARN"
    row["flags"] = flags
    return row


# ------------------------------------------------------------ group outliers
OUTLIER_METRICS = ("fa_mean", "fa_p95", "md_median", "frac_nonphysical",
                   "frac_fa_saturated", "frac_fa_gt_0p4", "nde_median")


def robust_z(values: np.ndarray) -> np.ndarray:
    """(x - median) / scale, scale = max(1.4826*MAD, 2 % of |median|, 0.005).

    The floor stops near-identical groups (MAD ~ 0) from turning tiny
    differences into huge z-scores.
    """
    v = np.asarray(values, dtype=float)
    med = np.nanmedian(v)
    mad = np.nanmedian(np.abs(v - med)) * 1.4826
    scale = max(mad, 0.02 * abs(med), 0.005) if np.isfinite(med) else np.nan
    if not np.isfinite(scale):
        return np.zeros_like(v)
    return (v - med) / scale


def flag_outliers(rows: List[dict], z_thresh: float = OUTLIER_Z) -> None:
    """Add ``z_<metric>`` and outlier flags in place. Needs
    ``MIN_GROUP_FOR_OUTLIERS`` subjects sharing the same workflow label
    prefix (compared within workflow)."""
    by_wf: Dict[str, List[dict]] = {}
    for r in rows:
        by_wf.setdefault(r.get("workflow", ""), []).append(r)
    for grp in by_wf.values():
        if len(grp) < MIN_GROUP_FOR_OUTLIERS:
            continue
        for m in OUTLIER_METRICS:
            vals = np.array([g.get(m, np.nan) for g in grp], dtype=float)
            if np.isfinite(vals).sum() < MIN_GROUP_FOR_OUTLIERS:
                continue
            z = robust_z(vals)
            for g, zi in zip(grp, z):
                g[f"z_{m}"] = float(zi)
                if abs(zi) > z_thresh:
                    g["flags"].append(f"group outlier: {m} z={zi:+.1f}")
                    if g["status"] == "PASS":
                        g["status"] = "WARN"


# ------------------------------------------------------------------- output
def _workflow(label: str) -> str:
    return label.split("/", 1)[0] if "/" in label else ""


def write_tsv(rows: List[dict], out: Path) -> Path:
    cols: List[str] = []
    for r in rows:
        for k in r:
            if k not in cols and k != "flags":
                cols.append(k)
    cols.append("flags")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            cells = []
            for c in cols:
                v = r.get(c, "")
                if c == "flags":
                    v = "; ".join(v)
                elif isinstance(v, float):
                    v = "" if math.isnan(v) else f"{v:.6g}"
                cells.append(str(v))
            f.write("\t".join(cells) + "\n")
    return out


def write_html(rows: List[dict], out: Path) -> Path:
    colour = {"PASS": "#d9f2d9", "WARN": "#fff2cc", "FAIL": "#f8d0d0"}
    show = ["status", "n_voxels", "snr", "frac_nonphysical", "fa_mean",
            "md_median", "frac_fa_saturated", "nde_median"]
    th = "".join(f"<th>{html.escape(c)}</th>" for c in ["subject"] + show
                 + ["flags"])
    body = []
    for r in rows:
        tds = [f"<td>{html.escape(r['label'])}</td>"]
        for c in show:
            v = r.get(c, "")
            if isinstance(v, float):
                v = "" if math.isnan(v) else f"{v:.4g}"
            style = f' style="background:{colour[r["status"]]}"' \
                if c == "status" else ""
            tds.append(f"<td{style}>{html.escape(str(v))}</td>")
        tds.append("<td>" + html.escape("; ".join(r["flags"])) + "</td>")
        body.append("<tr>" + "".join(tds) + "</tr>")
    n = {s: sum(r["status"] == s for r in rows)
         for s in ("PASS", "WARN", "FAIL")}
    doc = (
        "<!doctype html><meta charset='utf-8'><title>NERVES QC</title>"
        "<style>body{font:14px system-ui;margin:2rem}"
        "table{border-collapse:collapse}td,th{border:1px solid #ccc;"
        "padding:4px 8px;text-align:left}th{background:#f3f3f3}</style>"
        f"<h1>NERVES tensor-fit QC</h1><p>nerves_qc {__version__} — "
        f"{len(rows)} set(s): {n['PASS']} pass, {n['WARN']} warn, "
        f"{n['FAIL']} fail.</p><table><tr>{th}</tr>{''.join(body)}</table>"
        "<p>Thresholds and interpretation: docs/reliability.md.</p>")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(doc)
    return out


def run(root: Path, out: Path, snr: Optional[float] = None,
        n_dirs: int = 30, b: float = 1000.0, method: str = "ols") -> List[dict]:
    root = Path(root)
    sets = discover(root)
    rows: List[dict] = []
    for label, paths in sets.items():
        report = _find_report(paths, root)
        s = snr
        if s is None and report:
            s = (report.get("snr") or {}).get("snr_median")
        try:
            row = analyse_subject(paths, s, n_dirs, b, method)
        except Exception as exc:                     # keep going per subject
            row = {"status": "FAIL", "flags": [f"could not analyse: {exc}"],
                   "n_voxels": 0}
        row["label"] = label
        row["workflow"] = _workflow(label)
        rows.append(row)
    flag_outliers(rows)
    out = Path(out)
    write_tsv(rows, out / "nerves_qc_summary.tsv")
    (out / "nerves_qc_summary.json").write_text(
        json.dumps(rows, indent=2, default=float))
    write_html(rows, out / "nerves_qc_report.html")
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("root", type=Path, help="Derivatives directory")
    ap.add_argument("--out", type=Path, default=None,
                    help="Output directory (default ROOT/nerves_qc)")
    ap.add_argument("--snr", type=float, default=None,
                    help="b=0 SNR to use for all subjects (otherwise taken "
                         "from each NERVES report if present)")
    ap.add_argument("--dirs", type=int, default=30,
                    help="Number of diffusion directions (default 30)")
    ap.add_argument("--b", type=float, default=1000.0)
    ap.add_argument("--method", choices=["ols", "wls"], default="ols",
                    help="Tensor fit used to compute expected error")
    a = ap.parse_args(argv)

    rows = run(a.root, a.out or (a.root / "nerves_qc"), a.snr, a.dirs, a.b,
               a.method)
    if not rows:
        print("No tensor scalar maps (FA+MD) found.", file=sys.stderr)
        return 1
    for r in rows:
        extra = f"  [{'; '.join(r['flags'])}]" if r["flags"] else ""
        print(f"{r['status']:4s} {r['label']}{extra}")
    return 3 if any(r["status"] == "FAIL" for r in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
