#!/usr/bin/env python3
"""
nerves_reliability.py
=====================
How much can you trust tensor-derived scalars (FA, MD, tensor mode, NDE) at a
given SNR, number of directions and b-value?

Monte-Carlo (Rician noise, same tensor fit as the pipeline) gives the
expected **bias** and **voxelwise SD** of each metric for a tissue-like
tensor, plus the fraction of voxels whose fit yields non-physical
(non-positive) eigenvalues. Used three ways:

* as a protocol-planning calculator (CLI below),
* by ``nerves_qc.py`` to attach "expected error at this subject's SNR" to
  the QC report,
* as a library (``expected_error``, ``reliability_table``,
  ``min_snr_for_bias``).

CLI examples::

    python nerves_reliability.py --snr 20 --dirs 30 --b 1000
    python nerves_reliability.py --sweep-snr 10 15 20 30 40 --dirs 30
    python nerves_reliability.py --required FA:0.02 --tissue GM --dirs 30

Caveats (also in docs/reliability.md): single-shell b-value, Rician noise
model with one tensor orientation distribution (uniform random), a single
tissue tensor per call. It does not model motion, eddy residuals, partial
volume, or crossing fibres.
"""

from __future__ import annotations

import argparse
import functools
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Union

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validation import noise_simulation as ns  # noqa: E402

METRICS = ("FA", "MD", "MODE", "NDE")

# Representative principal diffusivities (1e-3 mm^2/s) for 3 T adult brain.
TISSUES: Dict[str, Sequence[float]] = {
    "WM": (1.7, 0.35, 0.30),          # coherent white matter, FA ~0.78
    "WM_LOWFA": (1.2, 0.55, 0.40),    # crossing / low-coherence WM, FA ~0.53
    "GM": (0.95, 0.75, 0.65),         # cortical / deep grey, FA ~0.19
}


def _resolve_tissue(tissue: Union[str, Sequence[float]]):
    if isinstance(tissue, str):
        if tissue not in TISSUES:
            raise ValueError(f"Unknown tissue '{tissue}'. "
                             f"Choose from {sorted(TISSUES)} or give three "
                             "principal diffusivities (1e-3 mm^2/s).")
        return tuple(float(x) for x in TISSUES[tissue])
    lam = tuple(float(x) for x in tissue)
    if len(lam) != 3 or min(lam) <= 0:
        raise ValueError("tissue must be 3 positive diffusivities")
    return lam


@functools.lru_cache(maxsize=256)
def _expected_error_cached(snr: float, n_dirs: int, b: float, lam: tuple,
                           method: str, n: int, seed: int, n_b0: int):
    rng = np.random.default_rng(seed)
    bvals, bvecs = ns.make_scheme(n_dirs=n_dirs, n_b0=n_b0, b=b)
    sim = ns.simulate_metrics(lam, snr, n, bvals, bvecs, rng, method)
    truth = ns.true_metrics(lam)
    out = {"invalid_pct": float(100.0 * sim["frac_invalid"]), "truth": truth}
    for k in METRICS:
        v = sim[k][np.isfinite(sim[k])]
        if v.size < 10:
            out[k] = {"truth": truth[k], "bias": float("nan"),
                      "sd": float("nan"), "rmse": float("nan")}
            continue
        bias = float(v.mean() - truth[k])
        sd = float(v.std())
        out[k] = {"truth": truth[k], "bias": bias, "sd": sd,
                  "rmse": float(np.sqrt(bias ** 2 + sd ** 2))}
    return out


def expected_error(snr: float, n_dirs: int = 30, b: float = 1000.0,
                   tissue: Union[str, Sequence[float]] = "WM",
                   method: str = "ols", n: int = 4000, seed: int = 0,
                   n_b0: int = 3) -> dict:
    """Expected bias / SD / RMSE of FA, MD (1e-3 mm^2/s), MODE and NDE.

    Returns ``{"invalid_pct", "truth", "FA": {"truth","bias","sd","rmse"}, ...}``.
    Deterministic for a given seed; results are cached.
    """
    if not np.isfinite(snr) or snr <= 0:
        raise ValueError("snr must be a positive number")
    if n_dirs < 6:
        raise ValueError("n_dirs must be >= 6 for a tensor fit")
    if method not in ("ols", "wls"):
        raise ValueError("method must be 'ols' or 'wls'")
    return _expected_error_cached(float(snr), int(n_dirs), float(b),
                                  _resolve_tissue(tissue), method, int(n),
                                  int(seed), int(n_b0))


def reliability_table(snrs: Iterable[float], n_dirs: int = 30,
                      b: float = 1000.0,
                      tissue: Union[str, Sequence[float]] = "WM",
                      method: str = "ols", n: int = 4000,
                      seed: int = 0) -> List[dict]:
    """One flat row per SNR (bias and SD per metric, invalid %)."""
    rows = []
    for s in snrs:
        r = expected_error(s, n_dirs, b, tissue, method, n, seed)
        row = {"snr": float(s), "invalid_pct": r["invalid_pct"]}
        for k in METRICS:
            row[f"{k}_bias"] = r[k]["bias"]
            row[f"{k}_sd"] = r[k]["sd"]
        rows.append(row)
    return rows


def min_snr_for_bias(metric: str, max_abs_bias: float, n_dirs: int = 30,
                     b: float = 1000.0,
                     tissue: Union[str, Sequence[float]] = "WM",
                     method: str = "ols", n: int = 4000, seed: int = 0,
                     grid: Sequence[float] = (5, 7.5, 10, 12.5, 15, 20, 25,
                                              30, 40, 50, 60, 80, 100, 150)
                     ) -> Optional[float]:
    """Smallest grid SNR at which |bias| and SD of ``metric`` are both
    <= ``max_abs_bias`` (None if no grid SNR qualifies).

    Requiring SD as well as bias makes the answer "the metric is within
    +-``max_abs_bias`` of the truth in a typical voxel".
    """
    metric = metric.upper()
    if metric not in METRICS:
        raise ValueError(f"metric must be one of {METRICS}")
    for s in grid:
        r = expected_error(s, n_dirs, b, tissue, method, n, seed)[metric]
        if abs(r["bias"]) <= max_abs_bias and r["sd"] <= max_abs_bias:
            return float(s)
    return None


def roi_sd(voxel_sd: float, n_voxels: int) -> float:
    """SD of an ROI mean of independent voxels (optimistic: ignores
    spatial correlation and does not shrink bias)."""
    return float(voxel_sd / np.sqrt(max(int(n_voxels), 1)))


# --------------------------------------------------------------------- CLI
def _fmt_table(rows: List[dict]) -> str:
    head = ["SNR", "invalid %"] + [f"{k} bias±sd" for k in METRICS]
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for r in rows:
        cells = [f"{r['snr']:g}", f"{r['invalid_pct']:.1f}"]
        for k in METRICS:
            cells.append(f"{r[k + '_bias']:+.3f}±{r[k + '_sd']:.3f}")
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--snr", type=float, help="Single SNR (b=0 SNR) to report")
    ap.add_argument("--sweep-snr", type=float, nargs="+",
                    help="Report a table over these SNRs")
    ap.add_argument("--required", metavar="METRIC:TOL",
                    help="Minimum SNR so |bias| and SD of METRIC <= TOL "
                         "(e.g. FA:0.02)")
    ap.add_argument("--dirs", type=int, default=30,
                    help="Number of diffusion directions (default 30)")
    ap.add_argument("--b", type=float, default=1000.0)
    ap.add_argument("--tissue", default="WM",
                    help=f"One of {sorted(TISSUES)} or 'l1,l2,l3' in "
                         "1e-3 mm^2/s (default WM)")
    ap.add_argument("--method", choices=["ols", "wls"], default="ols")
    ap.add_argument("--n", type=int, default=4000, help="Monte-Carlo voxels")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)

    tissue: Union[str, Sequence[float]] = a.tissue
    try:
        if "," in a.tissue:
            tissue = [float(x) for x in a.tissue.split(",")]
        lam = _resolve_tissue(tissue)
    except ValueError as exc:
        ap.error(f"--tissue: {exc}")
    if a.dirs < 6:
        ap.error("--dirs must be >= 6")
    if a.b <= 0:
        ap.error("--b must be positive")
    if a.snr is not None and a.snr <= 0 or \
            any(x <= 0 for x in (a.sweep_snr or [])):
        ap.error("SNR values must be positive")
    t = ns.true_metrics(lam)
    print(f"Tissue {a.tissue} λ={lam} → truth FA={t['FA']:.3f} "
          f"MD={t['MD']:.3f} MODE={t['MODE']:.3f} NDE={t['NDE']:.3f}; "
          f"{a.dirs} dirs, b={a.b:g}, {a.method.upper()}\n")

    if a.required:
        metric, _, tol = a.required.partition(":")
        try:
            tol_f = float(tol)
            if tol_f <= 0:
                raise ValueError
            metric = metric.upper()
            if metric not in METRICS:
                raise ValueError
        except ValueError:
            ap.error("--required must look like METRIC:TOL with METRIC in "
                     f"{'/'.join(METRICS)} and TOL > 0, e.g. FA:0.02")
        s = min_snr_for_bias(metric, tol_f, a.dirs, a.b, tissue,
                             a.method, a.n, a.seed)
        print(f"Minimum SNR for |bias|,SD of {metric} ≤ {tol_f:g}: "
              + (f"{s:g}" if s is not None else "not reached on grid (>150)"))
        return 0 if s is not None else 1

    snrs = a.sweep_snr or ([a.snr] if a.snr else [10, 15, 20, 30, 40])
    print(_fmt_table(reliability_table(snrs, a.dirs, a.b, tissue, a.method,
                                       a.n, a.seed)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
