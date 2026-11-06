#!/usr/bin/env python3
"""
nerves_bids.py
==============
Run NERVES over a BIDS dataset.

For every ``sub-*/[ses-*/]dwi/*_dwi.nii[.gz]`` with its ``.bval``/``.bvec``
this finds the participant's T1w, reads ``PhaseEncodingDirection`` and
``TotalReadoutTime`` from the JSON sidecar, picks an opposite-direction
``fmap/*_epi`` image for topup if one exists, and runs ``nde_pipeline.py``
with outputs in ``<bids_dir>/derivatives/nerves/<sub>[_<ses>][_<entities>]``.
Finally it writes ``group_summary.csv`` (see ``nerves_group_summary.py``).

    python nerves_bids.py /data/bids [--out DIR] [--participant-label 01 02]
                          [--dry_run] [-- extra nde_pipeline.py options]

Unrecognised options (e.g. ``--wls --max_b 1500 --n_threads 8 --resume``)
are passed straight to ``nde_pipeline.py``.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Optional

HERE = Path(__file__).resolve().parent


def opposite(pe: str) -> str:
    """Opposite BIDS phase-encode direction ('j-' <-> 'j')."""
    return pe[:-1] if pe.endswith("-") else pe + "-"


def _read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {}


def _strip_nii(p: Path) -> str:
    n = p.name
    for ext in (".nii.gz", ".nii"):
        if n.endswith(ext):
            return n[: -len(ext)]
    return n


def find_t1w(sub_dir: Path, ses_dir: Optional[Path]) -> Optional[Path]:
    for base in (ses_dir, sub_dir):
        if base is None:
            continue
        hits = sorted((base / "anat").glob("*_T1w.nii*"))
        if hits:
            return hits[0]
    return None


def find_rpe(base: Path, pe: str) -> Optional[Path]:
    """First fmap EPI image whose PE direction is opposite to ``pe``."""
    for img in sorted((base / "fmap").glob("*_epi.nii*")):
        meta = _read_json(img.parent / (_strip_nii(img) + ".json"))
        if meta.get("PhaseEncodingDirection") == opposite(pe):
            return img
    return None


def discover(bids_dir: Path, participants: Optional[list] = None) -> list:
    """List of job dicts (one per DWI series)."""
    bids_dir = Path(bids_dir)
    jobs = []
    for sub_dir in sorted(bids_dir.glob("sub-*")):
        if not sub_dir.is_dir():
            continue
        label = sub_dir.name[4:]
        if participants and label not in participants:
            continue
        ses_dirs = sorted(d for d in sub_dir.glob("ses-*") if d.is_dir())
        for ses_dir in (ses_dirs or [None]):
            base = ses_dir or sub_dir
            for dwi in sorted((base / "dwi").glob("*_dwi.nii*")):
                stem = _strip_nii(dwi)
                bval = dwi.parent / (stem + ".bval")
                bvec = dwi.parent / (stem + ".bvec")
                if not (bval.exists() and bvec.exists()):
                    continue
                meta = _read_json(dwi.parent / (stem + ".json"))
                pe = meta.get("PhaseEncodingDirection")
                jobs.append({
                    "name": stem[: -len("_dwi")] if stem.endswith("_dwi")
                    else stem,
                    "dwi": dwi, "bval": bval, "bvec": bvec,
                    "t1": find_t1w(sub_dir, ses_dir),
                    "pe_dir": pe,
                    "readout_time": meta.get("TotalReadoutTime"),
                    "rpe_b0": find_rpe(base, pe) if pe else None,
                })
    return jobs


def build_command(job: dict, out_root: Path, extra: Optional[list] = None,
                  python: str = sys.executable) -> list:
    cmd = [python, str(HERE / "nde_pipeline.py"),
           "--dwi", str(job["dwi"]), "--bval", str(job["bval"]),
           "--bvec", str(job["bvec"]), "--out", str(Path(out_root) / job["name"])]
    if job["t1"] is not None:
        cmd += ["--t1", str(job["t1"])]
    if job["pe_dir"]:
        cmd += ["--pe_dir", job["pe_dir"]]
    if job["readout_time"]:
        cmd += ["--readout_time", str(job["readout_time"])]
    if job["rpe_b0"] is not None:
        cmd += ["--rpe_b0", str(job["rpe_b0"])]
    return cmd + list(extra or [])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("bids_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--participant-label", nargs="*", default=None)
    ap.add_argument("--dry_run", action="store_true",
                    help="Print the commands without running them")
    a, extra = ap.parse_known_args(argv)
    extra = [e for e in extra if e != "--"]

    out_root = a.out or (a.bids_dir / "derivatives" / "nerves")
    jobs = discover(a.bids_dir, a.participant_label)
    if not jobs:
        print("No DWI series with bval/bvec found.", file=sys.stderr)
        return 1
    results = []
    for job in jobs:
        if job["t1"] is None:
            print(f"[{job['name']}] no T1w found: the full pipeline needs a "
                  "T1, skipping this series.", file=sys.stderr)
            results.append((job["name"], "skipped (no T1w)"))
            continue
        if job["readout_time"] is None:
            print(f"[{job['name']}] WARNING: no TotalReadoutTime in sidecar; "
                  "the 0.05 s placeholder will be used.", file=sys.stderr)
        cmd = build_command(job, out_root, extra)
        print(" ".join(cmd))
        if a.dry_run:
            results.append((job["name"], "dry run"))
            continue
        rc = subprocess.run(cmd).returncode
        results.append((job["name"], f"exit {rc}"))
    print("\nSummary:")
    for name, status in results:
        print(f"  {name}: {status}")
    if not a.dry_run:
        try:
            import nerves_group_summary as gs
            gs.write_csv(gs.collect(out_root), out_root / "group_summary.csv")
        except Exception as e:     # never fail the run on summary problems
            print("Group summary not written:", e, file=sys.stderr)
    bad = [s for _, s in results if s.startswith("exit") and s not in
           ("exit 0", "exit 2")]
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
