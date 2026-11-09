#!/usr/bin/env python3
"""
nerves_group_summary.py
=======================
Collect the ``nde_pipeline_report.json`` of every NERVES run under a
directory into one CSV (one row per run) for group-level inspection and
outlier screening.

    python nerves_group_summary.py DERIVATIVES_DIR [--out group_summary.csv]

Columns: run (relative directory), run_mode, NDE statistics, QC results,
shape/SNR statistics when present, software versions, timestamp.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def _flat(prefix: str, d: dict) -> dict:
    return {f"{prefix}{k}": v for k, v in d.items()
            if isinstance(v, (int, float, bool, str)) or v is None}


def collect(root: Path) -> list:
    """One flat dict per report found under ``root``."""
    root = Path(root)
    rows = []
    for rpt in sorted(root.rglob("nde_pipeline_report.json")):
        try:
            data = json.loads(rpt.read_text())
        except (OSError, ValueError):
            continue
        row = {"run": str(rpt.parent.relative_to(root)) or ".",
               "timestamp": data.get("timestamp"),
               "run_mode": data.get("run_mode")}
        row.update(_flat("nde_", {k.replace("nde_", "", 1)
                                  if k.startswith("nde_") else k: v
                                  for k, v in
                                  data.get("nde_statistics", {}).items()}))
        row.update(_flat("qc_", data.get("qc_results", {})))
        row.update(_flat("shape_", data.get("shape_statistics", {})))
        row.update(_flat("", data.get("snr", {})))
        v = data.get("versions", {})
        row["nerves_version"] = v.get("nerves")
        row["git_commit"] = v.get("git_commit")
        row["fsl_version"] = v.get("fsl")
        rows.append(row)
    return rows


def write_csv(rows: list, out: Path) -> Path:
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    rows = collect(a.root)
    out = write_csv(rows, a.out or (a.root / "group_summary.csv"))
    print(f"{len(rows)} run(s) -> {out}")
    return 0 if rows else 1


if __name__ == "__main__":
    raise SystemExit(main())
