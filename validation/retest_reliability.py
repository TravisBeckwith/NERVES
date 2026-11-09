#!/usr/bin/env python3
"""
validation/retest_reliability.py
-------------------------
Test-retest reliability of NERVES metrics (NDE vs FA vs tensor mode ...).

Input: a CSV manifest, one row per subject x session, with columns
    subject, session, mask, <metric columns...>
where each metric column is the path to that metric's NIfTI for that scan
(e.g. nde, fa, mode, md) and ``mask`` is the ROI (e.g. a WM mask in the same
space). Exactly two sessions per subject are used.

For each metric the ROI mean (NaN-aware) is taken per scan, then:
    ICC(3,1)  consistency,  ICC(2,1)  absolute agreement,
    wCV %     within-subject coefficient of variation,
    between-subject SD.
A metric with higher ICC and lower wCV is more reproducible. Use it to check
whether NDE is actually more reliable than FA/mode before relying on it.

    python validation/retest_reliability.py manifest.csv [--out table.csv]
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import nibabel as nib
import numpy as np


def icc_two_sessions(x: np.ndarray) -> dict:
    """x: (n_subjects, 2) -> ICC(2,1), ICC(3,1), wCV% and between-subject SD."""
    x = np.asarray(x, dtype=float)
    n, k = x.shape
    if n < 3 or k != 2:
        raise ValueError("Need >= 3 subjects with exactly 2 sessions each.")
    grand = x.mean()
    msr = k * ((x.mean(1) - grand) ** 2).sum() / (n - 1)
    msc = n * ((x.mean(0) - grand) ** 2).sum() / (k - 1)
    sse = ((x - x.mean(1, keepdims=True) - x.mean(0, keepdims=True)
            + grand) ** 2).sum()
    mse = sse / ((n - 1) * (k - 1))
    icc3 = (msr - mse) / (msr + (k - 1) * mse)
    icc2 = (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
    within_var = ((x[:, 0] - x[:, 1]) ** 2 / 2).mean()
    wcv = 100 * np.sqrt(within_var) / abs(grand) if grand else float("nan")
    return {"ICC2_1": float(icc2), "ICC3_1": float(icc3), "wCV_pct": float(wcv),
            "between_sd": float(x.mean(1).std(ddof=1)), "n": n}


def roi_mean(metric_path: Path, mask_path: Path) -> float:
    m = np.asarray(nib.load(str(mask_path)).dataobj) > 0
    v = np.asarray(nib.load(str(metric_path)).dataobj, dtype=float)[m]
    v = v[np.isfinite(v)]
    return float(v.mean())


def analyse(manifest: Path) -> dict:
    rows = list(csv.DictReader(open(manifest, newline="")))
    metrics = [c for c in rows[0] if c not in ("subject", "session", "mask")]
    by_sub = defaultdict(list)
    base = Path(manifest).parent
    for r in rows:
        by_sub[r["subject"]].append(r)
    out = {}
    for met in metrics:
        pairs = []
        for sub, rs in by_sub.items():
            if len(rs) != 2:
                continue
            rs = sorted(rs, key=lambda r: r["session"])
            pairs.append([roi_mean(base / r[met], base / r["mask"])
                          for r in rs])
        out[met] = icc_two_sessions(np.array(pairs))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("manifest", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    res = analyse(a.manifest)
    print("| metric | n | ICC(3,1) | ICC(2,1) | wCV % | between-subj SD |")
    print("|---|---|---|---|---|---|")
    for m, r in res.items():
        print(f"| {m} | {r['n']} | {r['ICC3_1']:.3f} | {r['ICC2_1']:.3f} | "
              f"{r['wCV_pct']:.2f} | {r['between_sd']:.4f} |")
    if a.out:
        with open(a.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["metric", "n", "ICC3_1", "ICC2_1", "wCV_pct",
                        "between_sd"])
            for m, r in res.items():
                w.writerow([m, r["n"], r["ICC3_1"], r["ICC2_1"],
                            r["wCV_pct"], r["between_sd"]])
    return res


if __name__ == "__main__":
    main()
