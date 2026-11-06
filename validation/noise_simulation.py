#!/usr/bin/env python3
"""
validation/noise_simulation.py
------------------------------
Monte-Carlo comparison of NDE with FA, MD and tensor mode under realistic
acquisition noise. Answers two questions with numbers instead of claims:

1. Precision/bias: how do the metrics behave versus SNR?
2. Discriminability: for three kinds of tissue change, how well does each
   metric separate "before" from "after" voxels (d' = |mean diff| / pooled SD)?

Run:   python validation/noise_simulation.py [--snr 25] [--n 20000]
Uses the same tensor fit as the pipeline (nerves_dti.fit_tensor_signals).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import nerves_dti as nd  # noqa: E402

METRICS = ("FA", "MD", "NDE", "MODE")


def fibonacci_directions(n: int) -> np.ndarray:
    """~Uniform unit vectors on the sphere."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi), np.cos(phi)])


def make_scheme(n_dirs: int = 30, n_b0: int = 3, b: float = 1000.0):
    bvals = np.r_[np.zeros(n_b0), np.full(n_dirs, b)]
    bvecs = np.r_[np.zeros((n_b0, 3)), fibonacci_directions(n_dirs)]
    return bvals, bvecs


def random_rotations(n: int, rng) -> np.ndarray:
    q, r = np.linalg.qr(rng.standard_normal((n, 3, 3)))
    return q * np.sign(np.diagonal(r, axis1=1, axis2=2))[:, None, :]


def simulate_signals(lam_1e3, snr: float, n: int, bvals, bvecs, rng):
    """Rician-noise DWI signals for ``n`` randomly oriented tensors."""
    lam = np.asarray(lam_1e3, dtype=float) * 1e-3
    Q = random_rotations(n, rng)
    D = np.einsum("nij,j,nkj->nik", Q, lam, Q)
    q = np.einsum("mi,nij,mj->nm", bvecs, D, bvecs)
    sig = np.exp(-bvals[None, :] * q)
    s = 1.0 / snr
    return np.sqrt((sig + s * rng.standard_normal(sig.shape)) ** 2
                   + (s * rng.standard_normal(sig.shape)) ** 2)


def simulate_metrics(lam_1e3, snr, n, bvals, bvecs, rng, method="ols"):
    """Fit and return metric arrays (valid voxels only) + invalid fraction."""
    sig = simulate_signals(lam_1e3, snr, n, bvals, bvecs, rng)
    x = nd.fit_tensor_signals(sig, bvals, bvecs, method)
    eig = nd.tensor_eigenvalues(x) * 1e3          # back to 1e-3 mm^2/s
    ok = (eig > 0).all(axis=1)
    eig = eig[ok]
    fa, md = nd.fa_md(eig)
    return {"FA": fa, "MD": md, "NDE": nd.nde_from_eigenvalues(eig),
            "MODE": nd.tensor_mode(eig), "frac_invalid": 1 - ok.mean()}


def true_metrics(lam_1e3):
    e = np.asarray(lam_1e3, dtype=float)[None, :]
    fa, md = nd.fa_md(e)
    return {"FA": float(fa[0]), "MD": float(md[0]),
            "NDE": float(nd.nde_from_eigenvalues(e)[0]),
            "MODE": float(nd.tensor_mode(e)[0])}


def dprime(a, b) -> float:
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    return float(abs(a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2))


def matched_fa_partner(base):
    """Tensor with the same trace and FA as ``base`` but l2 == l3."""
    from scipy.optimize import brentq
    base = np.asarray(base, dtype=float)
    target = float(nd.fa_md(base[None, :])[0][0])
    tr = base.sum()

    def mk(l1):
        rest = tr - l1
        return np.array([l1, rest / 2, rest / 2])

    f = lambda l1: float(nd.fa_md(mk(l1)[None, :])[0][0]) - target  # noqa: E731
    return mk(brentq(f, tr / 3 + 1e-6, tr - 1e-6))


SCENARIOS = {
    "A: same FA & trace, different λ2/λ3 split": ([1.6, 0.45, 0.20], None),
    "B: symmetric radial diffusivity increase": ([1.7, 0.35, 0.30],
                                                 [1.6, 0.45, 0.40]),
    "C: λ2/λ3 redistribution, similar trace": ([1.7, 0.40, 0.25],
                                               [1.7, 0.30, 0.35]),
}


def run_precision(snrs=(10, 15, 25, 40), n=20000, seed=7, method="ols"):
    rng = np.random.default_rng(seed)
    bvals, bvecs = make_scheme()
    truth_lam = [1.7, 0.35, 0.30]
    out = {"truth": true_metrics(truth_lam), "rows": []}
    for snr in snrs:
        r = simulate_metrics(truth_lam, snr, n, bvals, bvecs, rng, method)
        row = {"snr": snr, "invalid_pct": 100 * r["frac_invalid"]}
        for k in ("FA", "NDE", "MODE"):
            v = r[k][np.isfinite(r[k])]
            row[k] = (float(v.mean()), float(v.std()))
        out["rows"].append(row)
    return out


def run_discrimination(snr=25, n=20000, seed=11, method="ols"):
    rng = np.random.default_rng(seed)
    bvals, bvecs = make_scheme()
    res = {}
    for name, (base, alt) in SCENARIOS.items():
        alt = matched_fa_partner(base) if alt is None else alt
        a = simulate_metrics(base, snr, n, bvals, bvecs, rng, method)
        b = simulate_metrics(alt, snr, n, bvals, bvecs, rng, method)
        res[name] = {k: dprime(a[k], b[k]) for k in METRICS}
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--snr", type=float, default=25.0)
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--method", choices=["ols", "wls"], default="ols")
    a = ap.parse_args(argv)

    prec = run_precision(n=a.n, seed=a.seed, method=a.method)
    t = prec["truth"]
    print(f"Truth (1.7/0.35/0.30): FA={t['FA']:.3f} NDE={t['NDE']:.3f} "
          f"MODE={t['MODE']:.3f}\n")
    print("| SNR | invalid % | FA mean±sd | NDE mean±sd | MODE mean±sd |")
    print("|---|---|---|---|---|")
    for r in prec["rows"]:
        cells = " | ".join(f"{r[k][0]:.3f}±{r[k][1]:.3f}"
                           for k in ("FA", "NDE", "MODE"))
        print(f"| {r['snr']} | {r['invalid_pct']:.1f} | {cells} |")
    print(f"\nDiscriminability d' at SNR {a.snr:g} ({a.method.upper()} fit)\n")
    print("| Scenario | " + " | ".join(METRICS) + " |")
    print("|---|" + "---|" * len(METRICS))
    for name, d in run_discrimination(a.snr, a.n, a.seed + 4,
                                      a.method).items():
        print(f"| {name} | " + " | ".join(f"{d[k]:.2f}" for k in METRICS)
              + " |")


if __name__ == "__main__":
    main()
