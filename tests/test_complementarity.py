"""
tests/test_complementarity.py
-----------------------------
Quantifies how much information NDE carries beyond FA and tensor mode.

NDE is a function of the three (normalised) eigenvalues, i.e. of two free
shape parameters. FA and tensor mode (Ennis & Kindlmann 2006) are two
orthogonal shape invariants, so (FA, mode) determine NDE completely. These
tests pin that down so documentation cannot over-claim "complementarity".
"""

import numpy as np

RNG = np.random.default_rng(1)


def _sample(n=100_000):
    l1 = RNG.uniform(1.0, 1.9, n)
    l2 = RNG.uniform(0.2, 1.0, n)
    l3 = RNG.uniform(0.1, 1.0, n)
    return np.sort(np.stack([l1, l2, l3], 1), 1)[:, ::-1]


def _invariants(L):
    p = L / L.sum(1, keepdims=True)
    nde = -(p * np.log(p)).sum(1) / np.log(3)
    d = L - L.mean(1, keepdims=True)
    fa = np.sqrt(1.5 * (d ** 2).sum(1) / (L ** 2).sum(1))
    dn = d / np.sqrt((d ** 2).sum(1, keepdims=True))
    mode = 3 * np.sqrt(6) * np.prod(dn, axis=1)   # = 3*sqrt(6)*det of diag(dn)
    return nde, fa, mode


def _poly(a, b, deg=5):
    return np.stack([a ** i * b ** j for i in range(deg + 1)
                     for j in range(deg + 1 - i)], 1)


def _r2(y, X):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    return 1 - (y - X @ beta).var() / y.var()


def test_nde_is_nearly_a_function_of_fa_and_mode():
    nde, fa, mode = _invariants(_sample())
    assert _r2(nde, _poly(fa, mode)) > 0.999


def test_nde_mostly_explained_by_fa_alone():
    nde, fa, _ = _invariants(_sample())
    assert _r2(nde, _poly(fa, np.zeros_like(fa))) > 0.95


def test_matched_fa_different_lambda2_lambda3_gives_different_nde():
    """The legitimate (and only) added contrast: at fixed FA, NDE varies
    with the shape of the minor-eigenvalue split (i.e. with tensor mode)."""
    from scipy.optimize import brentq

    def fa_of(L):
        d = np.array(L) - np.mean(L)
        return np.sqrt(1.5 * (d ** 2).sum() / (np.array(L) ** 2).sum())

    def nde_of(L):
        p = np.array(L) / sum(L)
        return -(p * np.log(p)).sum() / np.log(3)

    # Fix trace = 2.4 and FA = 0.5; vary split of the minor eigenvalues.
    target = 0.5
    results = []
    for l2_frac in (0.5, 0.7):       # l2 as a fraction of (l2 + l3)
        def f(l1):
            rest = 2.4 - l1
            return fa_of([l1, rest * l2_frac, rest * (1 - l2_frac)]) - target
        l1 = brentq(f, 0.81, 2.39)
        rest = 2.4 - l1
        L = [l1, rest * l2_frac, rest * (1 - l2_frac)]
        results.append((fa_of(L), nde_of(L)))
    (fa_a, nde_a), (fa_b, nde_b) = results
    assert abs(fa_a - fa_b) < 1e-6
    assert abs(nde_a - nde_b) > 0.005
