"""
tests/test_nde_math.py
----------------------
Unit tests for the core NDE computation.
No FSL, MRtrix3, or NIfTI files required.
"""

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Inline NDE implementation for testing (mirrors nde_pipeline.py logic)
# ---------------------------------------------------------------------------

def nde_scalar(l1: float, l2: float, l3: float) -> float:
    """Compute NDE for a single voxel with scalar eigenvalues."""
    ev = np.array([l1, l2, l3], dtype=np.float64)
    p = ev / ev.sum()
    log_p = np.where(p > 0.0, np.log(p), 0.0)
    H = -np.dot(p, log_p)
    return float(np.clip(H / np.log(3.0), 0.0, 1.0))


def nde_vectorised(l1: np.ndarray, l2: np.ndarray,
                   l3: np.ndarray) -> np.ndarray:
    """Vectorised NDE computation matching nde_pipeline.py."""
    eigs = np.stack([l1, l2, l3], axis=-1).astype(np.float64)
    eig_sum = eigs.sum(axis=-1, keepdims=True)
    p = eigs / eig_sum
    log_p = np.where(p > 0.0, np.log(p), 0.0)
    H = -np.sum(p * log_p, axis=-1)
    return np.clip(H / np.log(3.0), 0.0, 1.0)


def fa_scalar(l1: float, l2: float, l3: float) -> float:
    """FA for a single voxel."""
    ev = np.array([l1, l2, l3], dtype=np.float64)
    lm = ev.mean()
    return float(np.sqrt(1.5) * np.sqrt(np.sum((ev - lm) ** 2))
                 / np.sqrt(np.sum(ev ** 2)))


# ---------------------------------------------------------------------------
# Boundary conditions
# ---------------------------------------------------------------------------

class TestBoundaryConditions:

    def test_isotropic_nde_equals_one(self):
        """Perfectly isotropic tensor → NDE = 1."""
        for ev in [(1.0, 1.0, 1.0), (0.5, 0.5, 0.5), (2.0, 2.0, 2.0)]:
            assert abs(nde_scalar(*ev) - 1.0) < 1e-10, \
                f"Expected NDE=1.0 for isotropic {ev}, got {nde_scalar(*ev)}"

    def test_nde_range(self):
        """NDE must always be in [0, 1]."""
        rng = np.random.default_rng(0)
        l1 = rng.uniform(0.5, 2.5, 10000)
        l2 = rng.uniform(0.1, 1.0, 10000)
        l3 = rng.uniform(0.05, 0.5, 10000)
        # Enforce l1 >= l2 >= l3
        ev = np.sort(np.stack([l1, l2, l3], axis=-1), axis=-1)[:, ::-1]
        result = nde_vectorised(ev[:, 0], ev[:, 1], ev[:, 2])
        assert result.min() >= 0.0
        assert result.max() <= 1.0

    def test_scale_invariance(self):
        """NDE must be scale-invariant: multiplying all λ by constant is no-op."""
        base = (1.5, 0.75, 0.25)
        nde_base = nde_scalar(*base)
        for scale in [0.1, 0.5, 2.0, 10.0, 1000.0]:
            scaled = tuple(v * scale for v in base)
            assert abs(nde_scalar(*scaled) - nde_base) < 1e-10, \
                f"Scale invariance failed at scale={scale}"

    def test_near_zero_minor_eigenvalue(self):
        """Very small minor eigenvalues should not produce NaN or error."""
        result = nde_scalar(1.5, 1e-10, 1e-10)
        assert np.isfinite(result)
        assert 0.0 <= result <= 1.0

    def test_symmetry_of_minor_eigenvalues(self):
        """Swapping λ₂ and λ₃ should not change NDE (entropy is symmetric)."""
        assert abs(nde_scalar(1.5, 0.75, 0.25) -
                   nde_scalar(1.5, 0.25, 0.75)) < 1e-10


# ---------------------------------------------------------------------------
# Manuscript Table 1 values (verified against corrected manuscript)
# ---------------------------------------------------------------------------

class TestManuscriptTable1:

    TABLE1 = [
        ("A", 1.5, 0.50, 0.50, 0.8650),
        ("B", 1.5, 0.75, 0.25, 0.8173),
        ("C", 1.2, 0.60, 0.60, 0.9464),
        ("D", 1.2, 0.80, 0.40, 0.9207),
    ]

    @pytest.mark.parametrize("name,l1,l2,l3,expected", TABLE1)
    def test_nde_value(self, name, l1, l2, l3, expected):
        """NDE must match manuscript Table 1 to 3 decimal places."""
        got = nde_scalar(l1, l2, l3)
        assert abs(got - expected) < 0.001, \
            f"Config {name}: expected {expected:.4f}, got {got:.4f}"

    def test_ab_absolute_difference(self):
        """Corrected manuscript: |ΔNDE| > |ΔFA| between Configs A and B."""
        nde_A = nde_scalar(1.5, 0.50, 0.50)
        nde_B = nde_scalar(1.5, 0.75, 0.25)
        fa_A  = fa_scalar(1.5, 0.50, 0.50)
        fa_B  = fa_scalar(1.5, 0.75, 0.25)
        assert abs(nde_A - nde_B) > abs(fa_A - fa_B), \
            "NDE absolute difference should exceed FA absolute difference"

    def test_prolate_higher_nde_than_asymmetric(self):
        """Prolate configs (equal minor eigenvalues) have higher NDE."""
        assert nde_scalar(1.5, 0.50, 0.50) > nde_scalar(1.5, 0.75, 0.25)
        assert nde_scalar(1.2, 0.60, 0.60) > nde_scalar(1.2, 0.80, 0.40)


# ---------------------------------------------------------------------------
# Gradient property
# ---------------------------------------------------------------------------

class TestGradientProperty:

    def test_gradient_steepest_near_zero(self):
        """d/dp[-p ln p] should be largest in magnitude near p=0."""
        # Gradient at p=0.01 should far exceed gradient at p=0.9
        p_small = 0.01
        p_large = 0.90
        grad_small = abs(-np.log(p_small) - 1)
        grad_large = abs(-np.log(p_large) - 1)
        assert grad_small > grad_large * 3, \
            "Gradient near p=0 should be much steeper than near p=1"

    def test_gradient_zero_at_one_over_e(self):
        """Gradient of -p ln p is zero at p = 1/e."""
        p = 1.0 / np.e
        grad = -np.log(p) - 1
        assert abs(grad) < 1e-10


# ---------------------------------------------------------------------------
# Vectorised implementation
# ---------------------------------------------------------------------------

class TestVectorised:

    def test_matches_scalar(self):
        """Vectorised and scalar implementations must agree."""
        rng = np.random.default_rng(42)
        n = 1000
        ev = rng.uniform(0.1, 2.0, (n, 3))
        ev = np.sort(ev, axis=-1)[:, ::-1]

        scalar_results = np.array([
            nde_scalar(ev[i, 0], ev[i, 1], ev[i, 2]) for i in range(n)
        ])
        vec_results = nde_vectorised(ev[:, 0], ev[:, 1], ev[:, 2])

        np.testing.assert_allclose(scalar_results, vec_results, atol=1e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
