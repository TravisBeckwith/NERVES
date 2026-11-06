"""Unit tests for nerves_dti (no external tools required)."""

import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import nerves_dti as nd
from validation.make_phantom import make_phantom, TISSUES
from validation.noise_simulation import make_scheme, simulate_metrics


def _save(path, arr, affine=None):
    nib.save(nib.Nifti1Image(np.asarray(arr, dtype=np.float32),
                             np.eye(4) if affine is None else affine), str(path))


# ------------------------------------------------------------------ b-values
class TestGradientTable:

    def test_read_bvecs_both_layouts(self, tmp_path):
        v = np.eye(3)
        np.savetxt(tmp_path / "a.bvec", v)            # 3xN
        np.savetxt(tmp_path / "b.bvec", v.T[:, :3])   # Nx3
        assert nd.read_bvecs(tmp_path / "a.bvec").shape == (3, 3)
        assert nd.read_bvecs(tmp_path / "b.bvec").shape == (3, 3)

    def test_shells_cluster(self):
        sh = nd.shells(np.array([0, 5, 995, 1000, 1005, 2000, 2010]))
        assert len(sh) == 3 and abs(sh[1] - 1000) < 10

    def test_volumes_up_to_keeps_b0(self):
        idx = nd.volumes_up_to(np.array([0, 1000, 2000, 0, 3000]), 1500)
        assert list(idx) == [0, 1, 3]

    def test_range_warnings(self):
        assert nd.check_bvalue_range(np.array([0, 1000] * 5)) == []
        w = nd.check_bvalue_range(np.array([0, 1000, 3000]))
        assert any("kurtosis" in m for m in w) and any("Multi-shell" in m for m in w)

    def _inputs(self, tmp_path, n=10, **kw):
        _save(tmp_path / "d.nii.gz", np.ones((2, 2, 2, n)))
        bvals = kw.get("bvals", np.r_[0, np.full(n - 1, 1000)])
        bvecs = kw.get("bvecs", np.tile([1.0, 0, 0], (n, 1)))
        np.savetxt(tmp_path / "d.bval", np.atleast_2d(bvals))
        np.savetxt(tmp_path / "d.bvec", np.asarray(bvecs).T)
        return tmp_path / "d.nii.gz", tmp_path / "d.bval", tmp_path / "d.bvec"

    def test_validate_ok(self, tmp_path):
        info = nd.validate_dwi_inputs(*self._inputs(tmp_path))
        assert info["n_b0"] == 1 and info["n_dw"] == 9

    def test_validate_length_mismatch(self, tmp_path):
        d, bval, bvec = self._inputs(tmp_path)
        np.savetxt(bval, np.atleast_2d(np.r_[0, np.full(5, 1000)]))
        with pytest.raises(ValueError, match="does not match"):
            nd.validate_dwi_inputs(d, bval, bvec)

    def test_validate_no_b0(self, tmp_path):
        with pytest.raises(ValueError, match="No b=0"):
            nd.validate_dwi_inputs(*self._inputs(tmp_path, bvals=np.full(10, 1000)))

    def test_validate_too_few_directions(self, tmp_path):
        with pytest.raises(ValueError, match=">= 6"):
            nd.validate_dwi_inputs(*self._inputs(
                tmp_path, n=5, bvals=np.r_[0, np.full(4, 1000)]))

    def test_validate_non_unit_bvecs(self, tmp_path):
        bad = np.tile([3.0, 0, 0], (10, 1))
        with pytest.raises(ValueError, match="unit length"):
            nd.validate_dwi_inputs(*self._inputs(tmp_path, bvecs=bad))

    def test_select_shells(self, tmp_path):
        bvals = np.array([0] + [1000] * 6 + [2000] * 3 + [0, 1000.0])
        n = len(bvals)
        data = np.arange(n)[None, None, None, :] * np.ones((2, 2, 2, n))
        _save(tmp_path / "d.nii.gz", data)
        bvecs = np.random.default_rng(0).standard_normal((n, 3))
        np.savetxt(tmp_path / "d.bval", bvals[None])
        np.savetxt(tmp_path / "d.bvec", bvecs.T)
        sel = nd.select_shells(tmp_path / "d.nii.gz", tmp_path / "d.bvec",
                               tmp_path / "d.bval", 1000, tmp_path / "sel")
        out = nib.load(str(sel["dwi"])).get_fdata()
        assert list(out[0, 0, 0]) == [0, 1, 2, 3, 4, 5, 6, 10, 11]
        assert list(nd.read_bvals(sel["bvals"])) == [0] + [1000] * 6 + [0, 1000]
        assert nd.read_bvecs(sel["bvecs"]).shape == (9, 3)

    def test_mean_b0_numpy(self, tmp_path):
        bvals = np.array([0, 1000, 0, 1000, 0.0])
        data = np.stack([np.full((2, 2, 2), v) for v in (2, 9, 4, 9, 6.0)], -1)
        _save(tmp_path / "d.nii.gz", data)
        np.savetxt(tmp_path / "d.bval", bvals[None])
        out = nd.mean_b0_numpy(tmp_path / "d.nii.gz", tmp_path / "d.bval",
                               tmp_path / "o" / "b0.nii.gz")
        assert nib.load(str(out)).get_fdata()[0, 0, 0] == pytest.approx(4.0)

    def test_select_shells_noop_and_too_few(self, tmp_path):
        d, bval, bvec = self._inputs(tmp_path)
        sel = nd.select_shells(d, bvec, bval, 1000, tmp_path / "s")
        assert sel["dwi"] == d
        with pytest.raises(ValueError, match="only 0"):
            nd.select_shells(d, bvec, bval, 10, tmp_path / "s2")


# ------------------------------------------------------------------ fitting
class TestTensorFit:

    def _truth(self, n=50, seed=0):
        rng = np.random.default_rng(seed)
        bvals, bvecs = make_scheme()
        lam = np.array([1.7, 0.35, 0.30]) * 1e-3
        Q = np.linalg.qr(rng.standard_normal((n, 3, 3)))[0]
        D = np.einsum("nij,j,nkj->nik", Q, lam, Q)
        sig = np.exp(-bvals[None] * np.einsum("mi,nij,mj->nm", bvecs, D, bvecs))
        return sig, bvals, bvecs, lam

    @pytest.mark.parametrize("method", ["ols", "wls"])
    def test_noise_free_recovery(self, method):
        sig, bvals, bvecs, lam = self._truth()
        x = nd.fit_tensor_signals(sig, bvals, bvecs, method)
        eig = nd.tensor_eigenvalues(x)
        np.testing.assert_allclose(eig, np.tile(lam, (len(eig), 1)), rtol=1e-6)

    def test_bad_method(self):
        with pytest.raises(ValueError):
            nd.fit_tensor_signals(np.ones((2, 33)), *make_scheme(), "nlls")

    def test_wls_is_better_than_ols_at_low_snr(self):
        bvals, bvecs = make_scheme()
        o = simulate_metrics([1.7, .35, .30], 10, 20000, bvals, bvecs,
                             np.random.default_rng(3), "ols")
        w = simulate_metrics([1.7, .35, .30], 10, 20000, bvals, bvecs,
                             np.random.default_rng(3), "wls")
        assert w["frac_invalid"] < o["frac_invalid"]
        assert w["NDE"].std() < o["NDE"].std()
        truth = 0.711
        assert abs(w["NDE"].mean() - truth) < abs(o["NDE"].mean() - truth)


# ------------------------------------------------------------------ metrics
class TestMetrics:

    def test_nde_known_values(self):
        e = np.array([[1, 1, 1], [1.5, .5, .5], [1.5, .75, .25]])
        np.testing.assert_allclose(nd.nde_from_eigenvalues(e),
                                   [1.0, 0.8650, 0.8173], atol=1e-3)

    def test_nde_scale_invariant_and_zero_safe(self):
        e = np.array([[1.5, .75, .25]])
        np.testing.assert_allclose(nd.nde_from_eigenvalues(e),
                                   nd.nde_from_eigenvalues(1000 * e))
        assert np.isfinite(nd.nde_from_eigenvalues(np.array([[1.0, 0, 0]]))).all()

    def test_mode_extremes(self):
        m = nd.tensor_mode(np.array([[1, 0.2, 0.2], [1, 1, 0.1],
                                     [2, 1, 0], [1, 1, 1]], dtype=float))
        assert m[0] == pytest.approx(1.0) and m[1] == pytest.approx(-1.0)
        assert m[2] == pytest.approx(0.0, abs=1e-9) and np.isnan(m[3])

    def test_westin_sums_to_one(self):
        e = np.random.default_rng(1).uniform(.1, 2, (100, 3))
        e = np.sort(e, axis=1)[:, ::-1]
        cl, cp, cs = nd.westin_indices(e)
        np.testing.assert_allclose(cl + cp + cs, 1.0)

    def test_ordering_and_trace_checks(self, tmp_path):
        shp = (3, 3, 3)
        l1, l2, l3 = np.full(shp, 1.5), np.full(shp, .6), np.full(shp, .4)
        mask = np.ones(shp, bool)
        assert nd.ordering_violations(l1, l2, l3, mask)["l2_lt_l3"] == 0
        assert nd.ordering_violations(l1, l3, l2, mask)["l2_lt_l3"] == 27
        # Dxx, Dxy, Dxz, Dyy, Dyz, Dzz with trace == 2.5
        t = np.stack([np.full(shp, v) for v in (1.0, 0, 0, 0.8, 0, 0.7)], -1)
        _save(tmp_path / "t.nii.gz", t)
        assert nd.trace_consistency(tmp_path / "t.nii.gz", l1, l2, l3, mask) \
            == pytest.approx(0.0, abs=1e-6)
        assert nd.trace_consistency(tmp_path / "t.nii.gz", l1 * 2, l2, l3, mask) > 0.1
        assert nd.trace_consistency(tmp_path / "none.nii.gz", l1, l2, l3, mask) is None

    def test_snr_map(self, tmp_path):
        shp = (3, 3, 3)
        _save(tmp_path / "b0.nii.gz", np.full(shp, 100.0))
        sig = np.full(shp, 5.0); sig[0, 0, 0] = 0
        _save(tmp_path / "n.nii.gz", sig)
        _save(tmp_path / "m.nii.gz", np.ones(shp))
        r = nd.compute_snr_map(tmp_path / "b0.nii.gz", tmp_path / "n.nii.gz",
                               tmp_path / "m.nii.gz", tmp_path / "snr.nii.gz")
        v = nib.load(str(r["path"])).get_fdata()
        assert v[1, 1, 1] == pytest.approx(20.0) and v[0, 0, 0] == 0
        assert r["stats"]["snr_median"] == pytest.approx(20.0)
        _save(tmp_path / "bad.nii.gz", np.ones((2, 2, 2)))
        with pytest.raises(ValueError):
            nd.compute_snr_map(tmp_path / "bad.nii.gz", tmp_path / "n.nii.gz",
                               tmp_path / "m.nii.gz", tmp_path / "x.nii.gz")


# ------------------------------------------------------------------ phantom
@pytest.fixture(scope="module")
def phantom(tmp_path_factory):  # noqa: D401
    d = tmp_path_factory.mktemp("phantom")
    return make_phantom(d, snr=40, seed=1), d

class TestPhantomRecovery:

    @pytest.mark.parametrize("method", ["ols", "wls"])
    def test_native_fit_recovers_truth_and_nde(self, phantom, method):
        p, d = phantom
        out = nd.estimate_tensor_numpy(p["dwi"], p["bvec"], p["bval"],
                                       p["mask"], d / f"fit_{method}", method)
        labels = nib.load(str(p["truth_labels"])).get_fdata()
        l = [nib.load(str(out[k])).get_fdata() for k in ("l1", "l2", "l3")]
        for lab, ev in TISSUES.items():
            sel = labels == lab
            est = np.stack([x[sel] for x in l], 1) * 1e3
            true_nde = nd.nde_from_eigenvalues(np.array([ev]))[0]
            est_nde = nd.nde_from_eigenvalues(np.clip(est, 1e-6, None)).mean()
            if lab == 3:
                # CSF at b=1000 sits near the Rician noise floor (signal
                # ~5 % of S0): eigenvalues are biased, but it must stay
                # clearly isotropic (high NDE).
                assert est_nde > 0.95
                continue
            assert abs(est.mean(0)[0] - ev[0]) < 0.12 * ev[0]
            assert abs(est_nde - true_nde) < 0.05, (lab, est_nde, true_nde)

    def test_shape_metrics_exclude_flagged_voxels(self, phantom):
        p, d = phantom
        out = nd.estimate_tensor_numpy(p["dwi"], p["bvec"], p["bval"],
                                       p["mask"], d / "shape_fit", "wls")
        flag = np.zeros(nib.load(str(p["mask"])).shape, np.uint8)
        flag[5, 5, 3] = 1
        _save(d / "flag.nii.gz", flag)
        r = nd.compute_shape_metrics(out["l1"], out["l2"], out["l3"],
                                     p["mask"], d / "shape", d / "flag.nii.gz")
        mode = nib.load(str(r["paths"]["mode"])).get_fdata()
        assert np.isnan(mode[5, 5, 3]) and np.isfinite(mode[6, 6, 3])
        assert r["stats"]["mode_mean"] > 0      # WM block is linear (mode ~ +1)
