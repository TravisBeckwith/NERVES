"""Tests for nerves_reliability.py and nerves_qc.py (no external tools)."""
import json
import subprocess
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import nerves_dti as nd  # noqa: E402
import nerves_qc as qc  # noqa: E402
import nerves_reliability as rel  # noqa: E402
from validation.make_phantom import make_phantom  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


# ------------------------------------------------------------ reliability
class TestReliability:
    def test_deterministic_and_cached(self):
        a = rel.expected_error(20, 30, 1000, "WM")
        b = rel.expected_error(20, 30, 1000, "WM")
        assert a == b

    def test_error_decreases_with_snr(self):
        lo = rel.expected_error(10, tissue="GM")
        hi = rel.expected_error(40, tissue="GM")
        assert hi["FA"]["sd"] < lo["FA"]["sd"]
        assert abs(hi["FA"]["bias"]) < abs(lo["FA"]["bias"])

    def test_low_snr_fa_bias_is_upward_in_low_fa_tissue(self):
        r = rel.expected_error(10, tissue="GM")
        assert r["FA"]["bias"] > 0.03          # known upward FA bias

    def test_md_is_nearly_unbiased(self):
        r = rel.expected_error(15, tissue="WM")
        assert abs(r["MD"]["bias"]) < 0.02

    def test_mode_unreliable_in_low_fa_tissue(self):
        gm = rel.expected_error(25, tissue="GM")
        wm = rel.expected_error(25, tissue="WM")
        assert gm["MODE"]["sd"] > 5 * wm["MODE"]["sd"]

    def test_more_directions_reduce_variance(self):
        a = rel.expected_error(20, n_dirs=15, tissue="WM")
        b = rel.expected_error(20, n_dirs=60, tissue="WM")
        assert b["FA"]["sd"] < a["FA"]["sd"]

    def test_min_snr_monotone_in_tolerance(self):
        tight = rel.min_snr_for_bias("FA", 0.02, tissue="WM")
        loose = rel.min_snr_for_bias("FA", 0.05, tissue="WM")
        assert loose is not None and (tight is None or tight >= loose)

    def test_min_snr_unreachable_returns_none(self):
        assert rel.min_snr_for_bias("FA", 1e-5) is None

    def test_custom_tissue_and_validation(self):
        assert rel.expected_error(20, tissue=[1.5, 0.4, 0.3])["FA"]["truth"] > 0.7
        with pytest.raises(ValueError):
            rel.expected_error(20, tissue="BONE")
        with pytest.raises(ValueError):
            rel.expected_error(0)
        with pytest.raises(ValueError):
            rel.expected_error(20, n_dirs=5)
        with pytest.raises(ValueError):
            rel.min_snr_for_bias("XX", 0.1)

    def test_roi_sd(self):
        assert rel.roi_sd(0.1, 100) == pytest.approx(0.01)
        assert rel.roi_sd(0.1, 0) == pytest.approx(0.1)

    def test_cli(self, capsys):
        assert rel.main(["--sweep-snr", "15", "30", "--n", "500"]) == 0
        out = capsys.readouterr().out
        assert "| 15 |" in out and "FA bias" in out
        assert rel.main(["--required", "FA:0.05", "--n", "500"]) == 0
        assert rel.main(["--required", "FA:1e-6", "--n", "500"]) == 1


# --------------------------------------------------------------------- QC
def _write(path, arr, affine=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(arr.astype(np.float32),
                             np.eye(4) if affine is None else affine),
             str(path))


def make_subject(root, sid, seed=0, wf="qsirecon-DSIStudio", corrupt=0.0,
                 mm2s=True, fsl_style=False, shape=(12, 12, 6), md_shift=1.0,
                 with_eig=True):
    """Write tensor scalar maps for random WM/GM voxels. ``corrupt`` is the
    fraction of voxels given non-physical values."""
    rng = np.random.default_rng(seed)
    n = int(np.prod(shape))
    lam = np.zeros((n, 3))
    wm = rng.random(n) < 0.5
    lam[wm] = [1.7, 0.35, 0.30]
    lam[~wm] = [0.95, 0.75, 0.65]
    lam *= md_shift * rng.normal(1.0, 0.03, (n, 1))
    lam.sort(axis=1)
    lam = lam[:, ::-1]
    if corrupt:
        bad = rng.random(n) < corrupt
        lam[bad, 2] = -0.1
    fa, md = nd.fa_md(np.clip(lam, 1e-9, None))
    md = lam.mean(axis=1)
    ad, rd = lam[:, 0], lam[:, 1:].mean(axis=1)
    k = 1e-3 if mm2s else 1.0
    maps = {"fa": fa, "md": md * k, "ad": ad * k, "rd": rd * k,
            "l1": lam[:, 0] * k, "l2": lam[:, 1] * k, "l3": lam[:, 2] * k}
    for p, v in maps.items():
        if p in ("l1", "l2", "l3") and not with_eig:
            continue
        if fsl_style:
            name = {"fa": "dti_FA", "md": "dti_MD", "l1": "dti_L1",
                    "l2": "dti_L2", "l3": "dti_L3"}.get(p)
            if name is None:
                continue
            _write(root / sid / f"{name}.nii.gz", v.reshape(shape))
        else:
            _write(root / wf / sid / "dwi" /
                   f"{sid}_space-ACPC_model-tensor_param-{p}_dwimap.nii.gz",
                   v.reshape(shape))
    return root


class TestDiscovery:
    def test_qsirecon_and_fsl_layouts(self, tmp_path):
        make_subject(tmp_path, "sub-01")
        make_subject(tmp_path, "sub-02", fsl_style=True)
        found = qc.discover(tmp_path)
        labels = sorted(found)
        assert any(l.startswith("qsirecon-DSIStudio/sub-01") for l in labels)
        assert "sub-02" in labels
        q = next(v for k, v in found.items() if k.startswith("qsirecon"))
        assert {"fa", "md", "ad", "rd"} <= set(q)

    def test_sets_without_fa_md_ignored(self, tmp_path):
        _write(tmp_path / "x" / "sub-1_param-ad_dwimap.nii.gz", np.ones((3, 3, 3)))
        assert qc.discover(tmp_path) == {}


class TestAnalyseSubject:
    def _paths(self, tmp_path, **kw):
        make_subject(tmp_path, "sub-01", **kw)
        return next(iter(qc.discover(tmp_path).values()))

    def test_clean_subject_passes(self, tmp_path):
        r = qc.analyse_subject(self._paths(tmp_path))
        assert r["status"] == "PASS", r["flags"]
        assert r["frac_nonphysical"] == 0
        assert r["consistency_md_ad_rd"] < 0.01
        assert 0.7 < r["md_median"] < 0.85       # rescaled to 1e-3 mm^2/s
        assert r["md_units_scaled"] is True
        assert "nde_median" in r

    def test_unit_autodetect(self, tmp_path):
        r = qc.analyse_subject(self._paths(tmp_path, mm2s=False))
        assert r["md_units_scaled"] is False
        assert 0.7 < r["md_median"] < 0.85

    def test_corrupt_subject_fails(self, tmp_path):
        r = qc.analyse_subject(self._paths(tmp_path, corrupt=0.2))
        assert r["status"] == "FAIL"
        assert r["frac_nonphysical"] > 0.1
        assert any("non-physical" in f for f in r["flags"])

    def test_empty_maps_fail(self, tmp_path):
        make_subject(tmp_path, "sub-01")
        paths = next(iter(qc.discover(tmp_path).values()))
        for p in paths.values():
            _write(p, np.zeros((12, 12, 6)))
        r = qc.analyse_subject(paths)
        assert r["status"] == "FAIL" and r["n_voxels"] == 0

    def test_inconsistent_md_flagged(self, tmp_path):
        paths = self._paths(tmp_path)
        md = np.asarray(nib.load(str(paths["md"])).dataobj) * 1.5
        _write(paths["md"], md)
        r = qc.analyse_subject(paths)
        assert any("inconsistent" in f for f in r["flags"])
        assert r["status"] == "WARN"

    def test_shape_mismatch_raises(self, tmp_path):
        paths = self._paths(tmp_path)
        _write(paths["md"], np.ones((5, 5, 5)))
        with pytest.raises(ValueError):
            qc.analyse_subject(paths)

    def test_snr_attaches_expected_error_and_warns_when_low(self, tmp_path):
        paths = self._paths(tmp_path)
        hi = qc.analyse_subject(paths, snr=40)
        lo = qc.analyse_subject(paths, snr=8)
        assert "exp_GM_FA_bias" in hi and hi["status"] == "PASS"
        assert lo["exp_GM_FA_bias"] > hi["exp_GM_FA_bias"]
        assert lo["status"] == "WARN"
        assert any("FA bias" in f for f in lo["flags"])

    def test_no_eigenvalue_maps_still_works(self, tmp_path):
        r = qc.analyse_subject(self._paths(tmp_path, with_eig=False))
        assert r["status"] == "PASS" and "nde_median" not in r

    def test_fsl_style_ordering_check(self, tmp_path):
        make_subject(tmp_path, "sub-01", fsl_style=True)
        paths = next(iter(qc.discover(tmp_path).values()))
        l1 = np.asarray(nib.load(str(paths["l1"])).dataobj)
        l3 = np.asarray(nib.load(str(paths["l3"])).dataobj)
        _write(paths["l1"], l3)             # swap L1 <-> L3
        _write(paths["l3"], l1)
        r = qc.analyse_subject(paths)
        assert r["frac_misordered_eigenvalues"] > 0.5
        assert any("not sorted" in f for f in r["flags"])


class TestOutliers:
    def test_robust_z(self):
        z = qc.robust_z(np.array([1, 1.1, 0.9, 1.0, 1.05, 10.0]))
        assert abs(z[-1]) > 3.5 and np.abs(z[:-1]).max() < 2
        assert (qc.robust_z(np.ones(5)) == 0).all()

    def test_outlier_subject_flagged_in_group(self, tmp_path):
        for i in range(7):
            make_subject(tmp_path, f"sub-{i:02d}", seed=i,
                         md_shift=3.0 if i == 6 else 1.0)
        rows = qc.run(tmp_path, tmp_path / "qc")
        bad = next(r for r in rows if "sub-06" in r["label"])
        assert any("group outlier" in f for f in bad["flags"])
        assert bad["status"] in ("WARN", "FAIL")
        good = [r for r in rows if "sub-06" not in r["label"]]
        assert all(not any("group outlier" in f for f in r["flags"])
                   for r in good)

    def test_small_group_not_flagged(self, tmp_path):
        for i in range(3):
            make_subject(tmp_path, f"sub-{i:02d}", seed=i,
                         md_shift=3.0 if i == 2 else 1.0)
        rows = qc.run(tmp_path, tmp_path / "qc")
        assert not any("group outlier" in f for r in rows for f in r["flags"])


class TestOutputsAndCli:
    def test_outputs_written(self, tmp_path):
        for i in range(2):
            make_subject(tmp_path, f"sub-{i:02d}", seed=i)
        out = tmp_path / "qc"
        qc.run(tmp_path, out, snr=25)
        tsv = (out / "nerves_qc_summary.tsv").read_text().splitlines()
        assert len(tsv) == 3
        assert "status" in tsv[0] and "flags" in tsv[0]
        rows = json.loads((out / "nerves_qc_summary.json").read_text())
        assert len(rows) == 2
        html = (out / "nerves_qc_report.html").read_text()
        assert "NERVES tensor-fit QC" in html and "sub-00" in html

    def test_cli_exit_codes(self, tmp_path):
        make_subject(tmp_path, "sub-00")
        r = subprocess.run([sys.executable, str(ROOT / "nerves_qc.py"),
                            str(tmp_path)], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        make_subject(tmp_path, "sub-bad", corrupt=0.3, seed=9)
        r = subprocess.run([sys.executable, str(ROOT / "nerves_qc.py"),
                            str(tmp_path)], capture_output=True, text=True)
        assert r.returncode == 3
        assert "FAIL" in r.stdout
        empty = tmp_path / "empty"
        empty.mkdir()
        r = subprocess.run([sys.executable, str(ROOT / "nerves_qc.py"),
                            str(empty)], capture_output=True, text=True)
        assert r.returncode == 1

    def test_unanalysable_subject_does_not_abort(self, tmp_path):
        make_subject(tmp_path, "sub-00")
        make_subject(tmp_path, "sub-01", seed=3)
        victim = next(p for p in tmp_path.rglob("*sub-01*param-md*"))
        _write(victim, np.ones((4, 4, 4)))
        rows = qc.run(tmp_path, tmp_path / "qc")
        assert {r["status"] for r in rows} >= {"FAIL"}
        assert any("could not analyse" in f for r in rows for f in r["flags"])


class TestPipelineIntegration:
    def test_qc_reads_pipeline_output(self, tmp_path):
        p = make_phantom(tmp_path / "ph", snr=40, seed=2)
        out = tmp_path / "o"
        r = subprocess.run(
            [sys.executable, str(ROOT / "nde_pipeline.py"),
             "--skip_dwi_preproc", "--tensor_fit", "python", "--wls",
             "--dwi", str(p["dwi"]), "--bval", str(p["bval"]),
             "--bvec", str(p["bvec"]), "--mask", str(p["mask"]),
             "--out", str(out)], capture_output=True, text=True)
        assert r.returncode in (0, 2), r.stderr
        rows = qc.run(out, tmp_path / "qc", snr=40)
        assert len(rows) == 1
        row = rows[0]
        assert row["n_voxels"] > 100
        assert row["frac_nonphysical"] < 0.05
        assert "nde_median" in row and "exp_GM_FA_bias" in row
