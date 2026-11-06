"""
tests/test_pipeline.py
-----------------------
Tests against the actual nde_pipeline.py functions (not a mirrored
reimplementation), covering code paths that previously had zero test
coverage: compute_nde as actually shipped, write_report, validate_tools,
CLI argument validation, the --readout_time fix, and a full
--skip_preproc run through the CLI end-to-end.

None of these require FSL or MRtrix3 to be installed — --skip_preproc
bypasses preprocess_dwi/preprocess_t1/estimate_tensor entirely, which is
exactly what makes a real end-to-end run possible here without those
external tools. See tests/test_nde_math.py for the pure-math boundary
condition and manuscript Table 1 tests.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import nibabel as nib
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

import nde_pipeline
from nde_pipeline import compute_nde, write_report, validate_tools, check_tool

PIPELINE_SCRIPT = Path(__file__).parent.parent / "nde_pipeline.py"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _save_nii(path: Path, arr: np.ndarray, affine=None):
    nib.save(nib.Nifti1Image(arr.astype(np.float32), affine if affine is not None else np.eye(4)), str(path))


def _make_skip_preproc_inputs(tmp_path: Path, shape=(6, 6, 4)):
    """Build a synthetic set of dtifit-style eigenvalue/FA/MD/PVE maps
    and a spherical brain mask, suitable for --skip_preproc.
    """
    rng = np.random.default_rng(0)
    zz, yy, xx = np.meshgrid(
        np.arange(shape[2]), np.arange(shape[1]), np.arange(shape[0]), indexing="ij"
    )
    center = np.array(shape) / 2
    dist = np.sqrt((xx - center[0]) ** 2 + (yy - center[1]) ** 2 + (zz - center[2]) ** 2)
    mask = (dist.transpose(2, 1, 0) < min(shape) / 2 - 0.5).astype(np.float32)

    # Plausible WM-like eigenvalues (L1 > L2 >= L3 > 0) inside the mask
    l1 = np.where(mask > 0, 1.2 + 0.1 * rng.standard_normal(shape), 0.0)
    l2 = np.where(mask > 0, 0.6 + 0.05 * rng.standard_normal(shape), 0.0)
    l3 = np.where(mask > 0, 0.4 + 0.05 * rng.standard_normal(shape), 0.0)
    l1, l2 = np.maximum(l1, l2 + 0.01), l2  # enforce L1 > L2
    l2, l3 = np.maximum(l2, l3 + 0.01), l3  # enforce L2 > L3

    fa = np.where(mask > 0, 0.45 + 0.05 * rng.standard_normal(shape), 0.0)
    md = np.where(mask > 0, 0.0008 + 0.0001 * rng.standard_normal(shape), 0.0)

    pve_wm = np.where(mask > 0, 0.95, 0.0).astype(np.float32)
    pve_gm = np.where(mask > 0, 0.02, 0.0).astype(np.float32)
    pve_csf = np.where(mask > 0, 0.02, 0.0).astype(np.float32)

    tensor_dir = tmp_path / "tensor_in"
    tensor_dir.mkdir()
    paths = {}
    for name, arr in [("dti_L1", l1), ("dti_L2", l2), ("dti_L3", l3),
                       ("dti_FA", fa), ("dti_MD", md), ("mask", mask)]:
        p = tensor_dir / f"{name}.nii.gz"
        _save_nii(p, arr)
        paths[name] = p

    pve_dir = tmp_path / "pve_in"
    pve_dir.mkdir()
    for name, arr in [("pve_wm", pve_wm), ("pve_gm", pve_gm), ("pve_csf", pve_csf)]:
        p = pve_dir / f"{name}.nii.gz"
        _save_nii(p, arr)
        paths[name] = p

    return paths


# ---------------------------------------------------------------------------
# compute_nde — against the real function, not a mirror
# ---------------------------------------------------------------------------

class TestComputeNdeReal:
    """Regression tests against the actual nde_pipeline.compute_nde.

    tests/test_nde_math.py checks a hand-mirrored copy of this logic;
    these tests import and call the real function directly so a change to
    compute_nde that isn't mirrored in the test file can't slip through.
    """

    TABLE1 = [
        ("A", 1.5, 0.50, 0.50, 0.8650),
        ("B", 1.5, 0.75, 0.25, 0.8173),
        ("C", 1.2, 0.60, 0.60, 0.9464),
        ("D", 1.2, 0.80, 0.40, 0.9207),
    ]

    def test_matches_manuscript_table1(self, tmp_path):
        shape = (len(self.TABLE1), 1, 1)
        l1 = np.array([[[v[1]]] for v in self.TABLE1])
        l2 = np.array([[[v[2]]] for v in self.TABLE1])
        l3 = np.array([[[v[3]]] for v in self.TABLE1])
        mask = np.ones(shape, dtype=np.float32)

        for name, arr in [("l1", l1), ("l2", l2), ("l3", l3), ("mask", mask)]:
            _save_nii(tmp_path / f"{name}.nii.gz", arr)

        result = compute_nde(
            tmp_path / "l1.nii.gz", tmp_path / "l2.nii.gz", tmp_path / "l3.nii.gz",
            tmp_path / "mask.nii.gz", tmp_path / "out",
        )
        nde_vol = nib.load(str(result["nde"])).get_fdata()
        for i, (name, _, _, _, expected) in enumerate(self.TABLE1):
            got = float(nde_vol[i, 0, 0])
            assert abs(got - expected) < 0.001, f"Config {name}: got {got}, expected {expected}"

    def test_nonphysical_eigenvalue_flagged(self, tmp_path):
        shape = (2, 1, 1)
        l1 = np.array([[[1.5]], [[-0.1]]])  # voxel 1: non-physical
        l2 = np.array([[[0.5]], [[0.3]]])
        l3 = np.array([[[0.5]], [[0.2]]])
        mask = np.ones(shape, dtype=np.float32)
        for name, arr in [("l1", l1), ("l2", l2), ("l3", l3), ("mask", mask)]:
            _save_nii(tmp_path / f"{name}.nii.gz", arr)

        result = compute_nde(
            tmp_path / "l1.nii.gz", tmp_path / "l2.nii.gz", tmp_path / "l3.nii.gz",
            tmp_path / "mask.nii.gz", tmp_path / "out",
        )
        nde_vol = nib.load(str(result["nde"])).get_fdata()
        qc_vol = nib.load(str(result["qc_flag"])).get_fdata()
        assert np.isnan(nde_vol[1, 0, 0])   # undefined, NOT 0 (=fully anisotropic)
        assert qc_vol[1, 0, 0] == 1.0
        assert qc_vol[0, 0, 0] == 0.0
        assert result["stats"]["n_nonphysical"] == 1

    def test_flagged_voxels_excluded_from_stats(self, tmp_path):
        """Non-physical voxels must not drag summary statistics toward 0."""
        shape = (3, 1, 1)
        l1 = np.array([[[1.2]], [[1.2]], [[-0.1]]])   # voxel 2 non-physical
        l2 = np.full(shape, 0.6)
        l3 = np.full(shape, 0.4)
        mask = np.ones(shape, dtype=np.float32)
        for name, arr in [("l1", l1), ("l2", l2), ("l3", l3), ("mask", mask)]:
            _save_nii(tmp_path / f"{name}.nii.gz", arr)
        result = compute_nde(
            tmp_path / "l1.nii.gz", tmp_path / "l2.nii.gz", tmp_path / "l3.nii.gz",
            tmp_path / "mask.nii.gz", tmp_path / "out",
        )
        st = result["stats"]
        assert st["n_valid_voxels"] == 2 and st["n_nonphysical"] == 1
        assert st["nde_min"] > 0.5            # would be 0.0 if flagged voxel leaked in
        assert abs(st["nde_mean"] - st["nde_max"]) < 1e-6

    def test_all_invalid_raises(self, tmp_path):
        shape = (2, 1, 1)
        for name, arr in [("l1", np.full(shape, -1.0)), ("l2", np.full(shape, 0.5)),
                          ("l3", np.full(shape, 0.5)), ("mask", np.ones(shape))]:
            _save_nii(tmp_path / f"{name}.nii.gz", arr)
        with pytest.raises(RuntimeError, match="No valid"):
            compute_nde(
                tmp_path / "l1.nii.gz", tmp_path / "l2.nii.gz",
                tmp_path / "l3.nii.gz", tmp_path / "mask.nii.gz", tmp_path / "out",
            )

    def test_stats_dict_keys(self, tmp_path):
        shape = (3, 1, 1)
        l1 = np.full(shape, 1.2)
        l2 = np.full(shape, 0.6)
        l3 = np.full(shape, 0.4)
        mask = np.ones(shape, dtype=np.float32)
        for name, arr in [("l1", l1), ("l2", l2), ("l3", l3), ("mask", mask)]:
            _save_nii(tmp_path / f"{name}.nii.gz", arr)

        result = compute_nde(
            tmp_path / "l1.nii.gz", tmp_path / "l2.nii.gz", tmp_path / "l3.nii.gz",
            tmp_path / "mask.nii.gz", tmp_path / "out",
        )
        expected_keys = {
            "n_mask_voxels", "n_valid_voxels", "n_nonphysical", "nonphysical_pct",
            "nde_mean", "nde_std", "nde_median", "nde_min", "nde_max",
            "nde_p5", "nde_p95",
        }
        assert expected_keys <= set(result["stats"].keys())


# ---------------------------------------------------------------------------
# write_report
# ---------------------------------------------------------------------------

class TestWriteReport:

    def test_writes_valid_json(self, tmp_path):
        write_report(
            out_dir=tmp_path,
            nde_stats={"nde_mean": 0.85, "nonphysical_pct": 0.5},
            qc_results={"csf_nde_pass": True},
            args_dict={"dwi": "x.nii.gz"},
        )
        report_path = tmp_path / "nde_pipeline_report.json"
        assert report_path.exists()
        report = json.loads(report_path.read_text())
        assert report["nde_statistics"]["nde_mean"] == 0.85
        assert report["qc_results"]["csf_nde_pass"] is True
        assert "qc_thresholds" in report
        assert report["qc_thresholds"]["csf_nde_mean_min"] == 0.90
        assert "reference" in report and "Fozouni" in report["reference"]
        assert "nde_fa_r_max" not in report["qc_thresholds"]   # r is informational
        assert "notes" in report


# ---------------------------------------------------------------------------
# run(): no shell, safe with awkward paths
# ---------------------------------------------------------------------------

class TestRun:

    def test_arguments_passed_verbatim(self, tmp_path):
        """Spaces and shell metacharacters must not be interpreted."""
        marker = tmp_path / "pwned"
        tricky = f"a b; touch {marker} #"
        res = nde_pipeline.run(
            [sys.executable, "-c", "import sys; print(sys.argv[1])", tricky])
        assert res.stdout.strip() == tricky
        assert not marker.exists()

    def test_failure_raises(self):
        with pytest.raises(RuntimeError):
            nde_pipeline.run([sys.executable, "-c", "raise SystemExit(3)"])

    def test_no_shell_pipelines_in_source(self):
        src = PIPELINE_SCRIPT.read_text()
        assert "shell=True" not in src


# ---------------------------------------------------------------------------
# validate_tools / check_tool
# ---------------------------------------------------------------------------

class TestValidateTools:

    def test_check_tool_true_for_python(self):
        # `python3` (or `python`) must exist in any environment running
        # this test, so use it as a known-present tool.
        assert check_tool(sys.executable.rsplit("/", 1)[-1]) or check_tool("python3")

    def test_check_tool_false_for_nonexistent(self):
        assert check_tool("definitely_not_a_real_binary_xyz123") is False

    def test_validate_tools_raises_when_missing(self, monkeypatch):
        monkeypatch.setattr(nde_pipeline, "check_tool", lambda name: False)
        with pytest.raises(EnvironmentError, match="Missing required tools"):
            validate_tools()

    def test_validate_tools_passes_when_present(self, monkeypatch):
        monkeypatch.setattr(nde_pipeline, "check_tool", lambda name: True)
        validate_tools()  # should not raise


# ---------------------------------------------------------------------------
# --readout_time threading (the fix for the previously-hardcoded value)
# ---------------------------------------------------------------------------

class TestReadoutTime:

    def test_custom_readout_time_written_to_acqparams(self, tmp_path, monkeypatch):
        """preprocess_dwi must write the caller-supplied readout_time into
        acqparams.txt, not the hardcoded 0.05 default.
        """
        out_dir = tmp_path / "preproc"
        out_dir.mkdir()

        # Pre-place the one file preprocess_dwi actually reads with nib
        # (dwi_degibbs.nii.gz, to get n_vols for index.txt) — every other
        # step just shells out via run(), which we no-op below.
        _save_nii(out_dir / "dwi_degibbs.nii.gz", np.zeros((4, 4, 4, 3)))

        def fake_run(cmd, cwd=None, check=True):
            class _Result:
                returncode = 0
                stdout = ""
                stderr = ""
            return _Result()

        monkeypatch.setattr(nde_pipeline, "run", fake_run)
        monkeypatch.setattr(nde_pipeline, "check_tool", lambda name: False)  # skip topup branch's tool check

        nde_pipeline.preprocess_dwi(
            dwi=tmp_path / "dwi.nii.gz",
            bval=tmp_path / "dwi.bval",
            bvec=tmp_path / "dwi.bvec",
            out_dir=out_dir,
            pe_dir="AP",
            rpe_b0=None,
            readout_time=0.077,
        )

        acqparams = (out_dir / "acqparams.txt").read_text()
        assert "0.077" in acqparams
        assert "0.05" not in acqparams  # old hardcoded default must be gone

    def test_default_readout_time_is_005(self, tmp_path, monkeypatch):
        out_dir = tmp_path / "preproc"
        out_dir.mkdir()
        _save_nii(out_dir / "dwi_degibbs.nii.gz", np.zeros((4, 4, 4, 3)))

        def fake_run(cmd, cwd=None, check=True):
            class _Result:
                returncode = 0
                stdout = ""
                stderr = ""
            return _Result()

        monkeypatch.setattr(nde_pipeline, "run", fake_run)
        monkeypatch.setattr(nde_pipeline, "check_tool", lambda name: False)

        nde_pipeline.preprocess_dwi(
            dwi=tmp_path / "dwi.nii.gz", bval=tmp_path / "dwi.bval",
            bvec=tmp_path / "dwi.bvec", out_dir=out_dir, pe_dir="AP", rpe_b0=None,
        )
        acqparams = (out_dir / "acqparams.txt").read_text()
        assert "0.05" in acqparams


# ---------------------------------------------------------------------------
# CLI argument validation (subprocess — real script invocation)
# ---------------------------------------------------------------------------

class TestCLIArgValidation:

    def test_skip_preproc_missing_required_exits_1(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT),
             "--skip_preproc",
             "--dwi", "x", "--bval", "x", "--bvec", "x", "--t1", "x",
             "--out", str(tmp_path / "out")],
            capture_output=True, text=True,
        )
        assert result.returncode == 1
        assert "--l1" in result.stderr or "--l1" in result.stdout

    def test_help_exits_0(self):
        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT), "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "--readout_time" in result.stdout


# ---------------------------------------------------------------------------
# Full end-to-end run via --skip_preproc (no FSL/MRtrix3 required)
# ---------------------------------------------------------------------------

class TestSkipPreprocEndToEnd:
    """Runs the real CLI end-to-end through compute_nde, run_qc, and
    write_report — the parts of the pipeline that don't need FSL/MRtrix3
    when --skip_preproc is used. This is what previously had zero coverage.
    """

    def test_full_run_produces_all_outputs(self, tmp_path):
        inputs = _make_skip_preproc_inputs(tmp_path)
        out_dir = tmp_path / "nerves_out"

        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT),
             "--skip_preproc",
             "--l1", str(inputs["dti_L1"]),
             "--l2", str(inputs["dti_L2"]),
             "--l3", str(inputs["dti_L3"]),
             "--mask", str(inputs["mask"]),
             "--pve_csf", str(inputs["pve_csf"]),
             "--pve_gm", str(inputs["pve_gm"]),
             "--pve_wm", str(inputs["pve_wm"]),
             "--dwi", "x", "--bval", "x", "--bvec", "x", "--t1", "x",
             "--out", str(out_dir)],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr

        assert (out_dir / "nde" / "nde.nii.gz").exists()
        assert (out_dir / "nde" / "nde_qc_flag.nii.gz").exists()
        assert (out_dir / "nde_pipeline_report.json").exists()

        report = json.loads((out_dir / "nde_pipeline_report.json").read_text())
        assert 0.0 <= report["nde_statistics"]["nde_mean"] <= 1.0
        assert report["nde_statistics"]["n_nonphysical"] == 0

        # FA/MD auto-detected via dtifit naming convention (dti_L1 -> dti_FA/dti_MD)
        # should have triggered QC figure generation.
        qc_dir = out_dir / "qc"
        assert (qc_dir / "qc_histograms.png").exists()
        assert (qc_dir / "qc_nde_fa_scatter.png").exists()
        assert (qc_dir / "qc_axial_montage.png").exists()
        assert (qc_dir / "qc_tissue_distributions.png").exists()

    def test_shape_mismatch_exits_1(self, tmp_path):
        shape = (4, 4, 3)
        in_dir = tmp_path / "in"
        in_dir.mkdir()
        for name, shp in [("eig1", shape), ("eig2", shape), ("eig3", (5, 4, 3)),
                          ("brainmask", shape)]:
            _save_nii(in_dir / f"{name}.nii.gz", np.ones(shp))
        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT), "--skip_preproc",
             "--l1", str(in_dir / "eig1.nii.gz"), "--l2", str(in_dir / "eig2.nii.gz"),
             "--l3", str(in_dir / "eig3.nii.gz"), "--mask", str(in_dir / "brainmask.nii.gz"),
             "--dwi", "x", "--bval", "x", "--bvec", "x", "--t1", "x",
             "--out", str(tmp_path / "o")],
            capture_output=True, text=True)
        assert result.returncode == 1
        assert "Shape mismatch" in result.stderr

    def test_qc_flagged_voxels_do_not_affect_report(self, tmp_path):
        """A flagged voxel must yield NaN in the map and no NaN in the JSON."""
        inputs = _make_skip_preproc_inputs(tmp_path)
        l1 = nib.load(str(inputs["dti_L1"])).get_fdata().astype(np.float32)
        mask = nib.load(str(inputs["mask"])).get_fdata() > 0
        idx = tuple(np.argwhere(mask)[0])
        l1[idx] = -0.5
        _save_nii(inputs["dti_L1"], l1)
        out_dir = tmp_path / "flagged_out"
        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT), "--skip_preproc",
             "--l1", str(inputs["dti_L1"]), "--l2", str(inputs["dti_L2"]),
             "--l3", str(inputs["dti_L3"]), "--mask", str(inputs["mask"]),
             "--pve_csf", str(inputs["pve_csf"]), "--pve_gm", str(inputs["pve_gm"]),
             "--pve_wm", str(inputs["pve_wm"]),
             "--dwi", "x", "--bval", "x", "--bvec", "x", "--t1", "x",
             "--out", str(out_dir)],
            capture_output=True, text=True)
        # 1 flagged voxel in a tiny mask exceeds 5% -> documented exit code 2
        assert result.returncode == 2, result.stderr
        rpt = json.loads((out_dir / "nde_pipeline_report.json").read_text())
        assert rpt["nde_statistics"]["n_nonphysical"] == 1
        assert np.isfinite(rpt["qc_results"]["nde_fa_pearson_r"])
        assert "nde_fa_anticorrelated" not in rpt["qc_results"]

    def test_fa_md_not_found_skips_qc_gracefully(self, tmp_path):
        """If FA/MD aren't alongside L1 (non-dtifit naming), the pipeline
        should still complete and just skip QC figures, not crash — and
        must not mistake the L1 file itself for the FA map (regression
        test for the str.replace() no-op bug).
        """
        shape = (4, 4, 3)
        mask = np.ones(shape, dtype=np.float32)
        l1 = np.full(shape, 1.2, dtype=np.float32)
        l2 = np.full(shape, 0.6, dtype=np.float32)
        l3 = np.full(shape, 0.4, dtype=np.float32)

        in_dir = tmp_path / "in"
        in_dir.mkdir()
        for name, arr in [("eig1", l1), ("eig2", l2), ("eig3", l3), ("brainmask", mask)]:
            _save_nii(in_dir / f"{name}.nii.gz", arr)

        out_dir = tmp_path / "nerves_out2"
        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT),
             "--skip_preproc",
             "--l1", str(in_dir / "eig1.nii.gz"),
             "--l2", str(in_dir / "eig2.nii.gz"),
             "--l3", str(in_dir / "eig3.nii.gz"),
             "--mask", str(in_dir / "brainmask.nii.gz"),
             "--dwi", "x", "--bval", "x", "--bvec", "x", "--t1", "x",
             "--out", str(out_dir)],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        assert (out_dir / "nde" / "nde.nii.gz").exists()
        assert not (out_dir / "qc").exists()
        assert "FA map not found" in result.stderr

    def test_explicit_fa_md_override(self, tmp_path):
        """--fa/--md should work even when filenames don't follow the
        dtifit L1/FA naming convention.
        """
        shape = (4, 4, 3)
        mask = np.ones(shape, dtype=np.float32)
        l1 = np.full(shape, 1.2, dtype=np.float32)
        l2 = np.full(shape, 0.6, dtype=np.float32)
        l3 = np.full(shape, 0.4, dtype=np.float32)
        fa = np.full(shape, 0.45, dtype=np.float32)
        md = np.full(shape, 0.0008, dtype=np.float32)

        in_dir = tmp_path / "in"
        in_dir.mkdir()
        for name, arr in [("eig1", l1), ("eig2", l2), ("eig3", l3),
                          ("brainmask", mask), ("myFA", fa), ("myMD", md)]:
            _save_nii(in_dir / f"{name}.nii.gz", arr)

        out_dir = tmp_path / "nerves_out3"
        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT),
             "--skip_preproc",
             "--l1", str(in_dir / "eig1.nii.gz"),
             "--l2", str(in_dir / "eig2.nii.gz"),
             "--l3", str(in_dir / "eig3.nii.gz"),
             "--mask", str(in_dir / "brainmask.nii.gz"),
             "--fa", str(in_dir / "myFA.nii.gz"),
             "--md", str(in_dir / "myMD.nii.gz"),
             "--dwi", "x", "--bval", "x", "--bvec", "x", "--t1", "x",
             "--out", str(out_dir)],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        assert (out_dir / "qc" / "qc_histograms.png").exists()
        assert (out_dir / "qc" / "qc_nde_fa_scatter.png").exists()


# ---------------------------------------------------------------------------
# Command construction (no FSL/MRtrix needed: run() is captured)
# ---------------------------------------------------------------------------

class TestCommandConstruction:

    def _capture(self, monkeypatch):
        calls = []

        def fake_run(cmd, cwd=None, check=True):
            calls.append(cmd)

            class _R:
                returncode = 0
                stdout = stderr = ""
            return _R()

        monkeypatch.setattr(nde_pipeline, "run", fake_run)
        return calls

    def test_dwi_commands_are_lists_with_spaces_preserved(self, tmp_path, monkeypatch):
        calls = self._capture(monkeypatch)
        monkeypatch.setattr(nde_pipeline, "check_tool", lambda n: True)
        out_dir = tmp_path / "out dir with spaces"
        out_dir.mkdir()
        _save_nii(out_dir / "dwi_degibbs.nii.gz", np.zeros((4, 4, 4, 3)))
        rpe = tmp_path / "rpe b0.nii.gz"
        _save_nii(rpe, np.zeros((4, 4, 4)))
        nde_pipeline.preprocess_dwi(
            dwi=tmp_path / "my dwi.nii.gz", bval=tmp_path / "a.bval",
            bvec=tmp_path / "a.bvec", out_dir=out_dir, rpe_b0=rpe)
        assert all(isinstance(c, list) for c in calls)
        flat = [str(a) for c in calls for a in c]
        assert str(tmp_path / "my dwi.nii.gz") in flat      # one argv element
        assert not any(a == "|" for a in flat)               # no shell pipes
        names = [c[0] for c in calls]
        assert "applytopup" in names                         # mask from corrected b0
        eddy_name = next(n for n in names if n.startswith("eddy"))
        assert names.index("applytopup") < names.index(eddy_name)
        eddy = next(c for c in calls if c[0] == eddy_name)
        assert any(a.startswith("--topup=") for a in eddy)

    def test_ants_uses_affine_only_and_is_labelled_so(self, tmp_path, monkeypatch):
        calls = self._capture(monkeypatch)
        monkeypatch.setattr(nde_pipeline, "check_tool", lambda n: True)
        nde_pipeline.preprocess_t1(tmp_path / "t1.nii.gz", tmp_path / "ref.nii.gz",
                                   tmp_path / "t1out", use_ants=True)
        reg = next(c for c in calls if c[0] == "antsRegistrationSyNQuick.sh")
        assert reg[reg.index("-t") + 1] == "a"
