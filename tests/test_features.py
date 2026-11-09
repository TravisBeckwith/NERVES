"""Tests for the v1.1 pipeline features. None need FSL / MRtrix3."""

import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import nde_pipeline as npl
import nerves_dti as nd
import nerves_bids
import nerves_group_summary as gs
from validation.make_phantom import make_phantom
from validation import noise_simulation as ns
from validation import retest_reliability as tr

SCRIPT = ROOT / "nde_pipeline.py"


def _save(path, arr):
    nib.save(nib.Nifti1Image(np.asarray(arr, dtype=np.float32), np.eye(4)),
             str(path))


def _cli(args, **kw):
    return subprocess.run([sys.executable, str(SCRIPT)] + [str(a) for a in args],
                          capture_output=True, text=True, **kw)


@pytest.fixture(scope="module")
def phantom(tmp_path_factory):
    d = tmp_path_factory.mktemp("ph")
    return make_phantom(d, snr=40, seed=2), d


# ----------------------------------------------------------- end to end
class TestEndToEndPhantom:
    """Full CLI run on a synthetic DWI with known tensors, no external tools."""

    def _run(self, phantom, out, *extra):
        p, _ = phantom
        return _cli(["--skip_dwi_preproc", "--tensor_fit", "python",
                     "--dwi", p["dwi"], "--bval", p["bval"], "--bvec", p["bvec"],
                     "--mask", p["mask"], "--out", out, *extra])

    def test_run_recovers_wm_nde(self, phantom, tmp_path):
        p, _ = phantom
        out = tmp_path / "o"
        r = self._run(phantom, out, "--wls", "--max_b", "1500")
        assert r.returncode in (0, 2), r.stderr
        nde = nib.load(str(out / "nde" / "nde.nii.gz")).get_fdata()
        lab = nib.load(str(p["truth_labels"])).get_fdata()
        true_wm = nd.nde_from_eigenvalues(np.array([[1.7, .35, .30]]))[0]
        assert abs(np.nanmean(nde[lab == 1]) - true_wm) < 0.03
        for m in ("mode", "cl", "cp", "cs"):
            assert (out / "shape" / f"{m}.nii.gz").exists()
        rpt = json.loads((out / "nde_pipeline_report.json").read_text())
        assert rpt["run_mode"] == "dwi_preprocessed"
        assert rpt["versions"]["nerves"] == npl.__version__
        assert rpt["tensor_fit"] == {"backend": "python", "wls": True}
        assert rpt["acquisition"]["n_b0"] == 3
        assert "mode_mean" in rpt["shape_statistics"]
        assert "nde_variance_explained_by_fa" in rpt["qc_results"]

    def test_nde_only_skips_shape_maps(self, phantom, tmp_path):
        out = tmp_path / "o"
        r = self._run(phantom, out, "--nde_only")
        assert r.returncode in (0, 2), r.stderr
        assert (out / "nde" / "nde.nii.gz").exists()
        assert not (out / "shape").exists()

    def test_resume_skips_tensor_fit(self, phantom, tmp_path):
        out = tmp_path / "o"
        self._run(phantom, out)
        r = self._run(phantom, out, "--resume")
        assert r.returncode in (0, 2)
        assert "Resume: skipping tensor fit" in r.stderr

    def test_bad_gradient_table_exits_1(self, phantom, tmp_path):
        p, _ = phantom
        bad = tmp_path / "bad.bval"
        bad.write_text("0 1000 1000\n")
        r = _cli(["--skip_dwi_preproc", "--tensor_fit", "python", "--dwi", p["dwi"],
                  "--bval", bad, "--bvec", p["bvec"], "--mask", p["mask"],
                  "--out", tmp_path / "o"])
        assert r.returncode == 1 and "Invalid DWI input" in r.stderr

    def test_full_mode_requires_t1(self, phantom, tmp_path):
        p, _ = phantom
        r = _cli(["--dwi", p["dwi"], "--bval", p["bval"], "--bvec", p["bvec"],
                  "--out", tmp_path / "o"])
        assert r.returncode == 1 and "--t1" in r.stderr

    @pytest.mark.skipif(shutil.which("bet") is not None, reason="FSL present")
    def test_missing_mask_without_bet_exits_1(self, phantom, tmp_path):
        p, _ = phantom
        r = _cli(["--skip_dwi_preproc", "--tensor_fit", "python", "--dwi", p["dwi"],
                  "--bval", p["bval"], "--bvec", p["bvec"], "--out", tmp_path / "o"])
        assert r.returncode == 1 and "--mask" in r.stderr

    def test_modes_mutually_exclusive(self, tmp_path):
        r = _cli(["--skip_preproc", "--skip_dwi_preproc", "--out", tmp_path / "o"])
        assert r.returncode == 1

    def test_min_snr_masks_voxels(self, phantom, tmp_path):
        p, d = phantom
        shp = nib.load(str(p["mask"])).shape
        snr = np.full(shp, 50.0)
        snr[6, 6, 4] = 2.0                       # one WM voxel below threshold
        _save(d / "snr.nii.gz", snr)
        out = tmp_path / "o"
        r = self._run(phantom, out, "--snr", d / "snr.nii.gz", "--min_snr", "10")
        assert r.returncode in (0, 2), r.stderr
        nde = nib.load(str(out / "nde" / "nde.nii.gz")).get_fdata()
        flag = nib.load(str(out / "nde" / "nde_qc_flag.nii.gz")).get_fdata()
        assert np.isnan(nde[6, 6, 4]) and flag[6, 6, 4] == 2
        rpt = json.loads((out / "nde_pipeline_report.json").read_text())
        assert rpt["nde_statistics"]["n_low_snr"] == 1
        mode = nib.load(str(out / "shape" / "mode.nii.gz")).get_fdata()
        assert np.isnan(mode[6, 6, 4])

    def test_version_flag(self):
        r = _cli(["--version", "--out", "x"])
        assert r.returncode == 0 and npl.__version__ in r.stdout


# ----------------------------------------------------------- command logic
class _Capture:
    def __init__(self, monkeypatch, tools=True):
        self.calls = []

        def fake_run(cmd, cwd=None, check=True):
            self.calls.append([str(c) for c in cmd])

            class R:
                returncode = 0
                stdout = stderr = ""
            return R()
        monkeypatch.setattr(npl, "run", fake_run)
        monkeypatch.setattr(npl, "check_tool", lambda n: tools)

    @property
    def names(self):
        return [c[0] for c in self.calls]


class TestCommandLogic:

    def test_find_eddy_order(self, monkeypatch):
        avail = {"eddy_openmp"}
        monkeypatch.setattr(npl, "check_tool", lambda n: n in avail)
        assert npl.find_eddy() == "eddy_openmp"
        avail.add("eddy_cpu")
        assert npl.find_eddy() == "eddy_cpu"
        avail.add("eddy_cuda")
        assert npl.find_eddy() == "eddy_cpu" and npl.find_eddy(gpu=True) == "eddy_cuda"
        avail.clear()
        assert npl.find_eddy() is None

    def test_validate_tools_modes(self, monkeypatch):
        monkeypatch.setattr(npl, "check_tool", lambda n: False)
        # native fit on preprocessed DWI without T1 needs nothing external
        npl.validate_tools(tensor_fit="python", skip_dwi_preproc=True, has_t1=False)
        with pytest.raises(EnvironmentError, match="dtifit"):
            npl.validate_tools(tensor_fit="dtifit", skip_dwi_preproc=True, has_t1=False)
        with pytest.raises(EnvironmentError, match="eddy"):
            npl.validate_tools()
        monkeypatch.setattr(npl, "check_tool", lambda n: n != "epi_reg")
        with pytest.raises(EnvironmentError, match="epi_reg"):
            npl.validate_tools(registration="epi_reg")

    def test_threads_and_eddy_options(self, tmp_path, monkeypatch):
        cap = _Capture(monkeypatch)
        out = tmp_path / "p"
        out.mkdir()
        _save(out / "dwi_degibbs.nii.gz", np.zeros((3, 3, 3, 4)))
        npl.preprocess_dwi(tmp_path / "d.nii.gz", tmp_path / "b.bval",
                           tmp_path / "b.bvec", out, n_threads=4)
        assert any("-nthreads" in c and c[c.index("-nthreads") + 1] == "4"
                   for c in cap.calls if c[0] == "dwidenoise")
        eddy = next(c for c in cap.calls if c[0].startswith("eddy"))
        assert "--nthr=4" in eddy

    def test_rpe_b0_is_averaged_and_bids_pe_codes(self, tmp_path, monkeypatch):
        cap = _Capture(monkeypatch)
        out = tmp_path / "p"
        out.mkdir()
        _save(out / "dwi_degibbs.nii.gz", np.zeros((3, 3, 3, 4)))
        rpe = tmp_path / "rpe.nii.gz"
        _save(rpe, np.zeros((3, 3, 3, 2)))
        npl.preprocess_dwi(tmp_path / "d.nii.gz", tmp_path / "b.bval",
                           tmp_path / "b.bvec", out, pe_dir="j-", rpe_b0=rpe,
                           readout_time=0.05)
        tmean = next(c for c in cap.calls if c[0] == "fslmaths" and "-Tmean" in c)
        assert str(rpe) in tmean
        rows = (out / "acqparams.txt").read_text().split("\n")
        assert rows[0].startswith("0 -1 0") and rows[1].startswith("0 1 0")
        assert npl.PE_VECTORS["j-"] == npl.PE_VECTORS["AP"]

    def test_resume_skips_completed_stages(self, tmp_path, monkeypatch):
        cap = _Capture(monkeypatch)
        out = tmp_path / "p"
        out.mkdir()
        _save(out / "dwi_degibbs.nii.gz", np.zeros((3, 3, 3, 4)))
        for n in ["dwi_denoised.mif", "noise_map.nii.gz", "dwi_degibbs.mif",
                  "dwi_degibbs.bvec", "dwi_degibbs.bval", "mean_b0.nii.gz",
                  "mean_b0_brain_mask.nii.gz", "dwi_eddy.nii.gz",
                  "dwi_eddy.eddy_rotated_bvecs", "mean_b0_eddy.nii.gz",
                  "mean_b0_final_brain_mask.nii.gz"]:
            if not (out / n).exists():
                (out / n).write_text("x")
        r = npl.preprocess_dwi(tmp_path / "d.nii.gz", tmp_path / "b.bval",
                               tmp_path / "b.bvec", out, resume=True)
        assert cap.calls == []
        assert r["noise_map"].name == "noise_map.nii.gz"
        # without --resume the same directory is reprocessed
        npl.preprocess_dwi(tmp_path / "d.nii.gz", tmp_path / "b.bval",
                           tmp_path / "b.bvec", out, resume=False)
        assert any(c[0] == "dwidenoise" for c in cap.calls)

    def test_epi_reg_registration(self, tmp_path, monkeypatch):
        cap = _Capture(monkeypatch)
        npl.preprocess_t1(tmp_path / "t1.nii.gz", tmp_path / "ref.nii.gz",
                          tmp_path / "t1", registration="epi_reg",
                          dwi_ref_brain=tmp_path / "refbrain.nii.gz")
        assert "epi_reg" in cap.names and "convert_xfm" in cap.names
        er = next(c for c in cap.calls if c[0] == "epi_reg")
        assert any(a.startswith("--wmseg=") for a in er)
        assert cap.names.index("convert_xfm") < cap.names.index("flirt")
        assert all("-dof" not in c for c in cap.calls)      # no plain FLIRT fit

    def test_epi_reg_without_brain_falls_back(self, tmp_path, monkeypatch):
        cap = _Capture(monkeypatch)
        npl.preprocess_t1(tmp_path / "t1.nii.gz", tmp_path / "ref.nii.gz",
                          tmp_path / "t1", registration="epi_reg")
        assert "epi_reg" not in cap.names
        assert any("-dof" in c for c in cap.calls)

    def test_unknown_registration(self, tmp_path, monkeypatch):
        _Capture(monkeypatch)
        with pytest.raises(ValueError):
            npl.preprocess_t1(tmp_path / "t.nii.gz", tmp_path / "r.nii.gz",
                              tmp_path / "o", registration="magic")

    def test_wls_flag_passed_to_dtifit(self, tmp_path, monkeypatch):
        cap = _Capture(monkeypatch)
        shp = (3, 3, 3)
        d = tmp_path / "tensor"
        d.mkdir()
        for n, v in (("L1", 1.5), ("L2", .6), ("L3", .4)):
            _save(d / f"dti_{n}.nii.gz", np.full(shp, v))
        _save(tmp_path / "m.nii.gz", np.ones(shp))
        npl.estimate_tensor(tmp_path / "d", tmp_path / "bvec", tmp_path / "bval",
                            tmp_path / "m.nii.gz", d, wls=True)
        dt = next(c for c in cap.calls if c[0] == "dtifit")
        assert "-w" in dt and "--save_tensor" in dt

    def test_collect_versions(self):
        v = npl.collect_versions()
        assert v["nerves"] == npl.__version__ and "numpy" in v and "nibabel" in v


# ----------------------------------------------------------- group tools
class TestBidsAndSummary:

    def _bids(self, root):
        for sub, ses in (("01", None), ("02", "1")):
            base = root / f"sub-{sub}" / (f"ses-{ses}" if ses else "")
            (base / "dwi").mkdir(parents=True)
            (base / "anat").mkdir()
            (base / "fmap").mkdir()
            stem = f"sub-{sub}" + (f"_ses-{ses}" if ses else "")
            for ext in (".nii.gz", ".bval", ".bvec"):
                (base / "dwi" / f"{stem}_dwi{ext}").write_text("x")
            (base / "dwi" / f"{stem}_dwi.json").write_text(json.dumps(
                {"PhaseEncodingDirection": "j-", "TotalReadoutTime": 0.062}))
            (base / "anat" / f"{stem}_T1w.nii.gz").write_text("x")
            (base / "fmap" / f"{stem}_dir-PA_epi.nii.gz").write_text("x")
            (base / "fmap" / f"{stem}_dir-PA_epi.json").write_text(
                json.dumps({"PhaseEncodingDirection": "j"}))
        # participant 03: no bvec -> ignored
        (root / "sub-03" / "dwi").mkdir(parents=True)
        (root / "sub-03" / "dwi" / "sub-03_dwi.nii.gz").write_text("x")

    def test_opposite(self):
        assert nerves_bids.opposite("j-") == "j" and nerves_bids.opposite("i") == "i-"

    def test_discover_and_command(self, tmp_path):
        self._bids(tmp_path)
        jobs = nerves_bids.discover(tmp_path)
        assert [j["name"] for j in jobs] == ["sub-01", "sub-02_ses-1"]
        j = jobs[0]
        assert j["pe_dir"] == "j-" and j["readout_time"] == 0.062
        assert j["rpe_b0"].name.endswith("_epi.nii.gz") and j["t1"] is not None
        cmd = nerves_bids.build_command(j, tmp_path / "out", ["--wls"])
        assert "--readout_time" in cmd and "0.062" in cmd and "--rpe_b0" in cmd
        assert cmd[-1] == "--wls" and cmd[cmd.index("--pe_dir") + 1] == "j-"
        assert len(nerves_bids.discover(tmp_path, ["02"])) == 1

    def test_dry_run(self, tmp_path, capsys):
        self._bids(tmp_path)
        rc = nerves_bids.main([str(tmp_path), "--dry_run", "--wls"])
        out = capsys.readouterr().out
        assert rc == 0 and out.count("nde_pipeline.py") == 2 and "--wls" in out

    def test_group_summary(self, tmp_path):
        for n, v in (("a", 0.8), ("b", 0.9)):
            d = tmp_path / n
            d.mkdir()
            (d / "nde_pipeline_report.json").write_text(json.dumps({
                "timestamp": "t", "run_mode": "full",
                "nde_statistics": {"nde_mean": v, "n_mask_voxels": 10},
                "qc_results": {"csf_nde_pass": True},
                "snr": {"snr_median": 25.0}, "versions": {"nerves": "2.0.0-dev"}}))
        rows = gs.collect(tmp_path)
        assert [r["run"] for r in rows] == ["a", "b"]
        out = gs.write_csv(rows, tmp_path / "g.csv")
        table = list(csv.DictReader(open(out)))
        assert float(table[1]["nde_mean"]) == 0.9 and table[0]["snr_median"] == "25.0"
        assert gs.main([str(tmp_path)]) == 0


# ----------------------------------------------------------- validation
class TestValidationTools:

    def test_matched_fa_partner(self):
        base = [1.6, .45, .20]
        alt = ns.matched_fa_partner(base)
        fa = lambda e: nd.fa_md(np.array([e]))[0][0]  # noqa: E731
        assert fa(alt) == pytest.approx(fa(base), abs=1e-6)
        assert sum(alt) == pytest.approx(sum(base))
        assert alt[1] == pytest.approx(alt[2])

    def test_mode_beats_nde_for_lambda2_lambda3_changes(self):
        d = ns.run_discrimination(snr=25, n=4000, seed=5)
        a = d["A: same FA & trace, different λ2/λ3 split"]
        assert a["MODE"] > 5 * a["NDE"] and a["NDE"] > a["FA"]
        b = d["B: symmetric radial diffusivity increase"]
        assert b["NDE"] > 0.7 * b["FA"] and b["MODE"] < b["FA"]

    def test_nde_biased_low_at_low_snr(self):
        r = ns.run_precision(snrs=(8, 40), n=6000, seed=1)
        truth = r["truth"]["NDE"]
        low, high = r["rows"][0], r["rows"][1]
        assert truth - low["NDE"][0] > truth - high["NDE"][0] > -0.01
        assert low["invalid_pct"] > high["invalid_pct"]

    def test_icc(self):
        rng = np.random.default_rng(0)
        subj = rng.normal(0.8, 0.05, 40)
        good = np.c_[subj + rng.normal(0, .002, 40), subj + rng.normal(0, .002, 40)]
        bad = np.c_[subj + rng.normal(0, .05, 40), subj + rng.normal(0, .05, 40)]
        g, b = tr.icc_two_sessions(good), tr.icc_two_sessions(bad)
        assert g["ICC3_1"] > 0.95 > b["ICC3_1"] and g["wCV_pct"] < b["wCV_pct"]
        with pytest.raises(ValueError):
            tr.icc_two_sessions(good[:2])

    def test_retest_manifest(self, tmp_path):
        rng = np.random.default_rng(2)
        shp = (3, 3, 3)
        _save(tmp_path / "mask.nii.gz", np.ones(shp))
        rows = []
        for s in range(6):
            base = 0.7 + 0.03 * s
            for ses in ("1", "2"):
                _save(tmp_path / f"nde_{s}_{ses}.nii.gz",
                      base + rng.normal(0, .001, shp))
                rows.append({"subject": f"s{s}", "session": ses, "mask": "mask.nii.gz",
                             "nde": f"nde_{s}_{ses}.nii.gz"})
        with open(tmp_path / "m.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        res = tr.analyse(tmp_path / "m.csv")
        assert res["nde"]["n"] == 6 and res["nde"]["ICC3_1"] > 0.95


# ----------------------------------------------------------- integration
@pytest.mark.integration
@pytest.mark.skipif(shutil.which("dtifit") is None, reason="FSL not installed")
def test_dtifit_agrees_with_native_fit(phantom, tmp_path):
    """With FSL available: dtifit and the NumPy fit give the same NDE."""
    p, _ = phantom
    nat = nd.estimate_tensor_numpy(p["dwi"], p["bvec"], p["bval"], p["mask"],
                                   tmp_path / "nat", "ols")
    fsl = npl.estimate_tensor(p["dwi"], p["bvec"], p["bval"], p["mask"],
                              tmp_path / "fsl", method="dtifit")
    m = nib.load(str(p["mask"])).get_fdata() > 0
    for k in ("l1", "l2", "l3"):
        a = nib.load(str(nat[k])).get_fdata()[m]
        b = nib.load(str(fsl[k])).get_fdata()[m]
        np.testing.assert_allclose(a, b, rtol=0.05, atol=5e-5)


# ----------------------------------------------------------- full-mode wiring
class TestFullModeOrchestration:
    """Drive the real main() in full mode with the FSL/MRtrix stages faked,
    to verify argument threading and data flow between stages."""

    def test_full_mode_wiring(self, phantom, tmp_path, monkeypatch):
        p, d = phantom
        seen = {}
        shp = nib.load(str(p["mask"])).shape
        lab = nib.load(str(p["truth_labels"])).get_fdata()
        monkeypatch.setattr(npl, "validate_tools", lambda **kw: seen.update(tools=kw))

        def fake_dwi(**kw):
            seen["dwi"] = kw
            pre = kw["out_dir"]
            pre.mkdir(parents=True, exist_ok=True)
            b0 = nd.mean_b0_numpy(p["dwi"], p["bval"], pre / "mean_b0_eddy.nii.gz")
            _save(pre / "noise_map.nii.gz", np.full(shp, 1000.0 / 40))
            return {"dwi_eddy": p["dwi"], "bvecs": p["bvec"], "bvals": p["bval"],
                    "mask": p["mask"], "mean_b0": b0,
                    "mean_b0_brain": b0, "noise_map": pre / "noise_map.nii.gz"}

        def fake_t1(**kw):
            seen["t1"] = kw
            o = kw["out_dir"]
            o.mkdir(parents=True, exist_ok=True)
            res = {}
            for name, l in (("csf", 3), ("gm", 2), ("wm", 1)):
                _save(o / f"pve_{name}.nii.gz", (lab == l).astype(float))
                res[f"pve_{name}"] = o / f"pve_{name}.nii.gz"
            return res

        real_est = npl.estimate_tensor

        def fake_est(**kw):
            seen["est"] = kw
            kw = dict(kw, method="python")
            return real_est(**kw)

        monkeypatch.setattr(npl, "preprocess_dwi", fake_dwi)
        monkeypatch.setattr(npl, "preprocess_t1", fake_t1)
        monkeypatch.setattr(npl, "estimate_tensor", fake_est)
        t1 = tmp_path / "t1.nii.gz"
        _save(t1, np.ones(shp))
        out = tmp_path / "o"
        try:
            npl.main(["--dwi", str(p["dwi"]), "--bval", str(p["bval"]),
                      "--bvec", str(p["bvec"]), "--t1", str(t1), "--out", str(out),
                      "--pe_dir", "j-", "--readout_time", "0.062", "--resume",
                      "--n_threads", "3", "--wls", "--max_b", "1500",
                      "--min_snr", "10", "--registration", "epi_reg"])
        except SystemExit as e:
            assert e.code == 2
        # arguments reach the stages
        assert seen["dwi"]["pe_dir"] == "j-" and seen["dwi"]["readout_time"] == 0.062
        assert seen["dwi"]["resume"] and seen["dwi"]["n_threads"] == 3
        assert seen["t1"]["registration"] == "epi_reg" and seen["t1"]["resume"]
        assert seen["t1"]["dwi_ref_brain"] is not None
        assert seen["est"]["wls"] is True
        assert seen["tools"]["registration"] == "epi_reg"
        # data flow: SNR map computed from the noise map and applied
        assert (out / "preproc" / "snr_map.nii.gz").exists()
        rpt = json.loads((out / "nde_pipeline_report.json").read_text())
        assert rpt["run_mode"] == "full" and rpt["snr"]["snr_median"] > 20
        assert rpt["nde_statistics"]["n_low_snr"] == 0     # SNR ~40 > 10
        assert (out / "qc" / "qc_tissue_distributions.png").exists()
        assert (out / "shape" / "mode.nii.gz").exists()
