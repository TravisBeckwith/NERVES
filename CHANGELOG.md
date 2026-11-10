# Changelog

All notable changes to NERVES will be documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
NERVES uses [Semantic Versioning](https://semver.org/).

---

## [Unreleased]

### Added (v2.0.0-dev)
- **`nerves_qc.py`**: tensor-fit health and scalar-reliability QC for QSIRecon
  (DSI Studio / TORTOISE / DIPY tensor scalars) and NERVES/FSL `dti_*` outputs.
  Non-physical voxel, FA-saturation, MD-plausibility and MD/AD/RD consistency checks;
  expected FA/MD/mode/NDE error at the subject's SNR; robust-z group outlier flags;
  TSV, JSON and HTML output; exit code 3 on failure.
- **`nerves_reliability.py`**: Monte-Carlo bias/SD calculator for FA, MD, tensor mode
  and NDE versus SNR, directions, b-value and tissue; `--required METRIC:TOL` finds the
  minimum SNR for a target accuracy.
- `docs/reliability.md`, `tests/test_qc.py` (29 tests).

### Changed (v2.0.0-dev)
- Project reframed around tensor-fit QC and reliability; NDE remains available as an
  optional, exploratory metric. Pipeline version is now 2.0.0-dev.

### Added (v1.1.0-dev)
- **Companion shape maps** written by default: tensor mode and Westin cl/cp/cs
  (`shape/`), defined on the same voxels as NDE (`--nde_only` disables). Motivated by
  the new simulation showing mode is far more sensitive than NDE to λ₂/λ₃ changes.
- **Native NumPy tensor fit** (`--tensor_fit python`, OLS or `--wls`) in the new
  `nerves_dti.py`; `--wls` also maps to `dtifit -w`. In simulation WLS cuts
  negative-eigenvalue voxels 2–5× and reduces NDE variance and bias at low SNR.
- **Input validation** of the gradient table (length, b=0, directions, unit vectors),
  b-value range / multi-shell warnings and `--max_b` single-shell selection.
- **SNR map** from the `dwidenoise` noise sigma, `--min_snr` masking (QC flag 2),
  `--snr` for external maps, SNR statistics in the report.
- **`--skip_dwi_preproc`** entry point for already-preprocessed DWI; `--dwi/--bval/
  --bvec/--t1` are now optional where not needed.
- **Eddy executable detection** (`eddy_cpu`, `eddy_openmp`, `eddy`, CUDA with
  `--eddy_gpu`), `--n_threads`, `--resume`, `--registration {flirt,epi_reg,ants}`,
  BIDS phase-encode codes for `--pe_dir`, `--version`.
- Reproducibility: software versions (NERVES, Python libs, FSL, MRtrix3, git commit),
  run mode and acquisition summary are recorded in the report.
- `nerves_bids.py` (run a BIDS dataset) and `nerves_group_summary.py` (reports → CSV).
- `validation/` tools: Monte-Carlo noise simulation, synthetic phantom, test-retest
  ICC; `docs/testing.md`.
- ~50 new tests, including a full CLI run on a phantom with known tensors and
  command-line checks for every external tool call.

### Fixed (v1.1.0-dev)
- Reverse phase-encode b0 volumes with more than one volume produced a mismatch
  between `acqparams.txt` rows and the topup input; they are now averaged.
- Eigenvalue ordering check now also tests L2 ≥ L3 and cross-checks trace against
  dtifit's tensor file.
- Single shared NDE implementation (`nerves_dti.nde_from_eigenvalues`).

### Changed (behaviour — may alter results)
- **Non-physical voxels (any eigenvalue ≤ 0) are now `NaN` in `nde.nii.gz`** instead
  of `0`, and are excluded from all summary statistics and QC. Previously they were
  written as 0 ("maximally anisotropic") *and* counted in `nde_mean/min/percentiles`,
  the WM/CSF checks and the NDE–FA correlation, biasing all of them. Downstream code
  must handle NaN (e.g. `np.nanmean`, `fslmaths -nan`).
- **QC thresholds revised.** The NDE–FA `r < −0.80` check and the absolute
  `WM NDE < 0.85` check were removed as pass/fail criteria (`nde_fa_anticorrelated`
  and `wm_nde_pass` no longer exist in the report): NDE and FA are functions of the
  same eigenvalues so the correlation cannot reveal a corrupted tensor, and the WM
  threshold was inconsistent with the documented FA range. Added
  `tissue_ordering_pass` (mean NDE CSF > GM > WM); `nde_fa_pearson_r`, `wm_nde_mean`
  and `gm_nde_mean` are informational.
- External commands now run **without a shell** (argument lists); paths with spaces
  or shell metacharacters are safe, and failures inside former pipelines
  (`dwiextract | mrmath`) are no longer masked.
- The eddy brain mask is now derived from the topup-corrected b0 (via `applytopup`)
  when topup is used.
- ANTs option is now correctly described as **rigid + affine** (it never ran a
  nonlinear SyN stage).

### Fixed
- Shape mismatches between `--skip_preproc` inputs now exit with a clear error.
- Removed unused `mask_ref` argument of `preprocess_t1` and the unused `seaborn` /
  `pandas` requirements.

### Documentation
- Clarified that NDE is a **tensor-eigenvalue entropy**, not the model-free Q-ball
  diffusion entropy of Fozouni et al. 2013 (now cited as related work only).
- Removed unsupported claims of clinical sensitivity and of "complementarity" with FA;
  added a quantified analysis (NDE is ≈ a function of FA and tensor mode) and
  `tests/test_complementarity.py`.
- Corrected the Table 1 narrative (A–B and C–D differ in FA too; |ΔNDE| > |ΔFA| holds
  only for A–B). QC reference ranges are labelled rules of thumb.
- Added `.github/workflows/ci.yml` to the repository.

### Added
- `--readout_time` CLI argument — total readout time for topup/eddy
  acqparams.txt was previously hardcoded (0.05s) and required editing the
  script directly; now configurable per-run
- `--fa` / `--md` CLI arguments for `--skip_preproc` mode — explicit
  overrides for QC plotting when eigenvalue files don't follow the
  dtifit `..._L1.nii.gz` naming convention
- GitHub Actions CI workflow (`.github/workflows/ci.yml`) actually added —
  the 1.0.0 entry below referenced CI that didn't exist in the repository
- Direct unit tests against `nde_pipeline.compute_nde` (previously the
  suite only exercised a hand-mirrored copy of the math in the test file
  itself, not the production function)
- Unit tests for `validate_tools`, `write_report`, and CLI argument
  validation (`--skip_preproc` without required inputs)
- A full end-to-end test running the real CLI via `--skip_preproc`
  (`tests/test_pipeline.py::TestSkipPreprocEndToEnd`) — this is what
  surfaced both bugs below; the previous suite only tested the entropy
  math in isolation and never ran the pipeline as a whole

### Fixed
- **`write_report()` crashed with `TypeError: Object of type bool is not
  JSON serializable` on essentially every real run that reached QC.**
  `nde_fa_anticorrelated = r < -0.80` compared scipy's `pearsonr` result
  (a numpy float64) directly, producing a `numpy.bool_` rather than a
  Python `bool`, which `json.dump()` rejects. Fixed the comparison and
  added a numpy-safe JSON encoder to `write_report()` as a general
  safeguard against the same failure mode elsewhere.
- **FA/MD auto-detection in `--skip_preproc` mode silently used the L1
  eigenvalue map as if it were the FA map** whenever the input filename
  didn't contain the literal substring `"L1"` (i.e. didn't follow dtifit
  naming). `l1_path.name.replace("L1", "FA")` no-ops on such filenames,
  so the "candidate FA path" ended up identical to the L1 path — which of
  course exists — so it was accepted as FA data instead of triggering the
  intended "FA map not found, skipping QC" fallback. QC figures and
  thresholds would have been computed against eigenvalue magnitudes
  mislabeled as FA. Now guarded, and `--fa`/`--md` provide an explicit
  escape hatch.
- `docs/theory.md` Table 1: Config D FA was listed as 0.490; the correct
  value from the same formula used for Configs A–C is 0.463
- `CHANGELOG.md` 1.0.0 entry had an unfilled placeholder date (2025-XX-XX)

## [1.0.0] - 2026-04-24

### Added
- Initial public release
- Full DWI preprocessing pipeline (MP-PCA denoising, Gibbs removal, topup, eddy)
- T1 preprocessing with FSL FAST tissue segmentation
- T1→DWI registration via FSL FLIRT (affine) or ANTs SyN (nonlinear)
- Diffusion tensor estimation via FSL dtifit with eigenvalue ordering verification
- Vectorised NDE computation with edge-case handling (non-physical eigenvalues,
  p=0 continuity, scale invariance)
- Four QC figures: histograms, tissue distributions, NDE–FA scatter, axial montage
- JSON summary report with automated pass/fail QC thresholds
- `--skip_preproc` mode for use with pre-computed eigenvalue maps
- `--use_ants` flag for improved T1→DWI registration
- 14-test unit suite covering boundary conditions and manuscript Table 1 values
- GitHub Actions CI across Python 3.8, 3.10, 3.12
