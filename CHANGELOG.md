# Changelog

All notable changes to NERVES will be documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
NERVES uses [Semantic Versioning](https://semver.org/).

---

## [Unreleased]

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
