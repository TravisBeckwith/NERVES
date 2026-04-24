# Changelog

All notable changes to NERVES will be documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
NERVES uses [Semantic Versioning](https://semver.org/).

---

## [Unreleased]

## [1.0.0] - 2025-XX-XX

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
