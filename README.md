<div align="center">

<img src="docs/assets/nerves_banner.svg" width="800" alt="NERVES banner"/>

# NERVES
### Normalized Entropy Representation of Voxelwise Eigenvalue Spectra

[![Python](https://img.shields.io/badge/python-3.8%2B-blue?logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.19738588-orange)](https://doi.org/10.5281/zenodo.19738588)
[![FSL](https://img.shields.io/badge/requires-FSL%206.0%2B-blueviolet)](https://fsl.fmrib.ox.ac.uk/)
[![MRtrix3](https://img.shields.io/badge/requires-MRtrix3%203.0%2B-blueviolet)](https://www.mrtrix.org/)

**A microstructural diffusion MRI pipeline for computing Normalized Diffusion Entropy (NDE) maps from T1 and DWI data.**

[Installation](#installation) · [Quick Start](#quick-start) · [Documentation](#documentation) · [Citation](#citation)

</div>

---

## Overview

NERVES provides an end-to-end pipeline for computing **Normalized Diffusion Entropy (NDE)** maps from raw diffusion-weighted MRI data. NDE is the Shannon entropy of the (normalised) diffusion tensor eigenvalue spectrum.

### What NDE is — and is not

> **Scope note.** NDE is a *scalar summary of the three tensor eigenvalues*. It is
> **not** the model-free (Q-ball / mean-free-path) diffusion entropy of
> Fozouni et al. 2013, which is cited as related work on entropy-based diffusion
> contrast but uses a different estimator and sign convention.

Because NDE depends only on the eigenvalue *shape*, it is a deterministic function of
two orthogonal tensor invariants, FA and tensor mode (Ennis & Kindlmann 2006). In
simulation (see `tests/test_complementarity.py`) NDE correlates with FA at about
r ≈ −0.94, FA alone explains ≈ 98 % of NDE variance, and FA + mode explain > 99.9 %.
Practically:

- At **fixed FA**, NDE still varies with how the minor eigenvalues split (λ₂ vs λ₃),
  i.e. with tensor mode. That is the only information NDE adds to FA.
- Tensor mode captures that same information directly and has an established
  interpretation. If you need the λ₂/λ₃ contrast, report mode alongside NDE.
- Claims that NDE is more sensitive than FA to Wallerian degeneration, crossing
  fibres, etc. are **hypotheses, not results of this repository**: no clinical or
  phantom validation data are included. Treat NDE as an exploratory metric.
- Like FA, NDE is biased by noise (eigenvalue sorting inflates λ₁ and deflates λ₃),
  so low-SNR acquisitions will shift values.

### Mathematical Definition

For eigenvalues λ₁ ≥ λ₂ ≥ λ₃ > 0 at each voxel, define normalised proportions:

$$p_i = \frac{\lambda_i}{\lambda_1 + \lambda_2 + \lambda_3}$$

Shannon entropy of the distribution:

$$H = -\sum_{i=1}^{3} p_i \ln(p_i)$$

Normalised to [0, 1]:

$$\text{NDE} = \frac{H}{\ln(3)}$$

- **NDE = 1** → perfectly isotropic (λ₁ = λ₂ = λ₃)  
- **NDE → 0** → maximally anisotropic (all diffusion along one eigenvector)
- **NDE = NaN** → undefined: at least one eigenvalue ≤ 0 (flagged in `nde_qc_flag.nii.gz`)

---

## Pipeline Overview

```
Raw DWI + T1
     │
     ├── DWI Preprocessing
     │     ├── MP-PCA denoising       (MRtrix3: dwidenoise)
     │     ├── Gibbs ringing removal  (MRtrix3: mrdegibbs)
     │     ├── Brain masking          (FSL: bet)
     │     ├── Susceptibility corr.   (FSL: topup, optional)
     │     └── Eddy/motion correction (FSL: eddy)
     │
     ├── T1 Preprocessing
     │     ├── Brain extraction       (FSL: bet)
     │     ├── Tissue segmentation    (FSL: fast)
     │     └── Registration → DWI    (FSL: flirt 6-DOF, or ANTs rigid+affine)
     │
     ├── Tensor Estimation            (FSL: dtifit)
     │
     ├── NDE Computation              (Python/NumPy — native)
     │
     └── QC + Report
           ├── Metric histograms
           ├── Tissue-specific distributions
           ├── NDE vs FA scatter
           ├── Axial slice montage
           └── JSON summary report
```

---

## Installation

### Prerequisites

| Dependency | Version | Purpose |
|---|---|---|
| [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation) | ≥ 6.0 | Preprocessing, tensor estimation |
| [MRtrix3](https://www.mrtrix.org/download/) | ≥ 3.0.3 | Denoising, Gibbs removal |
| Python | ≥ 3.8 | NDE computation, QC |
| [ANTs](https://github.com/ANTsX/ANTs) | ≥ 2.3 | Optional: improved T1→DWI registration |

### Python Environment

```bash
# Clone the repository
git clone https://github.com/TravisBeckwith/NERVES.git
cd NERVES

# Create and activate a virtual environment (recommended)
python -m venv nerves_env
source nerves_env/bin/activate        # Linux/macOS
# nerves_env\Scripts\activate         # Windows

# Install Python dependencies
pip install -r requirements.txt
```

### Verify Installation

```bash
python nde_pipeline.py --help
```

You should see the full argument list. If FSL or MRtrix3 tools are missing from PATH, an informative error will be raised before any processing begins.

---

## Quick Start

### Full Pipeline (Raw Data)

```bash
python nde_pipeline.py \
    --dwi  sub-01_dwi.nii.gz \
    --bval sub-01_dwi.bval   \
    --bvec sub-01_dwi.bvec   \
    --t1   sub-01_T1w.nii.gz \
    --out  ./sub-01_nerves   \
    --pe_dir AP
```

### With Susceptibility Distortion Correction

```bash
python nde_pipeline.py \
    --dwi    sub-01_dwi.nii.gz \
    --bval   sub-01_dwi.bval   \
    --bvec   sub-01_dwi.bvec   \
    --t1     sub-01_T1w.nii.gz \
    --rpe_b0 sub-01_b0_PA.nii.gz \
    --pe_dir AP \
    --out    ./sub-01_nerves
```

### Skip Preprocessing (Pre-computed Eigenvalues)

If you already have DTI eigenvalue maps (e.g., from FSL `dtifit`):

```bash
python nde_pipeline.py \
    --skip_preproc \
    --l1   dti_L1.nii.gz \
    --l2   dti_L2.nii.gz \
    --l3   dti_L3.nii.gz \
    --mask brain_mask.nii.gz \
    --dwi x --bval x --bvec x --t1 x \
    --out  ./nerves_output
```

### ANTs Registration (rigid + affine)

```bash
python nde_pipeline.py \
    --dwi  sub-01_dwi.nii.gz \
    --bval sub-01_dwi.bval   \
    --bvec sub-01_dwi.bvec   \
    --t1   sub-01_T1w.nii.gz \
    --out  ./sub-01_nerves   \
    --use_ants
```

---

## Output Structure

```
sub-01_nerves/
├── preproc/
│   ├── dwi_denoised.mif          # MP-PCA denoised DWI
│   ├── dwi_degibbs.mif           # Gibbs-corrected DWI
│   ├── dwi_eddy.nii.gz           # Eddy-corrected DWI
│   ├── mean_b0_final_brain.nii.gz
│   └── mean_b0_final_brain_mask.nii.gz
├── t1/
│   ├── T1_brain.nii.gz
│   ├── T1_fast_pve_0.nii.gz      # CSF partial volume
│   ├── T1_fast_pve_1.nii.gz      # GM partial volume
│   ├── T1_fast_pve_2.nii.gz      # WM partial volume
│   └── T1_pve_*_DWIspace.nii.gz  # PVEs registered to DWI space
├── tensor/
│   ├── dti_L1.nii.gz             # λ₁ (axial diffusivity)
│   ├── dti_L2.nii.gz             # λ₂
│   ├── dti_L3.nii.gz             # λ₃
│   ├── dti_FA.nii.gz             # Fractional anisotropy
│   ├── dti_MD.nii.gz             # Mean diffusivity
│   └── dti_RD.nii.gz             # Radial diffusivity = (λ₂+λ₃)/2
├── nde/
│   ├── nde.nii.gz                # ★ PRIMARY OUTPUT: NDE map (NaN where undefined)
│   └── nde_qc_flag.nii.gz        # Voxels with non-physical eigenvalues
├── qc/
│   ├── qc_histograms.png         # NDE / FA / MD distributions in WM
│   ├── qc_tissue_distributions.png  # NDE by tissue class
│   ├── qc_nde_fa_scatter.png     # NDE vs FA scatter plot
│   └── qc_axial_montage.png      # Axial slice comparison
└── nde_pipeline_report.json      # Full QC statistics and thresholds
```

---

## QC Thresholds

NERVES applies the following automated sanity checks and records them in the JSON report.
They are coarse screens, **not** validated quality criteria; always inspect the QC figures.

| Metric | Expected | Interpretation if failed |
|---|---|---|
| Non-physical eigenvalue rate | < 5% (exit code 2 otherwise) | Tensor fitting or masking issue |
| CSF mean NDE (`csf_nde_pass`) | ≥ 0.90 | Registration / segmentation / tensor fit problem |
| Tissue ordering (`tissue_ordering_pass`) | mean NDE CSF > GM > WM | Registration / segmentation / tensor fit problem |

The NDE–FA Pearson r and absolute WM mean NDE are **reported for information only**.
NDE and FA are computed from the same eigenvalues, so they stay strongly anticorrelated
even for badly corrupted tensors, and absolute WM NDE depends on b-value, SNR,
resolution and how WM is defined.

---

## Documentation

- [Installation Guide](docs/installation.md)
- [Pipeline Parameters](docs/parameters.md)
- [NDE Theory](docs/theory.md)
- [QC Interpretation](docs/qc_guide.md)
- [Troubleshooting](docs/troubleshooting.md)
- [Contributing](CONTRIBUTING.md)

---

## Citation

If you use NERVES in your research, please cite:

```bibtex
@software{beckwith2026nerves,
  title   = {{NERVES}: Normalized Entropy Representation of Voxelwise
             Eigenvalue Spectra},
  author  = {Beckwith, Travis J.},
  year    = {2026},
  url     = {https://github.com/TravisBeckwith/NERVES},
  doi     = {10.5281/zenodo.19738588}
}
```

Related work on entropy-based diffusion contrast (a different, model-free estimator —
cite it for context, not as the source of this implementation):

```bibtex
@article{fozouni2013,
  title   = {Characterizing brain structures and remodeling after {TBI}
             based on information content, diffusion entropy},
  author  = {Fozouni, N. and Chopp, M. and Nejad-Davarani, S.P. and
             Zhang, Z.G. and Lehman, N.L. and Gu, S. and Jiang, Q.},
  journal = {PLoS One},
  volume  = {8},
  number  = {10},
  pages   = {e76343},
  year    = {2013},
  doi     = {10.1371/journal.pone.0076343}
}
```

Closely related tensor-based scalar measures, which you may also wish to cite:
Özarslan E, Vemuri BC, Mareci TH. *Magn Reson Med* 2005;53(4):866–876
(trace, variance and entropy measures); Ennis DB, Kindlmann G. *Magn Reson Med*
2006;55(1):136–146 (orthogonal invariants, tensor mode).

---

## License

NERVES is released under the [MIT License](LICENSE).

---

## Acknowledgements

NERVES wraps and depends on [FSL](https://fsl.fmrib.ox.ac.uk/) and [MRtrix3](https://www.mrtrix.org/). Please ensure you comply with their respective license terms.
