# Pipeline Parameters

NERVES has three run modes. Which inputs are required depends on the mode.

| Mode | Command line | Needs FSL/MRtrix3? | Required inputs |
|---|---|---|---|
| **Full** (default) | *(none)* | Yes | `--dwi --bval --bvec --t1 --out` |
| **Preprocessed DWI** | `--skip_dwi_preproc` | Only for `dtifit` / T1 steps | `--dwi --bval --bvec --out` + `--mask` (or FSL `bet`) |
| **Eigenvalue maps** | `--skip_preproc` | No | `--l1 --l2 --l3 --mask --out` |

With `--skip_dwi_preproc --tensor_fit python --mask ...` and no `--t1`, no external
tools are needed at all.

---

## Inputs

| Argument | Description |
|---|---|
| `--dwi` | 4D DWI NIfTI. With `--skip_dwi_preproc` this must already be denoised and corrected for distortion, eddy currents and motion, with matching final bvecs |
| `--bval`, `--bvec` | FSL-format gradient files. Validated before any processing (length vs. volumes, ≥ 1 b=0, ≥ 6 directions, unit-length vectors) |
| `--t1` | 3D T1-weighted NIfTI. Required for full preprocessing; optional with `--skip_dwi_preproc` (enables tissue-specific QC) |
| `--out` | Output directory (created if absent) |
| `--mask` | Brain mask (`--skip_preproc` and `--skip_dwi_preproc`) |

---

## Acquisition and preprocessing

| Argument | Default | Description |
|---|---|---|
| `--pe_dir` | `AP` | Phase-encode direction: `AP`, `PA`, `LR`, `RL`, or the BIDS axis code `i`, `i-`, `j`, `j-`, `k`, `k-` (`j-` ≡ `AP`) |
| `--rpe_b0` | None | Reverse phase-encode b=0 volume(s) for topup. 4-D inputs are averaged |
| `--no_topup` | False | Skip topup even when `--rpe_b0` is supplied |
| `--readout_time` | `0.05` | Total readout time (s) for topup/eddy `acqparams.txt` (BIDS: `TotalReadoutTime`). The default is a placeholder, not a measured value — set it for real data |
| `--resume` | False | Skip preprocessing stages whose outputs already exist in `--out` (denoise, de-Gibbs, topup, eddy, masks, T1 steps, tensor fit). Delete a stage's outputs to force a rerun |
| `--n_threads` | 0 | Threads for MRtrix3 tools and eddy (`--nthr`). 0 = tool default |
| `--eddy_gpu` | False | Prefer a CUDA build of eddy if installed. The eddy executable is detected automatically (`eddy_cpu`, `eddy_openmp`, `eddy`) |

### Registration (T1 → DWI)

| Argument | Description |
|---|---|
| `--registration flirt` | (default) FSL FLIRT, 6 DOF, normalised mutual information |
| `--registration epi_reg` | FSL boundary-based registration using the WM segmentation. Generally more accurate for EPI data. Needs FSL `epi_reg` and `convert_xfm` |
| `--registration ants` | ANTs rigid + affine (no nonlinear stage is run) |
| `--use_ants` | Alias for `--registration ants` |

---

## Tensor estimation

| Argument | Default | Description |
|---|---|---|
| `--tensor_fit` | `dtifit` | `dtifit` (FSL) or `python` (native NumPy fit; no FSL needed) |
| `--wls` | False | Weighted least squares instead of ordinary least squares. In simulation (SNR 8–10) it cut the fraction of negative-eigenvalue voxels 2–5× and reduced NDE variance and low-SNR bias (`python validation/noise_simulation.py --method wls`). Recommended for noisy data |
| `--max_b` | None | Fit only volumes with b ≤ this value (b=0 volumes kept). Use for multi-shell data or to drop high-b shells: DTI assumes Gaussian diffusion, and above roughly b = 1000–1500 s/mm² kurtosis biases the fit. A warning is issued for higher b-values |

---

## Quality control and outputs

| Argument | Default | Description |
|---|---|---|
| `--min_snr` | 0 (off) | Set NDE to NaN, with QC flag value 2, where voxel SNR < this value. NDE is biased low at low SNR. In full mode SNR = mean b=0 / `dwidenoise` noise sigma (approximate: sigma is in acquisition space); otherwise supply `--snr` |
| `--snr` | None | Optional SNR map (NIfTI, same grid) for `--min_snr` |
| `--nde_only` | False | Skip the companion shape maps (tensor mode, Westin cl/cp/cs) |
| `--fa`, `--md` | auto | FA/MD maps for QC plots in `--skip_preproc` mode. If omitted, NERVES looks next to `--l1` for `..._FA` / `..._MD` (only if the filename contains `L1`) |
| `--pve_csf`, `--pve_gm`, `--pve_wm` | None | Partial-volume maps already in DWI space, enabling tissue-stratified QC |
| `--verbose` | False | DEBUG logging (subprocess stdout/stderr) |
| `--version` | | Print the version |

### `--skip_preproc` inputs

| Argument | Required | Description |
|---|---|---|
| `--l1`, `--l2`, `--l3` | Yes | Eigenvalue maps (any units; NDE is scale-invariant) |
| `--mask` | Yes | Binary brain mask |

`--dwi`, `--bval`, `--bvec`, `--t1` are not needed in this mode (older command lines that
pass placeholder values such as `x` still work).

---

## Outputs added by the optional features

| File | When |
|---|---|
| `shape/mode.nii.gz`, `cl`, `cp`, `cs` | unless `--nde_only`. Tensor mode (−1 planar … +1 linear; NaN for isotropic tensors) and Westin linear/planar/spherical coefficients, defined on exactly the same voxels as the NDE map |
| `preproc/snr_map.nii.gz` | full mode |
| `preproc/shells/` | when `--max_b` removes volumes |
| `nde_pipeline_report.json` → `versions`, `run_mode`, `acquisition`, `snr`, `shape_statistics`, `tensor_fit` | always (fields present when applicable) |

---

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Fatal error (missing/invalid inputs, tool not found, shape mismatch) |
| 2 | Success, but the non-physical eigenvalue rate exceeds 5% — review preprocessing (consider `--wls`) |

---

## Environment variables

| Variable | Effect |
|---|---|
| `FSLDIR` | FSL installation directory (also used to record the FSL version in the report) |
| `FSLOUTPUTTYPE` | Output format (default: `NIFTI_GZ`) |
| `MRTRIX_NTHREADS` | Default MRtrix3 threads (overridden by `--n_threads`) |

---

## Tuning constants

Set in `nde_pipeline.py` and possibly worth adjusting:

| Parameter | Location | Default | Notes |
|---|---|---|---|
| BET fractional threshold (DWI) | `preprocess_dwi()` | `0.3` | Lower for poor brain/background contrast |
| BET fractional threshold (T1) | `preprocess_t1()` | `0.35` | Raise if over-stripping; lower if under-stripping |
| FAST classes | `preprocess_t1()` | `3` | Do not change |
| FLIRT DOF | `preprocess_t1()` | `6` | Use `12` if T1 and DWI differ in scaling |
| `MAX_RECOMMENDED_B` | `nerves_dti.py` | `1500` | b-value above which a warning is issued |

---

## Companion tools

| Script | Purpose |
|---|---|
| `nerves_bids.py` | Run the pipeline over a BIDS dataset (reads PE direction and readout time from sidecars; finds the T1w and an opposite-PE fieldmap EPI). Extra options are passed through |
| `nerves_group_summary.py` | Collect every `nde_pipeline_report.json` below a directory into one CSV |
| `validation/noise_simulation.py` | Monte-Carlo comparison of NDE with FA/MD/mode versus SNR and for three change scenarios |
| `validation/make_phantom.py` | Synthetic DWI phantom with known tensors |
| `validation/retest_reliability.py` | ICC / within-subject CV for test-retest data |
