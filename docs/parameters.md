# Pipeline Parameters

## Required Arguments

| Argument | Type | Description |
|---|---|---|
| `--dwi` | path | 4D DWI NIfTI file (`.nii` or `.nii.gz`) |
| `--bval` | path | FSL-format b-value file |
| `--bvec` | path | FSL-format b-vector file |
| `--t1` | path | 3D T1-weighted NIfTI file |
| `--out` | path | Output directory (created if it does not exist) |

---

## Optional Arguments

### Acquisition Parameters

| Argument | Default | Description |
|---|---|---|
| `--pe_dir` | `AP` | Phase-encode direction: `AP`, `PA`, `LR`, or `RL` |
| `--rpe_b0` | None | Reverse phase-encode b=0 volume for topup distortion correction |
| `--no_topup` | False | Skip topup even when `--rpe_b0` is supplied |

> **Note on readout time:** The default readout time is set to 0.05 s in
> `nde_pipeline.py`. Update the `readout_time` variable near line 100 to
> match your acquisition. This value is scanner- and sequence-dependent
> and should be derived from your DICOM headers (typically
> `EffectiveEchoSpacing × (PhaseEncodingSteps − 1)`).

### Registration

| Argument | Default | Description |
|---|---|---|
| `--use_ants` | False | Use ANTs SyN nonlinear registration for T1→DWI (default: FSL FLIRT 6-DOF affine) |

ANTs registration is recommended when:
- T1 and DWI have significant field-of-view differences
- Subjects have pathology that distorts local anatomy
- High-resolution T1 data (≤ 0.8 mm isotropic) is available

### Skip Preprocessing Mode

Use this mode when DTI eigenvalue maps have already been computed (e.g.,
from a separate preprocessing pipeline or from FSL `dtifit`).

| Argument | Required if `--skip_preproc` | Description |
|---|---|---|
| `--skip_preproc` | — | Enable skip mode |
| `--l1` | Yes | λ₁ NIfTI map (axial diffusivity) |
| `--l2` | Yes | λ₂ NIfTI map |
| `--l3` | Yes | λ₃ NIfTI map |
| `--mask` | Yes | Brain mask NIfTI (binary) |

> When using `--skip_preproc`, the `--dwi`, `--bval`, `--bvec`, and `--t1`
> arguments are still required by the argument parser but are not used.
> Pass placeholder values (e.g., `x`) for these.

### Optional Pre-computed Tissue Maps

If you have PVE maps already registered to DWI space, you can supply them
directly to enable tissue-stratified QC figures:

| Argument | Description |
|---|---|
| `--pve_csf` | CSF partial volume estimate in DWI space |
| `--pve_gm` | Grey matter partial volume estimate in DWI space |
| `--pve_wm` | White matter partial volume estimate in DWI space |

If not supplied, the pipeline generates PVE maps from the T1 using FSL FAST
and registers them automatically.

### Logging

| Argument | Description |
|---|---|
| `--verbose` | Enable DEBUG-level logging (prints all subprocess stdout/stderr) |

---

## Exit Codes

| Code | Meaning |
|---|---|
| 0 | Success, all QC thresholds passed |
| 1 | Fatal error (missing inputs, tool not found, file I/O failure) |
| 2 | Success but non-physical eigenvalue rate exceeds 5% — review preprocessing |

---

## Environment Variables

NERVES respects standard FSL and MRtrix3 environment variables:

| Variable | Effect |
|---|---|
| `FSLDIR` | FSL installation directory |
| `FSLOUTPUTTYPE` | Output format (default: `NIFTI_GZ`) |
| `MRTRIX_NTHREADS` | Number of threads for MRtrix3 commands |

---

## Tuning Preprocessing Parameters

The following parameters are set as constants near the top of
`nde_pipeline.py` and may need adjustment for your data:

| Parameter | Location | Default | Notes |
|---|---|---|---|
| `readout_time` | `preprocess_dwi()` | `0.05` s | Update from DICOM headers |
| BET fractional threshold (DWI) | `preprocess_dwi()` | `0.3` | Lower for poor brain/background contrast |
| BET fractional threshold (T1) | `preprocess_t1()` | `0.35` | Raise if over-stripping; lower if under-stripping |
| FAST number of classes | `preprocess_t1()` | `3` | Do not change for standard brain segmentation |
| FLIRT DOF | `preprocess_t1()` | `6` (rigid) | Use `12` if T1 and DWI have scaling differences |
