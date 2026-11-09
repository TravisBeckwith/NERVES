# Testing & Validation

## Running the tests

```bash
pip install -r requirements.txt pytest
python -m pytest -m "not integration" -q        # no FSL / MRtrix3 required
```

The default suite (≈ 124 tests) covers:

* NDE maths and edge cases, including NaN handling and scale invariance
* gradient-table validation, shell selection, native tensor fitting (OLS/WLS)
* tensor mode and Westin indices, SNR maps, trace / ordering consistency checks
* the **full CLI end to end** on a synthetic phantom with known tensors
  (`--skip_dwi_preproc --tensor_fit python`), including recovery of WM NDE
* the exact external commands built for denoise/de-Gibbs/topup/eddy/dtifit/epi_reg
  (captured, not executed), including `--resume`, threading and BIDS PE codes
* the BIDS runner, group summary and noise simulation
* `nerves_reliability.py` (monotonicity in SNR/directions, known upward low-SNR FA bias,
  mode unreliability in low-FA tissue) and `nerves_qc.py` (QSIRecon and FSL layouts,
  unit detection, corrupted/empty/inconsistent maps, outlier flagging, exit codes, and
  a run on real pipeline output from the phantom)

## What the unit tests cannot cover

The FSL and MRtrix3 stages are only tested through **captured command lines**, not by
running the tools. `tests/test_features.py::test_dtifit_agrees_with_native_fit` is marked
`integration` and runs automatically when `dtifit` is on `PATH`:

```bash
python -m pytest -m integration -q
```

Please run it, and a full pipeline run on a real dataset, on a machine with FSL ≥ 6.0 and
MRtrix3 ≥ 3.0.3 before relying on the preprocessing stages. A container image bundling both
tools is the recommended way to make this repeatable.

## Synthetic phantom

```bash
python validation/make_phantom.py --out phantom/
python nde_pipeline.py --skip_dwi_preproc --tensor_fit python --wls \
    --dwi phantom/dwi.nii.gz --bval phantom/dwi.bval --bvec phantom/dwi.bvec \
    --mask phantom/mask.nii.gz --out phantom_out/
```

`phantom/truth_L1/L2/L3.nii.gz` hold the ground-truth eigenvalues; `truth_labels.nii.gz`
marks WM (1), GM (2) and CSF (3). Note that CSF at b = 1000 is close to the Rician noise
floor, so its eigenvalues are biased even though it stays clearly isotropic.

## Noise simulation

```bash
python validation/noise_simulation.py --snr 25 --n 20000 [--method wls]
```

Prints bias/precision versus SNR and d′ discriminability for three scenarios. See
`docs/theory.md` for interpretation.

## Test-retest reliability

Build a CSV manifest with one row per subject × session:

```
subject,session,mask,nde,fa,mode
sub-01,1,sub-01_wm.nii.gz,s01_1/nde/nde.nii.gz,s01_1/tensor/dti_FA.nii.gz,s01_1/shape/mode.nii.gz
sub-01,2,sub-01_wm.nii.gz,s01_2/nde/nde.nii.gz,...
```

```bash
python validation/retest_reliability.py manifest.csv --out reliability.csv
```

reporting ICC(3,1), ICC(2,1), within-subject CV and between-subject SD for every metric
column. Use it to check whether NDE is actually more reproducible than FA or mode in your
data before using it in a study.
