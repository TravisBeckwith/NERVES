# Troubleshooting

## Preprocessing Issues

### Eddy fails immediately with "cannot read file"
- Verify the `--bval` and `--bvec` paths are correct and the files are
  not empty.
- Check that `--dwi` is a 4D NIfTI (not a single volume).
- Run `mrinfo your_dwi.nii.gz` to confirm shape and dimensionality.

### Topup fails with "must have even number of volumes"
- The `b0_pair.nii.gz` must contain exactly two volumes: one in the
  primary PE direction and one in the reverse PE direction.
- If your reverse PE b=0 file contains multiple volumes, extract the
  first with: `fslroi rpe_b0.nii.gz rpe_b0_single.nii.gz 0 1`

### Brain mask looks too small / too large
- Adjust the BET threshold in `preprocess_dwi()`:
  - Increase `-f` (e.g., `0.4`) to shrink the mask (remove more)
  - Decrease `-f` (e.g., `0.2`) to expand the mask (keep more)
- For subjects with large ventricles or significant atrophy, consider
  running BET with `-R` (robust mode) or using a different masking tool
  such as `hdbet` or `optiBET`.

### dwidenoise produces no noise map
- Ensure the input DWI has at least 4 volumes per slice direction.
- If your acquisition has very few volumes (< 20 total), MP-PCA may
  not have sufficient data. Consider skipping denoising for n < 20.

---

## Tensor Estimation Issues

### High non-physical eigenvalue rate (> 5%)
Non-physical (negative) eigenvalues arise when tensor fitting is
ill-conditioned. Common causes:

1. **Wrong bvec sign convention** — FSL uses a specific sign convention.
   If your scanner uses a different convention, the bvecs may need to be
   negated. Visualise FA directional maps in `fsleyes` to check.
2. **Insufficient gradient directions** — Fewer than 12 non-collinear
   directions produces unreliable tensors. At least 30 directions is
   recommended.
3. **Low SNR** — Use more b=0 averages or reduce the b-value.
4. **Eddy correction failure** — Review the eddy outlier report:
   `cat dwi_eddy.eddy_outlier_report`

### FA map looks inverted or incorrect
- Confirm bvec convention (see above).
- Verify b-value is correct: `cat your_dwi.bval` — all DWI volumes
  should have b > 0, and b=0 volumes should have b = 0.
- Run `fslview_deprecated dti_V1.nii.gz dti_FA.nii.gz` and verify
  that primary eigenvectors (colour-coded) match known anatomy
  (corpus callosum = red, corticospinal tract = blue,
  superior longitudinal fasciculus = green).

---

## NDE Computation Issues

### All NDE values are zero
- The brain mask is empty or incorrectly formatted. Verify with:
  `fslstats brain_mask.nii.gz -V` — the first number should be > 0.
- Eigenvalue maps may be in the wrong units. dtifit outputs eigenvalues
  in mm²/s; all values should be > 0 within the brain mask.

### NDE values outside [0, 1]
- Should not occur: NDE is clipped to [0, 1] internally. If you see this
  in a custom script, verify the `np.log(3)` denominator is non-zero and
  that the `np.clip` call is present.

---

## Registration Issues

### T1 PVE maps misaligned with NDE map
- This is the most common failure mode when using default FLIRT.
  Re-run with `--use_ants` for improved accuracy.
- Visually inspect alignment in `fsleyes`:
  ```bash
  fsleyes nde/nde.nii.gz t1/T1_pve_wm_DWIspace.nii.gz -cm green
  ```
- If the T1 and DWI have a large field-of-view mismatch, pre-crop both
  images before running NERVES.

### FAST produces incorrect tissue classes
- FAST may fail if the T1 image has very low contrast (e.g., IR-prepared
  sequences with unusual TI). In this case, supply pre-computed PVE maps
  via `--pve_csf`, `--pve_gm`, `--pve_wm` and skip FAST.
- Alternatively, run FAST manually with the `-t 2` flag (T2 input type)
  if your structural image is not a standard T1.

---

## QC Figure Issues

### qc_tissue_distributions.png is missing
- This figure requires `--pve_csf`, `--pve_gm`, and `--pve_wm` maps
  (either supplied or generated from the T1). If T1 preprocessing was
  skipped, this figure will not be produced.

### matplotlib errors / figures not saved
- Ensure a non-interactive backend is set. NERVES forces `matplotlib.use("Agg")`
  at startup, but if you import matplotlib before running NERVES in a script,
  the backend may already be set. Run NERVES as a standalone subprocess or
  ensure `Agg` is set before any other matplotlib import.

---

## Getting Help

If you encounter an issue not covered here, please open a GitHub Issue
and include:

1. The full command you ran
2. The contents of `nde_pipeline_report.json` (if it was generated)
3. The last 50 lines of terminal output (run with `--verbose` if possible)
4. Your FSL version (`flirt -version`), MRtrix3 version (`mrconvert --version`),
   and Python version (`python --version`)
