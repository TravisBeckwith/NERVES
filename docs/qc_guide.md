# QC Interpretation Guide

NERVES generates four QC figures and a JSON report for every run.
This guide explains what to look for in each.

> **Important.** The automated checks are coarse sanity screens. The reference ranges
> quoted below are **rules of thumb, not validated norms**; they depend on b-value, SNR,
> resolution, and how tissue masks are defined. Absolute white-matter NDE is
> tied to FA: for a prolate tensor with MD = 0.8×10⁻³ mm²/s, FA 0.35 / 0.45 / 0.55 /
> 0.70 corresponds to NDE ≈ 0.96 / 0.93 / 0.89 / 0.80. Judge NDE against your own FA
> distribution, not against fixed numbers.

---

## qc_histograms.png

**What it shows:** NDE, FA, and MD distributions within the white matter mask.

**What to look for:**

| Metric | Typical in healthy WM (rule of thumb) | Red flag |
|---|---|---|
| NDE | Unimodal; peak depends on FA (see above), commonly ~0.85–0.95 for whole-WM masks | Bimodal; very broad |
| FA | Unimodal, peak ~0.35–0.55 | Very broad distribution; peak < 0.25 |
| MD | Peak ~0.7–0.9 ×10⁻³ mm²/s | Peak > 1.2 (over-smoothing or CSF contamination) |

A **right-shifted NDE distribution** (clearly higher than expected for the observed FA) suggests either:
- CSF contamination of the WM mask
- Poor tensor estimation (low SNR, insufficient directions)
- Inadequate eddy current correction inflating λ₂ and λ₃

---

## qc_tissue_distributions.png

**What it shows:** KDE curves of NDE within CSF, grey matter,
GM–WM boundary, and white matter (each segmented at >90% PVE purity).

**Expected ordering (high to low NDE):**

```
CSF (highest, close to 1) > Grey matter > WM (lowest)
```

**Red flags:**
- WM NDE overlapping or exceeding GM NDE: tensor estimation failure or
  registration error placing CSF/GM voxels in the WM mask
- CSF peak < 0.90: unusually anisotropic free-water compartment, likely
  due to a very short TE or significant pulsation artefact

---

## qc_nde_fa_scatter.png

**What it shows:** NDE vs FA scatter for a random 50,000-voxel subsample
within white matter, coloured by FA value.

**Expected:** Strong negative correlation. This is a mathematical consequence of both
metrics being computed from the same eigenvalues, so the Pearson r is **informational
only and is not a pass/fail criterion**. It cannot detect wrong bvecs or failed eddy
correction (the tensor can be badly wrong and NDE and FA will still be anticorrelated).

**What to look for:**
- The vertical spread of NDE at a given FA reflects the λ₂/λ₃ split (tensor mode). A
  visibly different spread between subjects/sessions can indicate noise or
  tensor-fit differences, but is not diagnostic by itself.
- Use directional FA colour maps and eddy QC to assess bvecs and motion correction.

---

## qc_axial_montage.png

**What it shows:** Equally-spaced axial slices (up to six) displaying NDE (hot
colourmap), FA (greyscale), and MD (jet colourmap) side by side.

**What to check:**
- NDE and FA maps should show **inverted contrast**: bright FA = dark NDE
  and vice versa, particularly in the corpus callosum and corticospinal tracts
- Check for **ringing artefacts** (concentric bands near bright structures) —
  these indicate insufficient Gibbs removal
- Check for **geometric distortion** in the NDE map relative to the T1:
  the brain outline should align with the T1 brain mask
- **NDE in deep white matter approaching or exceeding cortical values** suggests
  CSF partial-voluming or a tensor-fit problem

---

## nde_pipeline_report.json

The JSON report records the sanity-check results. Key fields:

```json
{
  "nde_statistics": {
    "n_mask_voxels": 850000,
    "n_valid_voxels": 848200,
    "n_nonphysical": 1800,
    "nonphysical_pct": 0.212,
    "nde_mean": 0.883,
    "nde_std": 0.062,
    "nde_median": 0.891,
    "nde_p5": 0.771,
    "nde_p95": 0.967
  },
  "qc_results": {
    "csf_nde_mean": 0.974,
    "csf_nde_pass": true,
    "gm_nde_mean": 0.91,
    "wm_nde_mean": 0.84,
    "tissue_ordering_pass": true,
    "nde_fa_pearson_r": -0.951
  }
}
```

(Values are illustrative.) `nde_statistics` are computed over **valid** voxels only;
voxels with a non-positive eigenvalue, or with SNR below `--min_snr`, are NaN in the map
and excluded (`nde_qc_flag.nii.gz`: 1 = non-physical, 2 = low SNR; counts in
`n_nonphysical` / `n_low_snr`). The WM mean, `nde_fa_pearson_r` and
`nde_variance_explained_by_fa` (= r²) are informational; the latter shows how much of
NDE is simply FA. In full mode the report also records `snr` (median and 5th/95th
percentile b=0 SNR), `shape_statistics` (mode and Westin means), the acquisition's
b-value shells, and software `versions`.

A run passes the automated screens when `csf_nde_pass` and `tissue_ordering_pass` are
`true` and `nonphysical_pct` is below 5%. This is a necessary, not sufficient, condition
for good data.

---

## Common Issues and Fixes

| Symptom | Likely cause | Fix |
|---|---|---|
| `nonphysical_pct` > 5% | Poor tensor fit; low SNR; wrong bvecs | Check bvec sign convention; increase b=0 volumes; review SNR |
| `csf_nde_pass: false` (CSF NDE < 0.90) | T1 registration failure; FAST misclassification | Re-run with `--use_ants`; inspect FAST segmentation manually |
| `tissue_ordering_pass: false` | Registration or FAST failure; tensor fit problem | Inspect PVE overlays on the b0; inspect eddy QC files (`*.eddy_outlier_report`) |
| Gibbs ringing in NDE map | Insufficient `mrdegibbs` | Ensure raw DWI is passed (not already Gibbs-corrected); check partial Fourier settings |
