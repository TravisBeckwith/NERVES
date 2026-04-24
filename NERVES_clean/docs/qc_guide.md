# QC Interpretation Guide

NERVES generates four QC figures and a JSON report for every run.
This guide explains what to look for in each.

---

## qc_histograms.png

**What it shows:** NDE, FA, and MD distributions within the white matter mask.

**What to look for:**

| Metric | Expected in healthy WM | Red flag |
|---|---|---|
| NDE | Unimodal, peak ~0.78–0.88 | Bimodal; peak > 0.92 or < 0.70 |
| FA | Unimodal, peak ~0.35–0.55 | Very broad distribution; peak < 0.25 |
| MD | Peak ~0.7–0.9 ×10⁻³ mm²/s | Peak > 1.2 (over-smoothing or CSF contamination) |

A **right-shifted NDE distribution** (mean > 0.90 in WM) suggests either:
- CSF contamination of the WM mask
- Poor tensor estimation (low SNR, insufficient directions)
- Inadequate eddy current correction inflating λ₂ and λ₃

---

## qc_tissue_distributions.png

**What it shows:** KDE curves of NDE within CSF, grey matter,
GM–WM boundary, and white matter (each segmented at >90% PVE purity).

**Expected ordering (high to low NDE):**

```
CSF (≈ 0.97–1.00) > Grey matter (≈ 0.88–0.94) > WM (≈ 0.78–0.88)
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

**Expected:** Strong negative correlation (Pearson r ≈ −0.92 to −0.98).

**What to look for:**
- **r < −0.80:** Acceptable — passes QC threshold
- **r > −0.80:** Unexpected; usually indicates a preprocessing failure
  that has corrupted the tensor (e.g., wrong bvecs, failed eddy run)
- **Yellow highlighted zone (FA 0.3–0.6):** This is the region of maximum
  NDE–FA dissociation, where the two metrics diverge most. A healthy dataset
  shows the widest vertical spread of NDE values in this zone. Narrow spread
  here suggests insufficient sensitivity differentiation.

---

## qc_axial_montage.png

**What it shows:** Six equally-spaced axial slices displaying NDE (hot
colourmap), FA (greyscale), and MD (jet colourmap) side by side.

**What to check:**
- NDE and FA maps should show **inverted contrast**: bright FA = dark NDE
  and vice versa, particularly in the corpus callosum and corticospinal tracts
- Check for **ringing artefacts** (concentric bands near bright structures) —
  these indicate insufficient Gibbs removal
- Check for **geometric distortion** in the NDE map relative to the T1:
  the brain outline should align with the T1 brain mask
- **Unusually high NDE in deep white matter** (approaching 0.95) relative to
  cortex suggests CSF partial-voluming of subcortical WM voxels

---

## nde_pipeline_report.json

The JSON report records all QC thresholds with pass/fail status.
Key fields:

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
    "wm_nde_mean": 0.821,
    "wm_nde_pass": true,
    "nde_fa_pearson_r": -0.951,
    "nde_fa_anticorrelated": true
  }
}
```

A run is considered **high quality** when all four boolean QC fields are `true`
and `nonphysical_pct` is below 5%.

---

## Common Issues and Fixes

| Symptom | Likely cause | Fix |
|---|---|---|
| `nonphysical_pct` > 5% | Poor tensor fit; low SNR; wrong bvecs | Check bvec sign convention; increase b=0 volumes; review SNR |
| `csf_nde_pass: false` (CSF NDE < 0.90) | T1 registration failure; FAST misclassification | Re-run with `--use_ants`; inspect FAST segmentation manually |
| `wm_nde_pass: false` (WM NDE ≥ 0.85) | CSF contamination of WM mask; eddy failure | Lower BET threshold; inspect eddy QC files (`*.eddy_outlier_report`) |
| NDE–FA r > −0.80 | Corrupted tensor | Check bvec directions; re-run eddy with `--repol`; verify b-value |
| Gibbs ringing in NDE map | Insufficient `mrdegibbs` | Ensure raw DWI is passed (not already Gibbs-corrected); check partial Fourier settings |
