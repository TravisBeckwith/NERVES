# Reliability and QC tools

NERVES v2 adds two tools that need only Python, NumPy, SciPy and nibabel (no FSL
or MRtrix3):

* `nerves_reliability.py` — how far can you trust FA, MD, tensor mode and NDE at a
  given SNR, direction count and b-value?
* `nerves_qc.py` — fit-health and reliability QC on tensor maps that already exist
  (QSIRecon workflows, or NERVES/FSL `dti_*` outputs), with group outlier flags.

## 1. Planning calculator — `nerves_reliability.py`

```bash
python nerves_reliability.py --sweep-snr 10 15 20 30 40 --tissue WM
python nerves_reliability.py --sweep-snr 10 15 20 30 40 --tissue GM
python nerves_reliability.py --required FA:0.03 --tissue GM --dirs 60
```

It runs a Rician-noise Monte Carlo (uniform random tensor orientations, the same
OLS/WLS fit the pipeline uses) and reports bias and voxelwise SD. Tissue presets,
in 10⁻³ mm²/s: `WM` (1.7/0.35/0.30, FA 0.78), `WM_LOWFA` (1.2/0.55/0.40, FA 0.53),
`GM` (0.95/0.75/0.65, FA 0.19); or pass `--tissue l1,l2,l3`.

Results at 30 directions, b = 1000, OLS (4000 voxels, seed 0):

**White matter (FA 0.78)**

| SNR | FA bias ± sd | MD bias ± sd | Mode bias ± sd | NDE bias ± sd |
|---|---|---|---|---|
| 10 | +0.002 ± 0.067 | −0.002 ± 0.076 | −0.060 ± 0.067 | −0.026 ± 0.097 |
| 20 | +0.002 ± 0.034 | +0.000 ± 0.038 | −0.015 ± 0.019 | −0.008 ± 0.046 |
| 40 | +0.000 ± 0.017 | +0.000 ± 0.019 | −0.004 ± 0.007 | −0.002 ± 0.022 |

**Grey matter (FA 0.19)**

| SNR | FA bias ± sd | MD bias ± sd | Mode bias ± sd | NDE bias ± sd |
|---|---|---|---|---|
| 10 | **+0.087** ± 0.083 | −0.000 ± 0.071 | −0.378 ± **0.563** | −0.017 ± 0.017 |
| 20 | +0.024 ± 0.046 | −0.000 ± 0.035 | −0.208 ± 0.471 | −0.004 ± 0.007 |
| 40 | +0.006 ± 0.024 | −0.000 ± 0.018 | −0.067 ± 0.290 | −0.001 ± 0.003 |

What this says:

* **FA is biased upward in low-FA tissue at low SNR**, and not in high-FA tissue
  (consistent with published in-vivo findings). At SNR 10 the GM FA bias (+0.087) is
  about 45 % of the true FA (0.19).
* **MD is essentially unbiased** at all SNRs, only noisier.
* **Tensor mode is not interpretable in low-FA tissue**: voxelwise SD is 0.3–0.6 on a
  [−1, 1] scale even at SNR 40. In WM it is stable. Do not interpret mode (or any
  shape metric derived from it) in GM or CSF-adjacent voxels.
* NDE is the least noise-biased of the four in low-FA tissue because it saturates near
  1 there, which is also why it carries little information in that regime.
* `--required FA:0.03` finds the lowest grid SNR at which both |bias| and SD are
  ≤ 0.03: 40 for GM vs 25 for WM at 30 directions (also 25 for GM with 60 directions).

Limits: one tissue tensor per call, no partial volume, crossing fibres, motion or
eddy residuals; single shell; independent Rician noise. Use it for relative
comparisons and sanity checks, not as exact error bars. SNR here is the b = 0 SNR
(mean b = 0 signal / noise σ), e.g. from `dwidenoise -noise`.

Library use:

```python
import nerves_reliability as rel
r = rel.expected_error(snr=18, n_dirs=32, b=1000, tissue="GM")
r["FA"]["bias"], r["FA"]["sd"], r["MODE"]["sd"]
rel.min_snr_for_bias("FA", 0.03, tissue="WM")
rel.roi_sd(r["FA"]["sd"], n_voxels=200)   # ROI mean, independent voxels (optimistic)
```

## 2. Tensor-fit QC — `nerves_qc.py`

```bash
python nerves_qc.py /path/to/derivatives --snr 20 --dirs 30 --b 1000
```

### What it reads

* **QSIRecon** tensor scalars: `*_model-tensor_param-{fa,md,ad,rd}_dwimap.nii.gz`
  (DSI Studio, TORTOISE, DIPY workflows). Each (workflow folder, filename prefix) is one
  set. Optional `l1/l2/l3` parameters enable eigenvalue checks and NDE/mode medians.
* **FSL / NERVES**: `<name>_FA`, `<name>_MD`, `<name>_L1/L2/L3` in the same directory
  (dtifit's default name is `dti`; any `-o` basename works). A set needs FA and MD with
  the same basename. A `nde_pipeline_report.json` in the same directory or its parent
  supplies the SNR automatically; `--snr` overrides it. The report's SNR is the median
  b = 0 SNR over the brain mask, which probably overstates tissue SNR slightly
  (CSF is bright at b = 0); this has not been quantified.

`work/`, `sourcedata/`, `log/`, `logs/`, `figures/`, `tmp/` and `.git/` folders under the
root are never searched, so intermediate pipeline copies are not mistaken for results.
The QSIRecon workflow is taken from the last `qsirecon-<workflow>` folder in the path.
**Not yet verified on real QSIRecon output**; the matching follows the documented
`..._model-tensor_param-<fa|md|ad|rd>_dwimap` naming.

A set needs at least FA and MD. The brain mask is inferred from voxels with finite,
non-zero MD. MD units (mm²/s vs 10⁻³ mm²/s) are detected automatically.

### Checks and thresholds

| Check | Warn | Fail |
|---|---|---|
| Non-physical voxels (FA outside [0,1], MD ≤ 0, RD < 0, AD < RD, eigenvalue ≤ 0) | > 1 % | > 5 % |
| FA ≥ 0.99 (noise floor / failed fit) | > 2 % of voxels | |
| MD > 4×10⁻³ mm²/s | > 5 % of voxels | |
| Median MD outside 0.3–2.0 ×10⁻³ mm²/s (wrong units / failed fit) | warn | |
| MD vs (AD+2·RD)/3 median relative error | > 5 % | |
| FA/MD/AD/RD maps on different grids or affines | | error → FAIL |
| Eigenvalues not sorted L1 ≥ L2 ≥ L3 | > 1 % | |
| Expected GM-like FA bias at the given SNR | > 0.05 | |
| Expected GM-like mode SD at the given SNR | > 0.5 | |
| Group outlier (robust z, per workflow, ≥ 5 sets) | \|z\| > 3.5 | |

Sets are compared only within a group: workflow plus filename prefix with
`sub-`/`ses-`/`run-` removed, so ACPC and MNI space, or tensor and DKI models, are never
mixed. Sessions of the same subject count as separate sets (this violates independence
in longitudinal data; interpret flags accordingly).

Group outlier metrics: mean FA, FA 95th percentile, median MD, non-physical fraction,
FA-saturation fraction, fraction FA > 0.4, and median NDE when eigenvalues are
available. z = (x − median) / max(1.4826·MAD, 2 % of |median|, 0.005); the floor
prevents near-identical groups from producing enormous z-scores.

### Outputs (`--out`, default `ROOT/nerves_qc`)

* `nerves_qc_summary.tsv` — one row per set, all metrics and flags.
* `nerves_qc_summary.json` — same, machine-readable.
* `nerves_qc_report.html` — colour-coded table.

Exit codes: 0 all pass or warn, 3 at least one FAIL, 1 nothing found.

### Interpreting results

* A WARN from expected error is about the **acquisition**, not a defect in the maps:
  it says low-FA tissue FA will be upward-biased and mode unreliable in that scan.
  Do not pool such subjects with high-SNR ones without accounting for SNR.
* Outlier flags say "unlike the rest of this group", not "wrong". Look at the
  subject before excluding.
* The tool judges the tensor maps as written. It does not detect misregistration,
  residual motion or eddy artefacts; use eddy QUAD/SQUAD, MRIQC or QSIPrep reports for
  those.
