# NDE Theory

## Background

Fractional Anisotropy (FA) is the most widely used scalar diffusion tensor metric.
It is the normalised standard deviation of the eigenvalues, i.e. a quadratic
(L2-type) summary.

**Normalized Diffusion Entropy (NDE)** as implemented here is the Shannon entropy of
the eigenvalue proportions pᵢ = λᵢ / Σλ, divided by ln 3. It is a purely
*tensor-based* quantity.

> **Relation to Fozouni et al. 2013.** That paper computes a *model-free* diffusion
> entropy from the probability distribution of path lengths estimated with Q-ball /
> Q-DTI sampling, motivated precisely by the limits of the single-tensor Gaussian
> model. It also reports *low* entropy for isotropic CSF, whereas NDE is *highest*
> for isotropic diffusion. NERVES is therefore **not** an implementation of that
> estimator; the paper is cited as related work on entropy-based diffusion contrast.
> Entropy-type scalar measures of the tensor itself have been discussed earlier
> (e.g. Özarslan, Vemuri & Mareci 2005).

---

## Mathematical Derivation

### Step 1: Eigenvalue Proportions

Given eigenvalues λ₁ ≥ λ₂ ≥ λ₃ > 0 from the diffusion tensor at a voxel,
normalise to a discrete probability distribution:

$$p_i = \frac{\lambda_i}{\lambda_1 + \lambda_2 + \lambda_3}, \quad i \in \{1, 2, 3\}$$

This normalisation is **scale-invariant**: multiplying all eigenvalues by
a constant leaves NDE unchanged. Only the relative distribution matters.

### Step 2: Shannon Entropy

$$H = -\sum_{i=1}^{3} p_i \ln(p_i)$$

where the convention 0 · ln(0) = 0 applies by continuity (L'Hôpital's rule).

### Step 3: Normalisation

The maximum entropy of a 3-element distribution is ln(3), achieved when
all proportions are equal (p₁ = p₂ = p₃ = 1/3, i.e., isotropic diffusion).

$$\text{NDE} = \frac{H}{\ln(3)} \in [0, 1]$$

---

## Sensitivity Analysis

### The log-weighting argument (and its limits)

The derivative of the entropy contribution function is:

$$\frac{d}{dp}(-p \ln p) = -\ln(p) - 1$$

| p value | Gradient |
|---|---|
| 0.01 | +3.61 |
| 0.10 | +1.30 |
| 1/e ≈ 0.37 | 0.00 |
| 0.90 | −0.89 |
| 1.00 | −1.00 |

The gradient grows without bound as p → 0, so *for very small proportions* NDE changes
rapidly with p. In diffusion tensors of brain tissue, however, the smallest proportion
p₃ is rarely below ~0.1 (a λ₃ of ≥ 0.3×10⁻³ mm²/s against a trace of ~2.4), where the
gradient is modest, and λ₃ is also the eigenvalue most affected by noise and by
eigenvalue-sorting bias. The divergence at p → 0 should therefore **not** be read as
evidence of practical sensitivity to subtle pathology.

### Numerical illustration (Table 1)

| Config | λ₁ | λ₂ | λ₃ | FA | NDE | Shape |
|---|---|---|---|---|---|---|
| A | 1.5 | 0.5 | 0.5 | 0.603 | 0.865 | Prolate |
| B | 1.5 | 0.75 | 0.25 | 0.643 | 0.817 | Asymmetric |
| C | 1.2 | 0.6 | 0.6 | 0.408 | 0.946 | Prolate |
| D | 1.2 | 0.8 | 0.4 | 0.463 | 0.921 | Asymmetric |

Eigenvalues in units of ×10⁻³ mm²/s.

These pairs differ in **both** FA and NDE, so they do not demonstrate contrast at fixed
FA. Between A and B, |ΔNDE| = 0.048 > |ΔFA| = 0.040; between C and D,
|ΔNDE| = 0.025 < |ΔFA| = 0.055. Absolute differences of two metrics are in any case
not directly comparable, since the metrics have different dynamic ranges.

### What NDE adds to FA (quantified)

NDE depends on the eigenvalues only through their normalised shape, which has two
degrees of freedom. FA and tensor mode (Ennis & Kindlmann 2006) are two orthogonal
shape invariants, so together they determine NDE. For eigenvalues sampled over a
white-matter-like range (`tests/test_complementarity.py`):

| Predictors of NDE | R² |
|---|---|
| FA alone | ≈ 0.98 |
| FA and tensor mode | > 0.999 |

At fixed FA, NDE does vary with the λ₂/λ₃ split, so it carries mode-like information
beyond FA. That information is available directly from tensor mode, which has an
established interpretation (prolate vs oblate). Any claim that NDE is more sensitive
than FA (or mode) to a particular pathology needs empirical validation that this
repository does not provide.

---

## Edge Cases

| Condition | NDE value | Handling |
|---|---|---|
| λ₁ = λ₂ = λ₃ (isotropic) | 1.0 | Computed normally |
| λ₁ >> λ₂ = λ₃ (prolate) | Lower (e.g. 0.865 for λ = 1.5/0.5/0.5) | Computed normally |
| Any λᵢ ≤ 0 | NaN (undefined) | Voxel written to QC flag mask; excluded from statistics and QC |
| Outside brain mask | 0 | Not computed; excluded from statistics |

NaN (rather than 0) is used for non-physical voxels because 0 is a legitimate NDE value
meaning "maximally anisotropic", and voxels with λ₃ ≤ 0 are typically noisy,
high-anisotropy voxels; writing 0 would bias summary statistics and QC.

---

## References

1. Fozouni N, Chopp M, Nejad-Davarani SP, et al. Characterizing brain
   structures and remodeling after TBI based on information content,
   diffusion entropy. *PLoS One*. 2013;8(10):e76343. (Related; model-free
   estimator, not the one implemented here.)

2. Özarslan E, Vemuri BC, Mareci TH. Generalized scalar measures for
   diffusion tensor imaging using trace, variance, and entropy.
   *Magn Reson Med*. 2005;53(4):866–876.

3. Ennis DB, Kindlmann G. Orthogonal tensor invariants and the analysis
   of diffusion tensor magnetic resonance images.
   *Magn Reson Med*. 2006;55(1):136–146.
