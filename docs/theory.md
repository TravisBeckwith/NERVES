# NDE Theory

## Background

Fractional Anisotropy (FA) is the most widely used scalar diffusion tensor metric,
but it applies a **quadratic (L2-norm) weighting** across all three eigenvalues.
This means FA treats changes in the dominant eigenvalue (λ₁) and the minor
eigenvalues (λ₂, λ₃) with equal mathematical weight.

In many pathological states — particularly those affecting radial diffusivity
asymmetrically, such as early Wallerian degeneration or crossing-fiber
configurations — the clinically meaningful change occurs in the **ratio of
the minor eigenvalues**, not in the dominant eigenvalue. FA is relatively
insensitive to this redistribution.

Normalized Diffusion Entropy (NDE) addresses this gap by applying a
**logarithmic weighting** via Shannon entropy, which is disproportionately
sensitive to changes in small probability values.

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

### Why NDE is sensitive to minor eigenvalue redistribution

The derivative of the entropy contribution function is:

$$\frac{d}{dp}(-p \ln p) = -\ln(p) - 1$$

| p value | Gradient magnitude |
|---|---|
| 0.01 | +3.61 |
| 0.10 | +1.30 |
| 1/e ≈ 0.37 | 0.00 (minimum) |
| 0.90 | −0.89 |
| 1.00 | −1.00 |

The gradient diverges to +∞ as p → 0. This means **small changes in minor
eigenvalue proportions produce large changes in NDE**, while the dominant
eigenvalue (large p) contributes with a gradient near −1. FA's L2-norm
weighting applies gradient magnitudes that are roughly symmetric across
all eigenvalues.

### Numerical demonstration (Manuscript Table 1)

| Config | λ₁ | λ₂ | λ₃ | FA | NDE | Shape |
|---|---|---|---|---|---|---|
| A | 1.5 | 0.5 | 0.5 | 0.603 | 0.865 | Prolate |
| B | 1.5 | 0.75 | 0.25 | 0.643 | 0.817 | Asymmetric |
| C | 1.2 | 0.6 | 0.6 | 0.408 | 0.946 | Prolate |
| D | 1.2 | 0.8 | 0.4 | 0.490 | 0.921 | Asymmetric |

Eigenvalues in units of ×10⁻³ mm²/s.

Between Configs A and B: |ΔNDE| = 0.048, |ΔFA| = 0.040 — NDE shows a
larger absolute difference on the shared [0,1] scale, demonstrating
complementary sensitivity to minor-eigenvalue redistribution.

---

## Relationship to FA and Tensor Mode

NDE is strongly anticorrelated with FA globally (Pearson r typically −0.92
to −0.98 in white matter), but **within narrow FA strata**, NDE spans a
substantial range driven by the λ₂/λ₃ ratio.

This is the key complementarity: NDE and FA together provide a richer
characterisation of the eigenvalue spectrum than either metric alone.

NDE does not map onto tensor mode (Ennis & Kindlmann, 2006) in a simple
way. Tensor mode distinguishes prolate (mode = +1) from oblate (mode = −1)
configurations, while NDE is sensitive to the magnitude of the eigenvalue
asymmetry regardless of its direction.

---

## Edge Cases

| Condition | NDE value | Handling |
|---|---|---|
| λ₁ = λ₂ = λ₃ (isotropic) | 1.0 | Computed normally |
| λ₁ >> λ₂ = λ₃ (prolate) | Lower, near 0.865 for typical WM | Computed normally |
| Any λᵢ ≤ 0 | 0 (flagged) | Voxel written to QC flag mask |
| Outside brain mask | 0 | Not computed |

---

## References

1. Fozouni N, Chopp M, Nejad-Davarani SP, et al. Characterizing brain
   structures and remodeling after TBI based on information content,
   diffusion entropy. *PLoS One*. 2013;8(10):e76343.

2. Özarslan E, Vemuri BC, Mareci TH. Generalized scalar measures for
   diffusion tensor imaging using trace, variance, and entropy.
   *Magn Reson Med*. 2005;53(4):866–876.

3. Ennis DB, Kindlmann G. Orthogonal tensor invariants and the analysis
   of diffusion tensor magnetic resonance images.
   *Magn Reson Med*. 2006;55(1):136–146.
