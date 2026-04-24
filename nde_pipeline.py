#!/usr/bin/env python3
"""
nde_pipeline.py
===============
Normalized Diffusion Entropy (NDE) Pipeline
--------------------------------------------
End-to-end pipeline for computing NDE maps from raw T1 and single-shell
DWI data.  Wraps FSL and MRtrix3 shell tools via subprocess, then
performs all NDE computation natively in Python/NumPy/nibabel.

Mathematical basis (Normalized Diffusion Entropy)
--------------------------------------------------
For eigenvalues λ₁ ≥ λ₂ ≥ λ₃ > 0 at each voxel, define normalised
proportions:

    pᵢ = λᵢ / (λ₁ + λ₂ + λ₃),   i ∈ {1,2,3}

Shannon entropy of the distribution:

    H = −Σ pᵢ ln(pᵢ)

Normalise to [0, 1] by dividing by the maximum entropy of a
three-element distribution, ln(3):

    NDE = H / ln(3)

NDE = 1  → perfectly isotropic  (λ₁ = λ₂ = λ₃)
NDE = 0  → maximally anisotropic (all diffusion along one eigenvector)

Key sensitivity property: the derivative d/dp(−p ln p) = −ln(p) − 1
diverges to +∞ as p → 0, making NDE disproportionately sensitive to
changes in the minor eigenvalue proportions (λ₂, λ₃) relative to the
dominant eigenvalue (λ₁).  FA applies a quadratic (L2-norm) weighting
that is symmetric across all eigenvalues.  NDE therefore carries
complementary information to FA: two voxels can share identical FA but
differ in NDE whenever the ratio λ₂/λ₃ differs.

Reference: Fozouni N, Chopp M, Nejad-Davarani SP, et al. Characterizing
brain structures and remodeling after TBI based on information content,
diffusion entropy. PLoS One. 2013;8(10):e76343.

Pipeline overview
-----------------
1. Validate inputs and external tools (FSL, MRtrix3, ANTs)
2. DWI preprocessing
   a. MP-PCA denoising          (MRtrix3: dwidenoise)
   b. Gibbs ringing removal      (MRtrix3: mrdegibbs)
   c. Brain masking of b=0       (FSL: bet)
   d. Topup (optional)           (FSL: topup)
   e. Eddy-current / motion      (FSL: eddy)
3. T1 preprocessing
   a. Brain extraction           (FSL: bet)
   b. Tissue segmentation        (FSL: fast)
   c. Registration T1 → DWI     (FSL: flirt / ANTs)
4. Diffusion tensor estimation   (FSL: dtifit)
5. NDE computation               (Python/NumPy)
6. Quality control               (Python/matplotlib)
7. Summary report

Usage
-----
    python nde_pipeline.py \\
        --dwi  sub-01_dwi.nii.gz \\
        --bval sub-01_dwi.bval   \\
        --bvec sub-01_dwi.bvec   \\
        --t1   sub-01_T1w.nii.gz \\
        --out  ./nde_output       \\
        [--pe_dir AP]             \\
        [--rpe_b0 sub-01_b0_PA.nii.gz]  \\
        [--no_topup]              \\
        [--use_ants]              \\
        [--skip_preproc]          \\
        [--l1 dti_L1.nii.gz --l2 dti_L2.nii.gz --l3 dti_L3.nii.gz \\
         --mask brain_mask.nii.gz]

Dependencies
------------
Python  : numpy, nibabel, scipy, matplotlib, seaborn, pandas  (pip)
External: FSL ≥ 6.0, MRtrix3 ≥ 3.0.3, ANTs ≥ 2.3 (optional)
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import nibabel as nib
from scipy.stats import gaussian_kde, pearsonr

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("nde_pipeline")


# ===========================================================================
# 1. Utility helpers
# ===========================================================================

def run(cmd: str | list, cwd: Optional[Path] = None,
        check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command, logging it, raising on failure."""
    if isinstance(cmd, str):
        cmd_str = cmd
    else:
        cmd_str = " ".join(str(c) for c in cmd)
    log.info("CMD: %s", cmd_str)
    result = subprocess.run(
        cmd_str, shell=True, capture_output=True, text=True, cwd=cwd
    )
    if result.stdout.strip():
        log.debug("STDOUT: %s", result.stdout.strip())
    if result.stderr.strip():
        log.debug("STDERR: %s", result.stderr.strip())
    if check and result.returncode != 0:
        log.error("Command failed (rc=%d): %s", result.returncode, cmd_str)
        log.error("STDERR: %s", result.stderr)
        raise RuntimeError(f"Command failed: {cmd_str}")
    return result


def check_tool(name: str) -> bool:
    return shutil.which(name) is not None


def validate_tools(use_ants: bool = False) -> None:
    """Check required external tools are on PATH."""
    required_fsl = ["bet", "eddy", "dtifit", "flirt", "fast",
                    "fslmerge", "fslmaths"]
    required_mrtrix = ["dwidenoise", "mrdegibbs", "mrconvert",
                       "dwiextract", "mrmath"]
    optional = ["topup", "applywarp"]
    ants = ["antsRegistrationSyNQuick.sh", "antsApplyTransforms"]

    missing = []
    for t in required_fsl + required_mrtrix:
        if not check_tool(t):
            missing.append(t)
    if missing:
        raise EnvironmentError(
            f"Missing required tools: {missing}\n"
            "Install FSL (https://fsl.fmrib.ox.ac.uk) and "
            "MRtrix3 (https://www.mrtrix.org)."
        )
    for t in optional:
        if not check_tool(t):
            log.warning("Optional tool not found: %s", t)
    if use_ants:
        for t in ants:
            if not check_tool(t):
                log.warning("ANTs tool not found: %s (needed for --use_ants)", t)


# ===========================================================================
# 2. DWI Preprocessing
# ===========================================================================

def preprocess_dwi(dwi: Path, bval: Path, bvec: Path,
                   out_dir: Path, pe_dir: str = "AP",
                   rpe_b0: Optional[Path] = None) -> dict:
    """
    Full DWI preprocessing:
      denoising → Gibbs removal → topup (optional) → eddy → bet mask

    Returns dict of output paths.
    """
    log.info("=== DWI Preprocessing ===")
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Convert input NIfTI → MIF for MRtrix3 ---
    dwi_mif = out_dir / "dwi_raw.mif"
    run(f"mrconvert {dwi} -fslgrad {bvec} {bval} {dwi_mif} -force")

    # --- 2a. MP-PCA denoising ---
    log.info("Step 2a: MP-PCA denoising")
    dwi_den = out_dir / "dwi_denoised.mif"
    noise_map = out_dir / "noise_map.nii.gz"
    run(f"dwidenoise {dwi_mif} {dwi_den} -noise {noise_map} -force")

    # --- 2b. Gibbs ringing removal ---
    log.info("Step 2b: Gibbs ringing removal")
    dwi_deg = out_dir / "dwi_degibbs.mif"
    run(f"mrdegibbs {dwi_den} {dwi_deg} -force")

    # Convert degibbs output to NIfTI for FSL tools
    dwi_deg_nii = out_dir / "dwi_degibbs.nii.gz"
    bvecs_deg   = out_dir / "dwi_degibbs.bvec"
    bvals_deg   = out_dir / "dwi_degibbs.bval"
    run(f"mrconvert {dwi_deg} {dwi_deg_nii} "
        f"-export_grad_fsl {bvecs_deg} {bvals_deg} -force")

    # --- 2c. Initial brain mask from mean b=0 ---
    log.info("Step 2c: Initial brain masking")
    mean_b0 = out_dir / "mean_b0.nii.gz"
    run(f"dwiextract {dwi_deg} -bzero - | "
        f"mrmath - mean {mean_b0} -axis 3 -force")
    mean_b0_brain = out_dir / "mean_b0_brain"
    run(f"bet {mean_b0} {mean_b0_brain} -m -f 0.3")
    mask_init = out_dir / "mean_b0_brain_mask.nii.gz"

    # --- 2d. Topup (optional) ---
    topup_results = None
    acqparams = out_dir / "acqparams.txt"

    readout_time = 0.05  # seconds — update to match your acquisition
    pe_map = {
        "AP": (0, -1, 0), "PA": (0, 1, 0),
        "LR": (-1, 0, 0), "RL": (1, 0, 0),
    }
    pe_vec = pe_map.get(pe_dir.upper(), (0, -1, 0))
    rpe_vec = tuple(-v for v in pe_vec)

    if rpe_b0 is not None and check_tool("topup"):
        log.info("Step 2d: Topup susceptibility distortion correction")
        rpe_b0_nii = out_dir / "rpe_b0.nii.gz"
        shutil.copy(rpe_b0, rpe_b0_nii)

        b0_pair = out_dir / "b0_pair.nii.gz"
        run(f"fslmerge -t {b0_pair} {mean_b0} {rpe_b0_nii}")

        with open(acqparams, "w") as f:
            f.write(f"{pe_vec[0]} {pe_vec[1]} {pe_vec[2]} {readout_time}\n")
            f.write(f"{rpe_vec[0]} {rpe_vec[1]} {rpe_vec[2]} {readout_time}\n")

        topup_results = out_dir / "topup_results"
        run(f"topup --imain={b0_pair} "
            f"--datain={acqparams} "
            f"--config=b02b0.cnf "
            f"--out={topup_results} "
            f"--fout={out_dir / 'field_map.nii.gz'}")
    else:
        if rpe_b0 is None:
            log.info("Step 2d: No reverse phase-encode b0 supplied; "
                     "skipping topup.")
        # Write single-line acqparams for eddy
        with open(acqparams, "w") as f:
            f.write(f"{pe_vec[0]} {pe_vec[1]} {pe_vec[2]} {readout_time}\n")

    # --- 2e. Eddy current and motion correction ---
    log.info("Step 2e: Eddy current / motion correction")
    n_vols = nib.load(str(dwi_deg_nii)).shape[3]
    index_file = out_dir / "index.txt"
    with open(index_file, "w") as f:
        f.write(" ".join(["1"] * n_vols) + "\n")

    eddy_out = out_dir / "dwi_eddy"
    eddy_cmd = (
        f"eddy --imain={dwi_deg_nii} "
        f"--mask={mask_init} "
        f"--acqp={acqparams} "
        f"--index={index_file} "
        f"--bvecs={bvecs_deg} "
        f"--bvals={bvals_deg} "
        f"--repol "
        f"--out={eddy_out}"
    )
    if topup_results is not None:
        eddy_cmd += f" --topup={topup_results}"
    run(eddy_cmd)

    eddy_nii    = Path(str(eddy_out) + ".nii.gz")
    eddy_bvecs  = Path(str(eddy_out) + ".eddy_rotated_bvecs")

    # --- Final brain mask from eddy-corrected data ---
    log.info("Step 2e (post): Final brain masking on eddy output")
    mean_b0_eddy = out_dir / "mean_b0_eddy.nii.gz"
    eddy_mif     = out_dir / "dwi_eddy.mif"
    run(f"mrconvert {eddy_nii} -fslgrad {eddy_bvecs} {bvals_deg} "
        f"{eddy_mif} -force")
    run(f"dwiextract {eddy_mif} -bzero - | "
        f"mrmath - mean {mean_b0_eddy} -axis 3 -force")
    mean_b0_final = out_dir / "mean_b0_final_brain"
    run(f"bet {mean_b0_eddy} {mean_b0_final} -m -f 0.3")
    mask_final = out_dir / "mean_b0_final_brain_mask.nii.gz"

    return {
        "dwi_eddy":  eddy_nii,
        "bvecs":     eddy_bvecs,
        "bvals":     bvals_deg,
        "mask":      mask_final,
        "mean_b0":   mean_b0_eddy,
    }


# ===========================================================================
# 3. T1 Preprocessing and Registration
# ===========================================================================

def preprocess_t1(t1: Path, dwi_ref: Path, mask_ref: Path,
                  out_dir: Path, use_ants: bool = False) -> dict:
    """
    T1 brain extraction, tissue segmentation, and registration to DWI space.

    Returns dict of output paths including tissue probability maps.
    """
    log.info("=== T1 Preprocessing ===")
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- 3a. Brain extraction ---
    log.info("Step 3a: T1 brain extraction")
    t1_brain = out_dir / "T1_brain"
    run(f"bet {t1} {t1_brain} -R -f 0.35 -g 0")
    t1_brain_nii = Path(str(t1_brain) + ".nii.gz")

    # --- 3b. Tissue segmentation with FAST ---
    log.info("Step 3b: FSL FAST tissue segmentation (3-class)")
    fast_base = out_dir / "T1_fast"
    run(f"fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 "
        f"-o {fast_base} {t1_brain_nii}")
    # FAST outputs: _pve_0 = CSF, _pve_1 = GM, _pve_2 = WM
    pve_csf = Path(str(fast_base) + "_pve_0.nii.gz")
    pve_gm  = Path(str(fast_base) + "_pve_1.nii.gz")
    pve_wm  = Path(str(fast_base) + "_pve_2.nii.gz")

    # --- 3c. Registration T1 → DWI space ---
    log.info("Step 3c: Registering T1 → DWI space")

    if use_ants and check_tool("antsRegistrationSyNQuick.sh"):
        log.info("  Using ANTs SyN registration")
        ants_prefix = out_dir / "T1_to_DWI_"
        run(f"antsRegistrationSyNQuick.sh -d 3 "
            f"-f {dwi_ref} -m {t1_brain_nii} "
            f"-o {ants_prefix} -t a")

        # Apply transform to PVE maps
        warped_pves = {}
        for name, pve in [("csf", pve_csf), ("gm", pve_gm), ("wm", pve_wm)]:
            out_pve = out_dir / f"T1_pve_{name}_DWIspace.nii.gz"
            run(f"antsApplyTransforms -d 3 "
                f"-i {pve} -r {dwi_ref} "
                f"-t {ants_prefix}0GenericAffine.mat "
                f"-n Linear -o {out_pve}")
            warped_pves[name] = out_pve
    else:
        if use_ants:
            log.warning("ANTs not found; falling back to FSL FLIRT")
        log.info("  Using FSL FLIRT affine registration")
        xfm = out_dir / "T1_to_DWI.mat"
        t1_in_dwi = out_dir / "T1_brain_DWIspace.nii.gz"
        run(f"flirt -in {t1_brain_nii} -ref {dwi_ref} "
            f"-out {t1_in_dwi} -omat {xfm} "
            f"-dof 6 -cost normmi")

        warped_pves = {}
        for name, pve in [("csf", pve_csf), ("gm", pve_gm), ("wm", pve_wm)]:
            out_pve = out_dir / f"T1_pve_{name}_DWIspace.nii.gz"
            run(f"flirt -in {pve} -ref {dwi_ref} "
                f"-applyxfm -init {xfm} "
                f"-out {out_pve} -interp trilinear")
            warped_pves[name] = out_pve

    return {
        "t1_brain": t1_brain_nii,
        "pve_csf":  warped_pves["csf"],
        "pve_gm":   warped_pves["gm"],
        "pve_wm":   warped_pves["wm"],
    }


# ===========================================================================
# 4. Diffusion Tensor Estimation
# ===========================================================================

def estimate_tensor(dwi: Path, bvecs: Path, bvals: Path,
                    mask: Path, out_dir: Path) -> dict:
    """
    Fit diffusion tensor and extract eigenvalue maps using FSL DTIFIT.

    Eigenvalues are output in descending order: L1 ≥ L2 ≥ L3
    (convention verified: L1 = axial diffusivity, RD = (L2+L3)/2).
    """
    log.info("=== Diffusion Tensor Estimation (FSL dtifit) ===")
    out_dir.mkdir(parents=True, exist_ok=True)

    dti_prefix = out_dir / "dti"
    run(f"dtifit "
        f"-k {dwi} "
        f"-m {mask} "
        f"-r {bvecs} "
        f"-b {bvals} "
        f"-o {dti_prefix} "
        f"--save_tensor")

    l1 = out_dir / "dti_L1.nii.gz"
    l2 = out_dir / "dti_L2.nii.gz"
    l3 = out_dir / "dti_L3.nii.gz"
    fa = out_dir / "dti_FA.nii.gz"
    md = out_dir / "dti_MD.nii.gz"
    rd = out_dir / "dti_RD.nii.gz"

    # Compute RD = (L2 + L3) / 2 if not already written
    if not rd.exists():
        run(f"fslmaths {l2} -add {l3} -div 2 {rd}")

    # Verify eigenvalue ordering: L1 should exceed L2 everywhere in mask
    l1_data = nib.load(str(l1)).get_fdata()
    l2_data = nib.load(str(l2)).get_fdata()
    mask_data = nib.load(str(mask)).get_fdata().astype(bool)
    violations = np.sum(l1_data[mask_data] < l2_data[mask_data])
    total = np.sum(mask_data)
    if violations > 0:
        log.warning("Eigenvalue ordering check: %d/%d voxels have L1 < L2 "
                    "(%.2f%%). Verify dtifit output.",
                    violations, total, 100 * violations / total)
    else:
        log.info("Eigenvalue ordering check passed: L1 ≥ L2 in all voxels.")

    return {"l1": l1, "l2": l2, "l3": l3, "fa": fa, "md": md, "rd": rd}


# ===========================================================================
# 5. NDE Computation (pure Python / NumPy)
# ===========================================================================

def compute_nde(l1_path: Path, l2_path: Path, l3_path: Path,
                mask_path: Path, out_dir: Path) -> dict:
    """
    Compute Normalized Diffusion Entropy (NDE) from DTI eigenvalue maps.

    NDE = −Σ pᵢ ln(pᵢ) / ln(3),  pᵢ = λᵢ / Σλⱼ

    Sensitivity property: the derivative d/dp(−p ln p) = −ln(p) − 1
    diverges to +∞ as p → 0, so NDE responds disproportionately to
    changes in the minor eigenvalues.  This gives complementary (not
    redundant) information relative to FA's quadratic (L2-norm) weighting.

    Edge-case handling
    ------------------
    - Any λᵢ ≤ 0  → voxel flagged in QC mask, NDE = 0
    - pᵢ = 0      → 0 · ln(0) defined as 0 (L'Hôpital / continuity)
    - All λᵢ equal → NDE = 1 (maximum entropy, perfectly isotropic)
    - Outside mask → NDE = 0 (not computed)
    """
    log.info("=== NDE Computation ===")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load eigenvalue images
    l1_img = nib.load(str(l1_path))
    l1 = l1_img.get_fdata().astype(np.float64)
    l2 = nib.load(str(l2_path)).get_fdata().astype(np.float64)
    l3 = nib.load(str(l3_path)).get_fdata().astype(np.float64)
    mask = nib.load(str(mask_path)).get_fdata().astype(bool)

    # Initialise outputs
    nde_vol   = np.zeros_like(l1, dtype=np.float32)
    qc_flag   = np.zeros_like(l1, dtype=np.uint8)

    # Voxel classification
    valid    = mask & (l1 > 0) & (l2 > 0) & (l3 > 0)
    non_phys = mask & ~valid           # inside mask but non-positive eigs

    qc_flag[non_phys] = 1

    # Stack eigenvalues for vectorised computation: shape (N_valid, 3)
    eigs = np.stack([l1[valid], l2[valid], l3[valid]], axis=-1)

    # Normalise to proportions (sum-normalisation)
    # Scale-invariant: multiplying all eigenvalues by a constant leaves
    # NDE unchanged — only the relative distribution matters.
    eig_sum = eigs.sum(axis=-1, keepdims=True)          # (N, 1)
    p       = eigs / eig_sum                             # (N, 3)

    # Shannon entropy: use np.where so that p=0 → 0·ln(0) = 0 (not NaN)
    log_p = np.where(p > 0.0, np.log(p), 0.0)
    H     = -np.sum(p * log_p, axis=-1)                 # (N,)

    # Normalise: NDE ∈ [0, 1]
    nde_vals = np.clip(H / np.log(3.0), 0.0, 1.0)
    nde_vol[valid] = nde_vals.astype(np.float32)

    # Save NDE map
    nde_path = out_dir / "nde.nii.gz"
    nib.save(
        nib.Nifti1Image(nde_vol, l1_img.affine, l1_img.header),
        str(nde_path)
    )

    # Save QC flag mask
    qc_path = out_dir / "nde_qc_flag.nii.gz"
    nib.save(
        nib.Nifti1Image(qc_flag, l1_img.affine, l1_img.header),
        str(qc_path)
    )

    # Summary statistics (brain mask)
    nde_brain = nde_vol[mask]
    n_mask    = int(mask.sum())
    n_valid   = int(valid.sum())
    n_nonphys = int(non_phys.sum())
    nonphys_pct = 100.0 * n_nonphys / n_mask if n_mask > 0 else 0.0

    stats = {
        "n_mask_voxels":      n_mask,
        "n_valid_voxels":     n_valid,
        "n_nonphysical":      n_nonphys,
        "nonphysical_pct":    round(nonphys_pct, 3),
        "nde_mean":           float(np.mean(nde_brain)),
        "nde_std":            float(np.std(nde_brain)),
        "nde_median":         float(np.median(nde_brain)),
        "nde_min":            float(nde_brain.min()),
        "nde_max":            float(nde_brain.max()),
        "nde_p5":             float(np.percentile(nde_brain, 5)),
        "nde_p95":            float(np.percentile(nde_brain, 95)),
    }

    log.info("NDE summary (brain mask):")
    log.info("  Voxels in mask:        %d", n_mask)
    log.info("  Valid (computed):      %d", n_valid)
    log.info("  Non-physical flagged:  %d  (%.2f%%)", n_nonphys, nonphys_pct)
    if nonphys_pct > 5.0:
        log.warning("Non-physical eigenvalue rate %.2f%% exceeds "
                    "recommended threshold of 5%%.", nonphys_pct)
    log.info("  NDE mean ± std:        %.4f ± %.4f",
             stats["nde_mean"], stats["nde_std"])
    log.info("  NDE range:             [%.4f, %.4f]",
             stats["nde_min"], stats["nde_max"])

    return {"nde": nde_path, "qc_flag": qc_path, "stats": stats}


# ===========================================================================
# 6. Quality Control and Visualisation
# ===========================================================================

def run_qc(nde_path: Path, fa_path: Path, md_path: Path,
           l1_path: Path, l2_path: Path, l3_path: Path,
           mask_path: Path, out_dir: Path,
           pve_csf: Optional[Path] = None,
           pve_gm:  Optional[Path] = None,
           pve_wm:  Optional[Path] = None) -> dict:
    """
    Quality control:
      - NDE / FA / MD histograms within WM mask
      - Tissue-specific NDE distributions (if PVE maps provided)
      - NDE vs FA scatter plot
      - Axial slice montage (NDE, FA, MD side-by-side)
      - Metric correlation matrix
      - Pastes pass/fail on QC thresholds
    """
    log.info("=== Quality Control ===")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    mask  = nib.load(str(mask_path)).get_fdata().astype(bool)
    nde   = nib.load(str(nde_path)).get_fdata()
    fa    = nib.load(str(fa_path)).get_fdata()
    md    = nib.load(str(md_path)).get_fdata()
    l1    = nib.load(str(l1_path)).get_fdata()
    l2    = nib.load(str(l2_path)).get_fdata()
    l3    = nib.load(str(l3_path)).get_fdata()

    # WM mask: use PVE if available, otherwise FA > 0.2 heuristic
    if pve_wm is not None:
        wm_mask = nib.load(str(pve_wm)).get_fdata() > 0.7
    else:
        wm_mask = mask & (fa > 0.2)
    wm_mask = wm_mask & mask

    qc_results = {}

    # ------------------------------------------------------------------
    # Figure 1: Histogram panel (NDE, FA, MD within WM)
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Metric distributions within white matter mask", fontsize=14)

    for ax, data, label, colour, xlim in [
        (axes[0], nde[wm_mask], "NDE", "#E05C5C", (0, 1)),
        (axes[1], fa[wm_mask],  "FA",  "#5C7BE0", (0, 1)),
        (axes[2], md[wm_mask],  "MD (×10⁻³ mm²/s)", "#5CBF5C", None),
    ]:
        data_plot = data * 1000 if label.startswith("MD") else data
        ax.hist(data_plot, bins=80, color=colour, alpha=0.75,
                edgecolor="none", density=True)
        ax.set_xlabel(label, fontsize=12)
        ax.set_ylabel("Density", fontsize=12)
        if xlim is not None:
            ax.set_xlim(*xlim)
        ax.axvline(np.median(data_plot), color="k", linewidth=1.5,
                   linestyle="--", label=f"Median={np.median(data_plot):.3f}")
        ax.legend(fontsize=10)

    fig.tight_layout()
    fig.savefig(out_dir / "qc_histograms.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: qc_histograms.png")

    # ------------------------------------------------------------------
    # Figure 2: Tissue-specific NDE distributions
    # ------------------------------------------------------------------
    if pve_csf is not None and pve_gm is not None and pve_wm is not None:
        csf_data = nib.load(str(pve_csf)).get_fdata()
        gm_data  = nib.load(str(pve_gm)).get_fdata()
        wm_data  = nib.load(str(pve_wm)).get_fdata()

        tissue_masks = {
            "CSF":         mask & (csf_data > 0.9),
            "Gray matter": mask & (gm_data > 0.9),
            "GM-WM boundary": mask & (gm_data > 0.3) & (gm_data < 0.7) &
                               (wm_data > 0.3) & (wm_data < 0.7),
            "White matter": mask & (wm_data > 0.9),
        }
        tissue_colours = ["#00BFFF", "#A0A0A0", "#FFA500", "#FFFFFF"]

        fig, ax = plt.subplots(figsize=(11, 6))
        ax.set_facecolor("#1A1A2E")
        fig.patch.set_facecolor("#1A1A2E")

        x_range = np.linspace(0, 1, 500)
        for (name, tmask), colour in zip(tissue_masks.items(), tissue_colours):
            vals = nde[tmask]
            if len(vals) > 200:
                try:
                    kde = gaussian_kde(vals, bw_method=0.08)
                    ax.fill_between(x_range, kde(x_range), alpha=0.25,
                                    color=colour)
                    ax.plot(x_range, kde(x_range), linewidth=2,
                            color=colour,
                            label=f"{name}  (n={len(vals):,},"
                                  f" μ={vals.mean():.3f})")
                except Exception:
                    pass  # Skip if KDE fails (e.g., all-zero tissue)

        ax.set_xlabel("NDE", fontsize=13, color="white")
        ax.set_ylabel("Density", fontsize=13, color="white")
        ax.set_title("NDE Distribution by Tissue Type", fontsize=14,
                     color="white")
        ax.tick_params(colors="white")
        ax.set_xlim(0, 1)
        leg = ax.legend(fontsize=10, facecolor="#2E2E4E",
                        labelcolor="white", framealpha=0.8)
        fig.tight_layout()
        fig.savefig(out_dir / "qc_tissue_distributions.png", dpi=200,
                    bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close(fig)
        log.info("Saved: qc_tissue_distributions.png")

        # QC threshold checks
        csf_vals = nde[tissue_masks["CSF"]]
        wm_vals  = nde[tissue_masks["White matter"]]
        if len(csf_vals) > 50:
            csf_mean = float(csf_vals.mean())
            qc_results["csf_nde_mean"] = csf_mean
            qc_results["csf_nde_pass"] = csf_mean >= 0.90
            if not qc_results["csf_nde_pass"]:
                log.warning("QC FAIL: Mean NDE in CSF = %.3f (expected ≥ 0.90)",
                            csf_mean)
        if len(wm_vals) > 50:
            wm_mean = float(wm_vals.mean())
            qc_results["wm_nde_mean"] = wm_mean
            qc_results["wm_nde_pass"] = wm_mean < 0.85
            if not qc_results["wm_nde_pass"]:
                log.warning("QC WARN: Mean NDE in WM = %.3f (expected < 0.85)."
                            " Check preprocessing quality.", wm_mean)

    # ------------------------------------------------------------------
    # Figure 3: NDE vs FA scatter (random subsample within WM)
    # ------------------------------------------------------------------
    nde_wm = nde[wm_mask]
    fa_wm  = fa[wm_mask]

    n_sample = min(50_000, len(nde_wm))
    rng = np.random.default_rng(42)
    idx = rng.choice(len(nde_wm), size=n_sample, replace=False)

    r, pval = pearsonr(fa_wm[idx], nde_wm[idx])
    qc_results["nde_fa_pearson_r"] = float(r)
    qc_results["nde_fa_anticorrelated"] = r < -0.80

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(fa_wm[idx], nde_wm[idx], s=1, alpha=0.15,
               c=fa_wm[idx], cmap="plasma", rasterized=True)
    ax.set_xlabel("FA", fontsize=13)
    ax.set_ylabel("NDE", fontsize=13)
    ax.set_title(f"NDE vs FA (white matter, n={n_sample:,})\n"
                 f"Pearson r = {r:.3f}  [expected −0.95 to −0.98]",
                 fontsize=13)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    # Annotate expected complementarity zone (intermediate FA)
    ax.axvspan(0.3, 0.6, alpha=0.07, color="yellow",
               label="Max dissociation zone (FA 0.3–0.6)")
    ax.legend(fontsize=10)
    fig.tight_layout()
    fig.savefig(out_dir / "qc_nde_fa_scatter.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: qc_nde_fa_scatter.png  (r=%.3f)", r)

    # ------------------------------------------------------------------
    # Figure 4: Axial slice montage (NDE | FA | MD)
    # ------------------------------------------------------------------
    _make_montage(nde, fa, md, mask, out_dir / "qc_axial_montage.png")
    log.info("Saved: qc_axial_montage.png")

    return qc_results


def _make_montage(nde: np.ndarray, fa: np.ndarray, md: np.ndarray,
                  mask: np.ndarray, out_path: Path,
                  n_slices: int = 6) -> None:
    """Produce an axial montage: NDE / FA / MD for n equally-spaced slices."""
    # Find slices with brain content
    brain_slices = np.where(mask.any(axis=(0, 1)))[0]
    if len(brain_slices) == 0:
        log.warning("No brain voxels found for montage.")
        return
    step = max(1, len(brain_slices) // (n_slices + 2))
    selected = brain_slices[step:-step:step][:n_slices]

    fig = plt.figure(figsize=(4 * len(selected), 4 * 3))
    gs = gridspec.GridSpec(3, len(selected), hspace=0.05, wspace=0.05)
    row_labels = ["NDE", "FA", "MD"]
    volumes    = [nde, fa, md * 1000]   # MD in ×10⁻³ mm²/s
    cmaps      = ["hot", "gray", "jet"]
    vlims      = [(0, 1), (0, 1), (0, 3)]

    for col, sl in enumerate(selected):
        for row, (vol, cmap, vlim, rlabel) in enumerate(
                zip(volumes, cmaps, vlims, row_labels)):
            ax = fig.add_subplot(gs[row, col])
            sldata = np.rot90(vol[:, :, sl])
            ax.imshow(sldata, cmap=cmap, vmin=vlim[0], vmax=vlim[1],
                      interpolation="nearest", aspect="equal")
            ax.axis("off")
            if col == 0:
                ax.set_ylabel(rlabel, fontsize=11)
                ax.text(-0.05, 0.5, rlabel, transform=ax.transAxes,
                        fontsize=12, va="center", ha="right", color="white",
                        fontweight="bold")
            if row == 0:
                ax.set_title(f"z={sl}", fontsize=10, color="white", pad=2)

    fig.patch.set_facecolor("black")
    fig.savefig(out_path, dpi=150, bbox_inches="tight",
                facecolor="black")
    plt.close(fig)


# ===========================================================================
# 7. Summary report
# ===========================================================================

def write_report(out_dir: Path, nde_stats: dict, qc_results: dict,
                 args_dict: dict) -> None:
    """Write a JSON summary report."""
    report = {
        "pipeline": "nde_pipeline.py",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "inputs": args_dict,
        "nde_statistics": nde_stats,
        "qc_results": qc_results,
        "qc_thresholds": {
            "nonphysical_eigenvalue_pct_max": 5.0,
            "csf_nde_mean_min": 0.90,
            "wm_nde_mean_max":  0.85,
            "nde_fa_r_max":     -0.80,
        },
        "reference": (
            "Fozouni N, Chopp M, Nejad-Davarani SP, Zhang ZG, Lehman NL, "
            "Gu S, et al. Characterizing brain structures and remodeling "
            "after TBI based on information content, diffusion entropy. "
            "PLoS One. 2013;8(10):e76343."
        ),
    }
    rpt_path = out_dir / "nde_pipeline_report.json"
    with open(rpt_path, "w") as f:
        json.dump(report, f, indent=2)
    log.info("Report written: %s", rpt_path)


# ===========================================================================
# 8. Main entry point
# ===========================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="NDE Pipeline: Normalized Diffusion Entropy from T1 + DWI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # --- Required inputs ---
    p.add_argument("--dwi",  required=True, type=Path,
                   help="4D DWI NIfTI (.nii / .nii.gz)")
    p.add_argument("--bval", required=True, type=Path,
                   help="FSL bval file")
    p.add_argument("--bvec", required=True, type=Path,
                   help="FSL bvec file")
    p.add_argument("--t1",   required=True, type=Path,
                   help="3D T1-weighted NIfTI")
    p.add_argument("--out",  required=True, type=Path,
                   help="Output directory (created if absent)")

    # --- Optional inputs ---
    p.add_argument("--pe_dir", default="AP",
                   choices=["AP", "PA", "LR", "RL"],
                   help="Phase-encode direction (default: AP)")
    p.add_argument("--rpe_b0", type=Path, default=None,
                   help="Reverse phase-encode b=0 volume for topup")
    p.add_argument("--no_topup", action="store_true",
                   help="Skip topup even if --rpe_b0 is supplied")
    p.add_argument("--use_ants", action="store_true",
                   help="Use ANTs SyN for T1→DWI registration (default: FLIRT)")

    # --- Skip preprocessing (use pre-computed eigenvalues) ---
    p.add_argument("--skip_preproc", action="store_true",
                   help="Skip preprocessing and tensor estimation; "
                        "requires --l1 --l2 --l3 --mask")
    p.add_argument("--l1",   type=Path, default=None,
                   help="Pre-computed λ₁ NIfTI (use with --skip_preproc)")
    p.add_argument("--l2",   type=Path, default=None,
                   help="Pre-computed λ₂ NIfTI (use with --skip_preproc)")
    p.add_argument("--l3",   type=Path, default=None,
                   help="Pre-computed λ₃ NIfTI (use with --skip_preproc)")
    p.add_argument("--mask", type=Path, default=None,
                   help="Brain mask NIfTI (use with --skip_preproc)")

    # --- Optional pre-computed PVE maps ---
    p.add_argument("--pve_csf", type=Path, default=None,
                   help="CSF partial volume estimate (DWI space)")
    p.add_argument("--pve_gm",  type=Path, default=None,
                   help="GM partial volume estimate (DWI space)")
    p.add_argument("--pve_wm",  type=Path, default=None,
                   help="WM partial volume estimate (DWI space)")

    p.add_argument("--verbose", action="store_true",
                   help="Enable DEBUG logging")

    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    out_dir = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("  NDE Pipeline — Normalized Diffusion Entropy")
    log.info("=" * 60)
    log.info("Output directory: %s", out_dir)

    # ------------------------------------------------------------------
    # Validate external tools
    # ------------------------------------------------------------------
    if not args.skip_preproc:
        validate_tools(use_ants=args.use_ants)

    # ------------------------------------------------------------------
    # Track paths that feed into downstream steps
    # ------------------------------------------------------------------
    l1_path   = args.l1
    l2_path   = args.l2
    l3_path   = args.l3
    mask_path = args.mask
    fa_path   = None
    md_path   = None
    pve_csf   = args.pve_csf
    pve_gm    = args.pve_gm
    pve_wm    = args.pve_wm

    if not args.skip_preproc:
        # --------------------------------------------------------------
        # Step 2: DWI preprocessing
        # --------------------------------------------------------------
        rpe_b0 = None if args.no_topup else args.rpe_b0
        preproc = preprocess_dwi(
            dwi=args.dwi, bval=args.bval, bvec=args.bvec,
            out_dir=out_dir / "preproc",
            pe_dir=args.pe_dir,
            rpe_b0=rpe_b0,
        )
        mask_path = preproc["mask"]

        # --------------------------------------------------------------
        # Step 3: T1 preprocessing and registration
        # --------------------------------------------------------------
        t1_results = preprocess_t1(
            t1=args.t1,
            dwi_ref=preproc["mean_b0"],
            mask_ref=mask_path,
            out_dir=out_dir / "t1",
            use_ants=args.use_ants,
        )
        pve_csf = t1_results["pve_csf"]
        pve_gm  = t1_results["pve_gm"]
        pve_wm  = t1_results["pve_wm"]

        # --------------------------------------------------------------
        # Step 4: Tensor estimation
        # --------------------------------------------------------------
        tensor_dir = out_dir / "tensor"
        tensor = estimate_tensor(
            dwi=preproc["dwi_eddy"],
            bvecs=preproc["bvecs"],
            bvals=preproc["bvals"],
            mask=mask_path,
            out_dir=tensor_dir,
        )
        l1_path = tensor["l1"]
        l2_path = tensor["l2"]
        l3_path = tensor["l3"]
        fa_path = tensor["fa"]
        md_path = tensor["md"]

    else:
        # Validate that all required pre-computed files were supplied
        missing = [n for n, p in
                   [("--l1", l1_path), ("--l2", l2_path),
                    ("--l3", l3_path), ("--mask", mask_path)]
                   if p is None]
        if missing:
            log.error("--skip_preproc requires: %s", ", ".join(missing))
            sys.exit(1)
        # Derive FA / MD paths if they exist alongside the eigenvalues
        # (assume dtifit naming convention)
        candidate_fa = l1_path.parent / l1_path.name.replace("L1", "FA")
        candidate_md = l1_path.parent / l1_path.name.replace("L1", "MD")
        fa_path = candidate_fa if candidate_fa.exists() else None
        md_path = candidate_md if candidate_md.exists() else None
        if fa_path is None:
            log.warning("FA map not found alongside L1; QC scatter plot "
                        "will be skipped.")

    # ------------------------------------------------------------------
    # Step 5: NDE computation
    # ------------------------------------------------------------------
    nde_dir = out_dir / "nde"
    nde_result = compute_nde(
        l1_path=l1_path,
        l2_path=l2_path,
        l3_path=l3_path,
        mask_path=mask_path,
        out_dir=nde_dir,
    )

    # ------------------------------------------------------------------
    # Step 6: Quality control
    # ------------------------------------------------------------------
    qc_results = {}
    if fa_path is not None and md_path is not None:
        qc_results = run_qc(
            nde_path=nde_result["nde"],
            fa_path=fa_path,
            md_path=md_path,
            l1_path=l1_path,
            l2_path=l2_path,
            l3_path=l3_path,
            mask_path=mask_path,
            out_dir=out_dir / "qc",
            pve_csf=pve_csf,
            pve_gm=pve_gm,
            pve_wm=pve_wm,
        )
    else:
        log.warning("FA and/or MD maps unavailable — QC figures skipped.")

    # ------------------------------------------------------------------
    # Step 7: Summary report
    # ------------------------------------------------------------------
    write_report(
        out_dir=out_dir,
        nde_stats=nde_result["stats"],
        qc_results=qc_results,
        args_dict={k: str(v) for k, v in vars(args).items()},
    )

    # ------------------------------------------------------------------
    # Final output summary
    # ------------------------------------------------------------------
    log.info("=" * 60)
    log.info("  Pipeline complete.")
    log.info("  NDE map : %s", nde_result["nde"])
    log.info("  QC flag : %s", nde_result["qc_flag"])
    log.info("  Report  : %s", out_dir / "nde_pipeline_report.json")
    log.info("=" * 60)

    # Non-zero exit code if non-physical rate is too high
    if nde_result["stats"]["nonphysical_pct"] > 5.0:
        log.warning("Exiting with code 2: non-physical eigenvalue rate "
                    "exceeds 5%%. Review preprocessing.")
        sys.exit(2)


if __name__ == "__main__":
    main()
