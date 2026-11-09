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

Relation to FA: NDE is a non-linear function of the same eigenvalues as
FA.  At fixed FA it still varies with the λ₂/λ₃ split (i.e. with tensor
mode), but in simulation FA alone explains ~98 % of its variance and
FA + mode > 99.9 %.  Any added sensitivity to pathology is a hypothesis,
not something this software demonstrates (see docs/theory.md).

Relation to prior work
----------------------
This is Shannon entropy of the *tensor eigenvalue spectrum* (a model-based,
Gaussian-diffusion quantity). It is related in spirit to, but is not the same
estimator as, the model-free (Q-ball / mean-free-path) diffusion entropy of:

    Fozouni N, Chopp M, Nejad-Davarani SP, et al. Characterizing brain
    structures and remodeling after TBI based on information content,
    diffusion entropy. PLoS One. 2013;8(10):e76343.

Because NDE is a function of the three eigenvalues only, it is (up to
noise) a function of FA and tensor mode; see docs/theory.md.

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
        [--readout_time 0.062]    \\
        [--registration flirt|epi_reg|ants]  \\
        [--wls] [--max_b 1000] [--min_snr 10]  \\
        [--resume] [--n_threads 8] [--nde_only]

    # Already-preprocessed (e.g. HCP-style) DWI: fit the tensor + NDE only
    python nde_pipeline.py --skip_dwi_preproc --dwi dwi.nii.gz \\
        --bval dwi.bval --bvec dwi.bvec --mask mask.nii.gz --out out

    # Pre-computed eigenvalue maps
    python nde_pipeline.py --skip_preproc \\
        --l1 dti_L1.nii.gz --l2 dti_L2.nii.gz --l3 dti_L3.nii.gz \\
        --mask brain_mask.nii.gz --out out

Dependencies
------------
Python  : numpy, nibabel, scipy, matplotlib  (pip)
External: FSL ≥ 6.0, MRtrix3 ≥ 3.0.3, ANTs ≥ 2.3 (optional);
          none are needed with --skip_preproc, or with --skip_dwi_preproc
          --tensor_fit python and a supplied --mask.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import shlex
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

import nerves_dti

__version__ = "1.1.0-dev"

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
    """Run an external command (no shell), logging it, raising on failure.

    ``cmd`` should be a list of arguments; each element is passed to the
    program verbatim, so paths containing spaces or shell metacharacters
    are safe. A plain string is tokenised with ``shlex.split`` for
    convenience (quote paths yourself in that case). Shell pipelines are
    deliberately unsupported -- write intermediate files instead, so a
    failure in any stage is detected.
    """
    argv = shlex.split(cmd) if isinstance(cmd, str) else [str(c) for c in cmd]
    cmd_str = " ".join(shlex.quote(a) for a in argv)
    log.info("CMD: %s", cmd_str)
    result = subprocess.run(argv, capture_output=True, text=True, cwd=cwd)
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


def find_eddy(gpu: bool = False) -> Optional[str]:
    """Name of the first available ``eddy`` executable.

    FSL >= 6.0.4 ships ``eddy_cpu`` / ``eddy_cuda*`` (older releases
    ``eddy_openmp`` or plain ``eddy``). GPU variants are only considered
    when ``gpu`` is true.
    """
    cands = (["eddy_cuda", "eddy_cuda10.2", "eddy_cuda11.0"] if gpu else [])
    cands += ["eddy_cpu", "eddy_openmp", "eddy"]
    for c in cands:
        if check_tool(c):
            return c
    return None


def validate_tools(use_ants: bool = False, registration: Optional[str] = None,
                   eddy_gpu: bool = False, tensor_fit: str = "dtifit",
                   skip_dwi_preproc: bool = False,
                   has_t1: bool = True) -> None:
    """Check required external tools are on PATH.

    The requirement set depends on the run mode: full preprocessing needs
    FSL and MRtrix3; ``--skip_dwi_preproc`` needs only what the chosen
    tensor fit and (optional) T1 processing use.
    """
    method = registration or ("ants" if use_ants else "flirt")
    required = []
    if not skip_dwi_preproc:
        required += ["bet", "fslmerge", "fslmaths", "dwidenoise", "mrdegibbs",
                     "mrconvert", "dwiextract", "mrmath"]
    if tensor_fit == "dtifit":
        required += ["dtifit"]
    if has_t1:
        required += ["bet", "flirt", "fast", "fslmaths"]
        if method == "epi_reg":
            required += ["epi_reg", "convert_xfm"]
    missing = sorted(set(t for t in required if not check_tool(t)))
    if not skip_dwi_preproc and find_eddy(eddy_gpu) is None:
        missing.append("eddy (eddy_cpu / eddy_openmp / eddy)")
    if missing:
        raise EnvironmentError(
            f"Missing required tools: {missing}\n"
            "Install FSL (https://fsl.fmrib.ox.ac.uk) and "
            "MRtrix3 (https://www.mrtrix.org)."
        )
    if not skip_dwi_preproc:
        for t in ["topup", "applytopup"]:
            if not check_tool(t):
                log.warning("Optional tool not found: %s", t)
    if method == "ants":
        for t in ["antsRegistrationSyNQuick.sh", "antsApplyTransforms"]:
            if not check_tool(t):
                log.warning("ANTs tool not found: %s (needed for ANTs "
                            "registration)", t)


def collect_versions() -> dict:
    """Software versions for the report (best effort, never raises)."""
    import platform
    import scipy
    v = {"nerves": __version__, "python": platform.python_version(),
         "numpy": np.__version__, "nibabel": nib.__version__,
         "scipy": scipy.__version__, "matplotlib": matplotlib.__version__}
    try:
        fsldir = os.environ.get("FSLDIR")
        vf = Path(fsldir) / "etc" / "fslversion" if fsldir else None
        if vf is not None and vf.exists():
            v["fsl"] = vf.read_text().strip().split(":")[0]
    except Exception:
        pass
    try:
        if check_tool("mrconvert"):
            r = subprocess.run(["mrconvert", "-version"], capture_output=True,
                               text=True, timeout=20)
            first = (r.stdout or r.stderr).strip().splitlines()
            if first:
                v["mrtrix3"] = first[0].strip()
    except Exception:
        pass
    try:
        r = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                           capture_output=True, text=True, timeout=10,
                           cwd=Path(__file__).resolve().parent)
        if r.returncode == 0 and r.stdout.strip():
            v["git_commit"] = r.stdout.strip()
    except Exception:
        pass
    return v


def _stage_done(resume: bool, label: str, *paths: Path) -> bool:
    """True (and logged) when --resume is on and every output already exists."""
    if resume and all(Path(p).exists() and Path(p).stat().st_size > 0
                      for p in paths):
        log.info("Resume: skipping %s (outputs exist)", label)
        return True
    return False


def _thr(n_threads: int) -> list:
    """MRtrix3 ``-nthreads`` arguments."""
    return ["-nthreads", str(n_threads)] if n_threads else []


# Phase-encode direction -> acqparams vector. Named directions (AP, PA, LR,
# RL) and BIDS axis codes (i, i-, j, j-, k, k-) are both accepted.
PE_VECTORS = {
    "AP": (0, -1, 0), "PA": (0, 1, 0), "LR": (-1, 0, 0), "RL": (1, 0, 0),
    "i": (1, 0, 0), "i-": (-1, 0, 0), "j": (0, 1, 0), "j-": (0, -1, 0),
    "k": (0, 0, 1), "k-": (0, 0, -1),
}
PE_CHOICES = list(PE_VECTORS)


# ===========================================================================
# 2. DWI Preprocessing
# ===========================================================================

def _mean_b0(dwi_mif: Path, out_nii: Path, tmp_mif: Path,
             n_threads: int = 0) -> None:
    """Mean of the b=0 volumes of ``dwi_mif`` -> ``out_nii``.

    Uses an intermediate file instead of a shell pipe so that a failure in
    ``dwiextract`` is not masked by a successful ``mrmath``.
    """
    run(["dwiextract", dwi_mif, "-bzero", tmp_mif, "-force"] + _thr(n_threads))
    run(["mrmath", tmp_mif, "mean", out_nii, "-axis", "3", "-force"]
        + _thr(n_threads))


def preprocess_dwi(dwi: Path, bval: Path, bvec: Path,
                   out_dir: Path, pe_dir: str = "AP",
                   rpe_b0: Optional[Path] = None,
                   readout_time: float = 0.05,
                   resume: bool = False, n_threads: int = 0,
                   eddy_gpu: bool = False) -> dict:
    """
    Full DWI preprocessing:
      denoising → Gibbs removal → topup (optional) → eddy → bet mask

    Parameters
    ----------
    readout_time : float
        Total readout time in seconds, used in acqparams.txt for topup and
        eddy. Scanner- and sequence-dependent — derive it from your DICOM
        headers as EffectiveEchoSpacing x (PhaseEncodingSteps - 1), or from
        the BIDS ``TotalReadoutTime`` field. The default (0.05s) is a
        placeholder, not a measured value for your acquisition.
    resume : bool
        Skip stages whose output files already exist (eddy can take hours).
    n_threads : int
        Threads for MRtrix3 tools and eddy (0 = tool default).
    eddy_gpu : bool
        Prefer a CUDA build of eddy if installed.

    Returns dict of output paths.
    """
    log.info("=== DWI Preprocessing ===")
    out_dir.mkdir(parents=True, exist_ok=True)
    thr = _thr(n_threads)

    dwi_mif = out_dir / "dwi_raw.mif"
    dwi_den = out_dir / "dwi_denoised.mif"
    noise_map = out_dir / "noise_map.nii.gz"
    dwi_deg = out_dir / "dwi_degibbs.mif"
    dwi_deg_nii = out_dir / "dwi_degibbs.nii.gz"
    bvecs_deg = out_dir / "dwi_degibbs.bvec"
    bvals_deg = out_dir / "dwi_degibbs.bval"

    # --- 2a. MP-PCA denoising (also yields the noise sigma map) ---
    log.info("Step 2a: MP-PCA denoising")
    if not _stage_done(resume, "denoising", dwi_den, noise_map):
        run(["mrconvert", dwi, "-fslgrad", bvec, bval, dwi_mif, "-force"]
            + thr)
        run(["dwidenoise", dwi_mif, dwi_den, "-noise", noise_map, "-force"]
            + thr)

    # --- 2b. Gibbs ringing removal ---
    log.info("Step 2b: Gibbs ringing removal")
    if not _stage_done(resume, "Gibbs removal", dwi_deg, dwi_deg_nii,
                       bvecs_deg, bvals_deg):
        # NOTE: mrdegibbs assumes slices are acquired along the third image
        # axis (axial acquisition); see docs/troubleshooting.md for -axes.
        run(["mrdegibbs", dwi_den, dwi_deg, "-force"] + thr)
        run(["mrconvert", dwi_deg, dwi_deg_nii,
             "-export_grad_fsl", bvecs_deg, bvals_deg, "-force"] + thr)

    # --- 2c. Initial brain mask from mean b=0 ---
    log.info("Step 2c: Initial brain masking")
    mean_b0 = out_dir / "mean_b0.nii.gz"
    mean_b0_brain = out_dir / "mean_b0_brain"
    mask_init = out_dir / "mean_b0_brain_mask.nii.gz"
    if not _stage_done(resume, "initial mask", mean_b0, mask_init):
        _mean_b0(dwi_deg, mean_b0, out_dir / "b0_vols_tmp.mif", n_threads)
        run(["bet", mean_b0, mean_b0_brain, "-m", "-f", "0.3"])

    # --- 2d. Topup (optional) ---
    topup_results = None
    acqparams = out_dir / "acqparams.txt"
    pe_key = pe_dir if pe_dir in PE_VECTORS else pe_dir.upper()
    pe_vec = PE_VECTORS.get(pe_key, (0, -1, 0))
    rpe_vec = tuple(-v for v in pe_vec)

    if rpe_b0 is not None and check_tool("topup"):
        log.info("Step 2d: Topup susceptibility distortion correction")
        rpe_b0_nii = out_dir / "rpe_b0.nii.gz"
        b0_pair = out_dir / "b0_pair.nii.gz"
        with open(acqparams, "w") as f:
            f.write(f"{pe_vec[0]} {pe_vec[1]} {pe_vec[2]} {readout_time}\n")
            f.write(f"{rpe_vec[0]} {rpe_vec[1]} {rpe_vec[2]} {readout_time}\n")
        topup_results = out_dir / "topup_results"
        topup_coef = Path(str(topup_results) + "_fieldcoef.nii.gz")
        mask_topup = out_dir / "mean_b0_topup_brain_mask.nii.gz"
        if not _stage_done(resume, "topup", topup_coef):
            # Average the reverse-PE b0(s) so there is exactly one volume,
            # matching the one acqparams row written for it.
            run(["fslmaths", rpe_b0, "-Tmean", rpe_b0_nii])
            run(["fslmerge", "-t", b0_pair, mean_b0, rpe_b0_nii])
            run(["topup", f"--imain={b0_pair}", f"--datain={acqparams}",
                 "--config=b02b0.cnf", f"--out={topup_results}",
                 f"--fout={out_dir / 'field_map.nii.gz'}"])

        # The eddy brain mask must be defined on the *distortion-corrected*
        # b0, otherwise it is misaligned with the data topup hands to eddy.
        if check_tool("applytopup"):
            if not _stage_done(resume, "topup-corrected mask", mask_topup):
                mean_b0_corr = out_dir / "mean_b0_topup.nii.gz"
                run(["applytopup", f"--imain={mean_b0}",
                     f"--datain={acqparams}", "--inindex=1",
                     f"--topup={topup_results}", f"--out={mean_b0_corr}",
                     "--method=jac"])
                run(["bet", mean_b0_corr, out_dir / "mean_b0_topup_brain",
                     "-m", "-f", "0.3"])
            mask_init = mask_topup
        else:
            log.warning("applytopup not found; eddy mask is defined on the "
                        "uncorrected b0.")
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
    eddy_nii = Path(str(eddy_out) + ".nii.gz")
    eddy_bvecs = Path(str(eddy_out) + ".eddy_rotated_bvecs")
    if not _stage_done(resume, "eddy", eddy_nii, eddy_bvecs):
        eddy_cmd = [
            find_eddy(eddy_gpu) or "eddy",
            f"--imain={dwi_deg_nii}", f"--mask={mask_init}",
            f"--acqp={acqparams}", f"--index={index_file}",
            f"--bvecs={bvecs_deg}", f"--bvals={bvals_deg}",
            "--repol", f"--out={eddy_out}",
        ]
        if n_threads:
            eddy_cmd.append(f"--nthr={n_threads}")
        if topup_results is not None:
            eddy_cmd.append(f"--topup={topup_results}")
        run(eddy_cmd)

    # --- Final brain mask from eddy-corrected data ---
    log.info("Step 2e (post): Final brain masking on eddy output")
    mean_b0_eddy = out_dir / "mean_b0_eddy.nii.gz"
    eddy_mif = out_dir / "dwi_eddy.mif"
    mean_b0_final = out_dir / "mean_b0_final_brain"
    mask_final = out_dir / "mean_b0_final_brain_mask.nii.gz"
    if not _stage_done(resume, "final mask", mean_b0_eddy, mask_final):
        run(["mrconvert", eddy_nii, "-fslgrad", eddy_bvecs, bvals_deg,
             eddy_mif, "-force"] + thr)
        _mean_b0(eddy_mif, mean_b0_eddy, out_dir / "b0_vols_tmp.mif",
                 n_threads)
        run(["bet", mean_b0_eddy, mean_b0_final, "-m", "-f", "0.3"])

    return {
        "dwi_eddy":  eddy_nii,
        "bvecs":     eddy_bvecs,
        "bvals":     bvals_deg,
        "mask":      mask_final,
        "mean_b0":   mean_b0_eddy,
        "mean_b0_brain": Path(str(mean_b0_final) + ".nii.gz"),
        "noise_map": noise_map,
    }


# ===========================================================================
# 3. T1 Preprocessing and Registration
# ===========================================================================

def preprocess_t1(t1: Path, dwi_ref: Path,
                  out_dir: Path, use_ants: bool = False,
                  registration: Optional[str] = None,
                  dwi_ref_brain: Optional[Path] = None,
                  resume: bool = False) -> dict:
    """
    T1 brain extraction, tissue segmentation, and registration to DWI space.

    ``registration``: ``"flirt"`` (6-DOF, normalised mutual information;
    default), ``"epi_reg"`` (FSL boundary-based registration using the WM
    segmentation; needs ``dwi_ref_brain``, a brain-extracted b=0) or
    ``"ants"`` (rigid + affine). ``use_ants`` is kept for compatibility.

    Returns dict of output paths including tissue probability maps.
    """
    log.info("=== T1 Preprocessing ===")
    out_dir.mkdir(parents=True, exist_ok=True)
    method = registration or ("ants" if use_ants else "flirt")
    if method not in ("flirt", "epi_reg", "ants"):
        raise ValueError(f"Unknown registration method: {method}")

    # --- 3a. Brain extraction ---
    log.info("Step 3a: T1 brain extraction")
    t1_brain = out_dir / "T1_brain"
    t1_brain_nii = Path(str(t1_brain) + ".nii.gz")
    if not _stage_done(resume, "T1 brain extraction", t1_brain_nii):
        run(["bet", t1, t1_brain, "-R", "-f", "0.35", "-g", "0"])

    # --- 3b. Tissue segmentation with FAST ---
    log.info("Step 3b: FSL FAST tissue segmentation (3-class)")
    fast_base = out_dir / "T1_fast"
    # FAST outputs: _pve_0 = CSF, _pve_1 = GM, _pve_2 = WM
    pve_csf = Path(str(fast_base) + "_pve_0.nii.gz")
    pve_gm  = Path(str(fast_base) + "_pve_1.nii.gz")
    pve_wm  = Path(str(fast_base) + "_pve_2.nii.gz")
    if not _stage_done(resume, "FAST", pve_csf, pve_gm, pve_wm):
        run(["fast", "-t", "1", "-n", "3", "-H", "0.1", "-I", "4", "-l", "20.0",
             "-o", fast_base, t1_brain_nii])

    # --- 3c. Registration T1 → DWI space ---
    log.info("Step 3c: Registering T1 → DWI space (%s)", method)
    warped_pves = {
        n: out_dir / f"T1_pve_{n}_DWIspace.nii.gz" for n in ("csf", "gm", "wm")
    }
    if _stage_done(resume, "registration", *warped_pves.values()):
        pass
    elif method == "ants" and check_tool("antsRegistrationSyNQuick.sh"):
        # '-t a' = rigid + affine only (no nonlinear SyN stage). Only the
        # affine matrix is applied below, which is what is actually wanted
        # for a within-subject T1 -> DWI alignment.
        log.info("  Using ANTs rigid + affine registration")
        ants_prefix = out_dir / "T1_to_DWI_"
        run(["antsRegistrationSyNQuick.sh", "-d", "3",
             "-f", dwi_ref, "-m", t1_brain_nii,
             "-o", ants_prefix, "-t", "a"])
        for name, pve in [("csf", pve_csf), ("gm", pve_gm), ("wm", pve_wm)]:
            run(["antsApplyTransforms", "-d", "3",
                 "-i", pve, "-r", dwi_ref,
                 "-t", f"{ants_prefix}0GenericAffine.mat",
                 "-n", "Linear", "-o", warped_pves[name]])
    else:
        if method == "ants":
            log.warning("ANTs not found; falling back to FSL FLIRT")
            method = "flirt"
        xfm = out_dir / "T1_to_DWI.mat"
        if method == "epi_reg" and dwi_ref_brain is not None:
            log.info("  Using FSL epi_reg (boundary-based registration)")
            wm_seg = out_dir / "T1_wmseg.nii.gz"
            run(["fslmaths", pve_wm, "-thr", "0.5", "-bin", wm_seg])
            epi_prefix = out_dir / "epi2t1"
            run(["epi_reg", f"--epi={dwi_ref_brain}", f"--t1={t1}",
                 f"--t1brain={t1_brain_nii}", f"--wmseg={wm_seg}",
                 f"--out={epi_prefix}"])
            run(["convert_xfm", "-omat", xfm, "-inverse",
                 Path(str(epi_prefix) + ".mat")])
        else:
            if method == "epi_reg":
                log.warning("epi_reg needs a brain-extracted b=0; falling "
                            "back to FLIRT")
            log.info("  Using FSL FLIRT affine registration")
            t1_in_dwi = out_dir / "T1_brain_DWIspace.nii.gz"
            run(["flirt", "-in", t1_brain_nii, "-ref", dwi_ref,
                 "-out", t1_in_dwi, "-omat", xfm,
                 "-dof", "6", "-cost", "normmi"])
        for name, pve in [("csf", pve_csf), ("gm", pve_gm), ("wm", pve_wm)]:
            run(["flirt", "-in", pve, "-ref", dwi_ref,
                 "-applyxfm", "-init", xfm,
                 "-out", warped_pves[name], "-interp", "trilinear"])

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
                    mask: Path, out_dir: Path,
                    wls: bool = False, method: str = "dtifit",
                    resume: bool = False) -> dict:
    """
    Fit the diffusion tensor and extract eigenvalue maps.

    ``method="dtifit"`` uses FSL DTIFIT (``--wls`` adds weighted least
    squares); ``method="python"`` uses the native NumPy fit in
    ``nerves_dti`` (OLS, or WLS when ``wls``) and needs no FSL.

    Eigenvalues are output in descending order: L1 ≥ L2 ≥ L3.
    """
    log.info("=== Diffusion Tensor Estimation (%s%s) ===",
             "FSL dtifit" if method == "dtifit" else "native NumPy",
             ", WLS" if wls else "")
    out_dir.mkdir(parents=True, exist_ok=True)

    l1 = out_dir / "dti_L1.nii.gz"
    l2 = out_dir / "dti_L2.nii.gz"
    l3 = out_dir / "dti_L3.nii.gz"
    fa = out_dir / "dti_FA.nii.gz"
    md = out_dir / "dti_MD.nii.gz"
    rd = out_dir / "dti_RD.nii.gz"

    if _stage_done(resume, "tensor fit", l1, l2, l3, fa, md):
        pass
    elif method == "python":
        nerves_dti.estimate_tensor_numpy(
            dwi, bvecs, bvals, mask, out_dir,
            method="wls" if wls else "ols")
    else:
        cmd = ["dtifit", "-k", dwi, "-m", mask, "-r", bvecs, "-b", bvals,
               "-o", out_dir / "dti", "--save_tensor"]
        if wls:
            cmd.append("-w")
        run(cmd)

    # Compute RD = (L2 + L3) / 2 if not already written
    if not rd.exists():
        run(["fslmaths", l2, "-add", l3, "-div", "2", rd])

    # Verify eigenvalue ordering (L1 >= L2 >= L3) inside the mask
    l1_data = nib.load(str(l1)).get_fdata()
    l2_data = nib.load(str(l2)).get_fdata()
    l3_data = nib.load(str(l3)).get_fdata()
    mask_data = nib.load(str(mask)).get_fdata().astype(bool)
    v = nerves_dti.ordering_violations(l1_data, l2_data, l3_data, mask_data)
    if v["l1_lt_l2"] or v["l2_lt_l3"]:
        tot = max(v["n_mask"], 1)
        log.warning("Eigenvalue ordering check: L1<L2 in %d and L2<L3 in %d "
                    "of %d voxels (%.2f%% / %.2f%%). Verify the tensor fit.",
                    v["l1_lt_l2"], v["l2_lt_l3"], v["n_mask"],
                    100 * v["l1_lt_l2"] / tot, 100 * v["l2_lt_l3"] / tot)
    else:
        log.info("Eigenvalue ordering check passed: L1 ≥ L2 ≥ L3 in all "
                 "voxels.")

    # Cross-check against dtifit's own tensor file (trace consistency)
    dev = nerves_dti.trace_consistency(
        out_dir / "dti_tensor.nii.gz", l1_data, l2_data, l3_data, mask_data)
    if dev is not None:
        if dev > 1e-2:
            log.warning("Tensor trace and L1+L2+L3 disagree (median relative "
                        "difference %.3g); eigenvalue and tensor files may "
                        "not belong together.", dev)
        else:
            log.info("Tensor trace consistency OK (median rel. diff %.2g).",
                     dev)

    return {"l1": l1, "l2": l2, "l3": l3, "fa": fa, "md": md, "rd": rd}


# ===========================================================================
# 5. NDE Computation (pure Python / NumPy)
# ===========================================================================

def compute_nde(l1_path: Path, l2_path: Path, l3_path: Path,
                mask_path: Path, out_dir: Path,
                snr_path: Optional[Path] = None,
                min_snr: float = 0.0) -> dict:
    """
    Compute Normalized Diffusion Entropy (NDE) from DTI eigenvalue maps.

    NDE = −Σ pᵢ ln(pᵢ) / ln(3),  pᵢ = λᵢ / Σλⱼ

    NDE depends only on the normalised eigenvalue shape, so it is a
    function of FA and tensor mode (see docs/theory.md).

    Edge-case handling
    ------------------
    - Any λᵢ ≤ 0  → voxel flagged in the QC mask and NDE = NaN (undefined).
                    NaN -- not 0 -- because 0 is a valid NDE meaning
                    "maximally anisotropic", and such voxels are typically
                    the noisy, high-anisotropy ones. They are excluded from
                    all summary statistics and QC.
    - SNR < ``min_snr`` (when ``snr_path`` and ``min_snr`` > 0 are given)
                  → flagged with value 2 in the QC mask and NDE = NaN.
                    NDE is biased low at low SNR (see docs/theory.md), so
                    unreliable voxels can be excluded.
    - pᵢ = 0      → 0 · ln(0) defined as 0 (L'Hôpital / continuity)
    - All λᵢ equal → NDE = 1 (maximum entropy, perfectly isotropic)
    - Outside mask → NDE = 0 (not computed; excluded from statistics)
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
    nde_vol[non_phys] = np.nan         # undefined, not "0 = fully anisotropic"

    n_low_snr = 0
    if snr_path is not None and min_snr > 0:
        snr = nib.load(str(snr_path)).get_fdata()
        if snr.shape != valid.shape:
            raise ValueError(f"SNR map shape {snr.shape} != eigenvalue grid "
                             f"{valid.shape}")
        low_snr = valid & ~(snr >= min_snr)
        n_low_snr = int(low_snr.sum())
        qc_flag[low_snr] = 2
        nde_vol[low_snr] = np.nan
        valid = valid & ~low_snr
    elif min_snr > 0:
        log.warning("--min_snr given but no SNR map is available; no SNR "
                    "masking applied.")

    # Stack eigenvalues for vectorised computation: shape (N_valid, 3)
    eigs = np.stack([l1[valid], l2[valid], l3[valid]], axis=-1)

    # Scale-invariant: multiplying all eigenvalues by a constant leaves
    # NDE unchanged -- only the relative distribution matters.
    nde_vals = nerves_dti.nde_from_eigenvalues(eigs)
    nde_vol[valid] = nde_vals.astype(np.float32)

    # Save NDE map
    nde_path = out_dir / "nde.nii.gz"
    nde_hdr = l1_img.header.copy()
    nde_hdr.set_data_dtype(np.float32)
    nib.save(nib.Nifti1Image(nde_vol, l1_img.affine, nde_hdr), str(nde_path))

    # Save QC flag mask
    qc_path = out_dir / "nde_qc_flag.nii.gz"
    flag_hdr = l1_img.header.copy()
    flag_hdr.set_data_dtype(np.uint8)
    nib.save(nib.Nifti1Image(qc_flag, l1_img.affine, flag_hdr), str(qc_path))

    # Summary statistics over VALID voxels only (flagged voxels are NaN)
    nde_brain = nde_vol[valid]
    n_mask    = int(mask.sum())
    n_valid   = int(valid.sum())
    n_nonphys = int(non_phys.sum())
    nonphys_pct = 100.0 * n_nonphys / n_mask if n_mask > 0 else 0.0

    if nde_brain.size == 0:
        raise RuntimeError("No valid (all eigenvalues > 0) voxels inside the "
                           "mask; cannot compute NDE statistics.")

    stats = {
        "n_mask_voxels":      n_mask,
        "n_valid_voxels":     n_valid,
        "n_nonphysical":      n_nonphys,
        "nonphysical_pct":    round(nonphys_pct, 3),
        "n_low_snr":          n_low_snr,
        "nde_mean":           float(np.mean(nde_brain)),
        "nde_std":            float(np.std(nde_brain)),
        "nde_median":         float(np.median(nde_brain)),
        "nde_min":            float(nde_brain.min()),
        "nde_max":            float(nde_brain.max()),
        "nde_p5":             float(np.percentile(nde_brain, 5)),
        "nde_p95":            float(np.percentile(nde_brain, 95)),
    }

    log.info("NDE summary (valid voxels in brain mask):")
    log.info("  Voxels in mask:        %d", n_mask)
    log.info("  Valid (computed):      %d", n_valid)
    log.info("  Non-physical flagged:  %d  (%.2f%%)", n_nonphys, nonphys_pct)
    if n_low_snr:
        log.info("  Low-SNR excluded:      %d  (SNR < %g)", n_low_snr, min_snr)
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
      - Sanity checks: CSF NDE high; tissue ordering CSF > GM > WM.
        The NDE-FA correlation is reported for information only: NDE and
        FA are both deterministic functions of the same eigenvalues, so a
        strong anticorrelation is expected even for badly corrupted
        tensors and cannot serve as a pass/fail quality criterion.
    """
    log.info("=== Quality Control ===")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    mask  = nib.load(str(mask_path)).get_fdata().astype(bool)
    nde   = nib.load(str(nde_path)).get_fdata()
    fa    = nib.load(str(fa_path)).get_fdata()
    md    = nib.load(str(md_path)).get_fdata()

    # Only voxels with a defined NDE take part in QC (flagged voxels are NaN)
    valid_nde = np.isfinite(nde)
    mask = mask & valid_nde

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
        ax.legend(fontsize=10, facecolor="#2E2E4E",
                        labelcolor="white", framealpha=0.8)
        fig.tight_layout()
        fig.savefig(out_dir / "qc_tissue_distributions.png", dpi=200,
                    bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close(fig)
        log.info("Saved: qc_tissue_distributions.png")

        # Sanity checks. Absolute WM NDE depends on acquisition (b-value,
        # SNR, resolution) and on how WM is defined, so it is reported but
        # not thresholded; the *relative* tissue ordering is checked instead.
        csf_vals = nde[tissue_masks["CSF"]]
        gm_vals  = nde[tissue_masks["Gray matter"]]
        wm_vals  = nde[tissue_masks["White matter"]]
        if len(csf_vals) > 50:
            csf_mean = float(csf_vals.mean())
            qc_results["csf_nde_mean"] = csf_mean
            qc_results["csf_nde_pass"] = bool(csf_mean >= 0.90)
            if not qc_results["csf_nde_pass"]:
                log.warning("QC FAIL: Mean NDE in CSF = %.3f (expected ≥ 0.90)",
                            csf_mean)
        if len(gm_vals) > 50:
            qc_results["gm_nde_mean"] = float(gm_vals.mean())
        if len(wm_vals) > 50:
            qc_results["wm_nde_mean"] = float(wm_vals.mean())
        if all(k in qc_results for k in ("csf_nde_mean", "gm_nde_mean",
                                          "wm_nde_mean")):
            order_ok = bool(qc_results["csf_nde_mean"]
                            > qc_results["gm_nde_mean"]
                            > qc_results["wm_nde_mean"])
            qc_results["tissue_ordering_pass"] = order_ok
            if not order_ok:
                log.warning("QC WARN: expected mean NDE CSF > GM > WM, got "
                            "%.3f / %.3f / %.3f. Check registration, "
                            "segmentation and tensor fit.",
                            qc_results["csf_nde_mean"],
                            qc_results["gm_nde_mean"],
                            qc_results["wm_nde_mean"])

    # ------------------------------------------------------------------
    # Figure 3: NDE vs FA scatter (random subsample within WM)
    # ------------------------------------------------------------------
    nde_wm = nde[wm_mask]
    fa_wm  = fa[wm_mask]

    n_sample = min(50_000, len(nde_wm))
    rng = np.random.default_rng(42)
    idx = rng.choice(len(nde_wm), size=n_sample, replace=False)

    r, pval = pearsonr(fa_wm[idx], nde_wm[idx])
    qc_results["nde_fa_pearson_r"] = float(r)   # informational only
    qc_results["nde_variance_explained_by_fa"] = float(r * r)

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(fa_wm[idx], nde_wm[idx], s=1, alpha=0.15,
               c=fa_wm[idx], cmap="plasma", rasterized=True)
    ax.set_xlabel("FA", fontsize=13)
    ax.set_ylabel("NDE", fontsize=13)
    ax.set_title(f"NDE vs FA (white matter, n={n_sample:,})\n"
                 f"Pearson r = {r:.3f}  (informational)",
                 fontsize=13)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    # Annotate expected complementarity zone (intermediate FA)
    ax.axvspan(0.3, 0.6, alpha=0.07, color="yellow",
               label="Intermediate FA (0.3–0.6)")
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
                 args_dict: dict, versions: Optional[dict] = None,
                 extra: Optional[dict] = None) -> None:
    """Write a JSON summary report."""

    class _NumpySafeEncoder(json.JSONEncoder):
        """Fallback encoder for numpy scalar types.

        nde_stats and qc_results are built from numpy computations
        throughout the pipeline (np.mean, np.std, scipy.stats.pearsonr,
        comparisons against numpy floats, ...). It's easy for a numpy
        scalar (np.float32/64, np.int32/64, np.bool_) to end up in one of
        these dicts undetected, since numpy scalars behave like Python
        scalars in almost every context except json.dump(), which raises
        TypeError on them. This encoder converts them instead of crashing
        report writing — see CHANGELOG for the nde_fa_anticorrelated bug
        this exact failure mode caused.
        """

        def default(self, obj):
            if isinstance(obj, np.bool_):
                return bool(obj)
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)

    report = {
        "pipeline": "nde_pipeline.py",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "inputs": args_dict,
        "nde_statistics": nde_stats,
        "qc_results": qc_results,
        "qc_thresholds": {
            "nonphysical_eigenvalue_pct_max": 5.0,
            "csf_nde_mean_min": 0.90,
            "tissue_ordering":  "CSF > GM > WM (mean NDE)",
        },
        "notes": (
            "NDE here is the normalised Shannon entropy of the diffusion "
            "tensor eigenvalue spectrum. It is NOT the model-free (Q-ball) "
            "diffusion entropy of Fozouni et al. 2013, which is cited for "
            "the general idea of entropy-based diffusion contrast. NDE-FA "
            "correlation is informational, not a quality criterion."
        ),
        "reference": (
            "Fozouni N, Chopp M, Nejad-Davarani SP, Zhang ZG, Lehman NL, "
            "Gu S, et al. Characterizing brain structures and remodeling "
            "after TBI based on information content, diffusion entropy. "
            "PLoS One. 2013;8(10):e76343 (related; different estimator)."
        ),
    }
    if versions:
        report["versions"] = versions
    if extra:
        report.update(extra)
    rpt_path = out_dir / "nde_pipeline_report.json"
    with open(rpt_path, "w") as f:
        json.dump(report, f, indent=2, cls=_NumpySafeEncoder)
    log.info("Report written: %s", rpt_path)


# ===========================================================================
# 8. Main entry point
# ===========================================================================

def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="NDE Pipeline: Normalized Diffusion Entropy from T1 + DWI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--version", action="version",
                   version=f"nerves {__version__}")

    # --- Inputs (requirements depend on the run mode; checked in main) ---
    p.add_argument("--dwi",  type=Path, default=None,
                   help="4D DWI NIfTI (.nii / .nii.gz). With "
                        "--skip_dwi_preproc: already corrected data whose "
                        "bvecs/bvals are final")
    p.add_argument("--bval", type=Path, default=None, help="FSL bval file")
    p.add_argument("--bvec", type=Path, default=None, help="FSL bvec file")
    p.add_argument("--t1",   type=Path, default=None,
                   help="3D T1-weighted NIfTI (required for full "
                        "preprocessing; optional with --skip_dwi_preproc, "
                        "enabling tissue-specific QC)")
    p.add_argument("--out",  required=True, type=Path,
                   help="Output directory (created if absent)")

    # --- Acquisition / preprocessing options ---
    p.add_argument("--pe_dir", default="AP", choices=PE_CHOICES,
                   help="Phase-encode direction: AP, PA, LR, RL or the BIDS "
                        "axis code (i, i-, j, j-, k, k-). Default: AP")
    p.add_argument("--rpe_b0", type=Path, default=None,
                   help="Reverse phase-encode b=0 volume(s) for topup "
                        "(4-D inputs are averaged)")
    p.add_argument("--no_topup", action="store_true",
                   help="Skip topup even if --rpe_b0 is supplied")
    p.add_argument("--readout_time", type=float, default=0.05,
                   help="Total readout time in seconds for topup/eddy "
                        "acqparams.txt (BIDS: TotalReadoutTime). Scanner- "
                        "and sequence-dependent — derive from DICOM headers "
                        "as EffectiveEchoSpacing x (PhaseEncodingSteps - 1). "
                        "Default 0.05s is a placeholder, not a measured "
                        "value; set this explicitly for real data. "
                        "(default: %(default)s)")
    p.add_argument("--registration", choices=["flirt", "epi_reg", "ants"],
                   default=None,
                   help="T1→DWI registration: flirt (6-DOF, default), "
                        "epi_reg (FSL boundary-based) or ants (rigid+affine)")
    p.add_argument("--use_ants", action="store_true",
                   help="Alias for --registration ants")
    p.add_argument("--resume", action="store_true",
                   help="Skip preprocessing stages whose outputs already "
                        "exist in --out")
    p.add_argument("--n_threads", type=int, default=0,
                   help="Threads for MRtrix3 tools and eddy (0 = default)")
    p.add_argument("--eddy_gpu", action="store_true",
                   help="Use a CUDA build of eddy if available")

    # --- Run modes ---
    p.add_argument("--skip_dwi_preproc", action="store_true",
                   help="Input DWI is already preprocessed (denoised, "
                        "distortion/eddy corrected; final bvecs): skip "
                        "straight to tensor fit. Needs --dwi --bval --bvec "
                        "and --mask (or FSL bet)")
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
                   help="Brain mask NIfTI (--skip_preproc / "
                        "--skip_dwi_preproc)")
    p.add_argument("--fa", type=Path, default=None,
                   help="Pre-computed FA NIfTI, for QC plots (use with "
                        "--skip_preproc). Overrides dtifit-naming "
                        "auto-detection from --l1; use this if your files "
                        "don't follow the ..._L1/..._FA convention.")
    p.add_argument("--md", type=Path, default=None,
                   help="Pre-computed MD NIfTI, for QC plots (use with "
                        "--skip_preproc). Same override behavior as --fa.")

    # --- Tensor estimation ---
    p.add_argument("--tensor_fit", choices=["dtifit", "python"],
                   default="dtifit",
                   help="Tensor fitting backend: FSL dtifit (default) or the "
                        "native NumPy fit (no FSL needed)")
    p.add_argument("--wls", action="store_true",
                   help="Weighted least squares instead of ordinary least "
                        "squares; fewer negative eigenvalues in noisy data")
    p.add_argument("--max_b", type=float, default=None,
                   help="Fit only volumes with b <= this value (b=0 kept). "
                        "Use for multi-shell data or to drop high-b shells "
                        "(DTI assumes b <= ~1000-1500)")

    # --- Quality control / outputs ---
    p.add_argument("--snr", type=Path, default=None,
                   help="Optional SNR map (NIfTI) for --min_snr when not "
                        "running full preprocessing (which derives it from "
                        "the dwidenoise noise map)")
    p.add_argument("--min_snr", type=float, default=0.0,
                   help="Set NDE to NaN (QC flag 2) where SNR < this value "
                        "(0 = off). NDE is biased low at low SNR")
    p.add_argument("--nde_only", action="store_true",
                   help="Write only the NDE map (skip tensor mode and "
                        "Westin shape maps)")
    p.add_argument("--pve_csf", type=Path, default=None,
                   help="CSF partial volume estimate (DWI space)")
    p.add_argument("--pve_gm",  type=Path, default=None,
                   help="GM partial volume estimate (DWI space)")
    p.add_argument("--pve_wm",  type=Path, default=None,
                   help="WM partial volume estimate (DWI space)")

    p.add_argument("--verbose", action="store_true",
                   help="Enable DEBUG logging")

    return p.parse_args(argv)


def _validate_shapes(named_paths: dict) -> None:
    """Exit with an error if supplied NIfTI inputs differ in grid shape."""
    shapes = {}
    for name, path in named_paths.items():
        if path is None:
            continue
        if not Path(path).exists():
            log.error("%s: file not found: %s", name, path)
            sys.exit(1)
        shapes[name] = tuple(nib.load(str(path)).shape)
    if not shapes:
        return
    ref_name = next(iter(shapes))
    for name, shp in shapes.items():
        if shp != shapes[ref_name]:
            log.error("Shape mismatch: %s is %s but %s is %s. All inputs "
                      "must be on the same voxel grid.",
                      name, shp, ref_name, shapes[ref_name])
            sys.exit(1)


def _fail(msg: str, *a) -> None:
    log.error(msg, *a)
    sys.exit(1)


def main(argv=None) -> None:
    args = parse_args(argv)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    out_dir = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("  NDE Pipeline — Normalized Diffusion Entropy  (v%s)",
             __version__)
    log.info("=" * 60)
    log.info("Output directory: %s", out_dir)

    if args.skip_preproc and args.skip_dwi_preproc:
        _fail("--skip_preproc and --skip_dwi_preproc are mutually exclusive.")
    run_mode = ("eigenvalues" if args.skip_preproc else
                "dwi_preprocessed" if args.skip_dwi_preproc else "full")
    registration = args.registration or ("ants" if args.use_ants else None)

    # ------------------------------------------------------------------
    # Per-mode input checks and tool validation
    # ------------------------------------------------------------------
    if run_mode != "eigenvalues":
        need = [("--dwi", args.dwi), ("--bval", args.bval),
                ("--bvec", args.bvec)]
        if run_mode == "full":
            need.append(("--t1", args.t1))
        missing = [n for n, v in need if v is None]
        if missing:
            _fail("This run mode requires: %s", ", ".join(missing))
        try:
            dwi_info = nerves_dti.validate_dwi_inputs(
                args.dwi, args.bval, args.bvec)
        except ValueError as e:
            _fail("Invalid DWI input: %s", e)
        validate_tools(
            use_ants=args.use_ants, registration=registration,
            eddy_gpu=args.eddy_gpu, tensor_fit=args.tensor_fit,
            skip_dwi_preproc=(run_mode == "dwi_preprocessed"),
            has_t1=(args.t1 is not None))
    else:
        dwi_info = None

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
    snr_path  = args.snr
    extra = {"run_mode": run_mode}
    if dwi_info is not None:
        extra["acquisition"] = {k: dwi_info[k] for k in
                                ("n_volumes", "n_b0", "n_dw", "shells")}

    if run_mode in ("full", "dwi_preprocessed"):
        preproc_dir = out_dir / "preproc"
        if run_mode == "full":
            # ----------------------------------------------------------
            # Step 2: DWI preprocessing
            # ----------------------------------------------------------
            rpe_b0 = None if args.no_topup else args.rpe_b0
            preproc = preprocess_dwi(
                dwi=args.dwi, bval=args.bval, bvec=args.bvec,
                out_dir=preproc_dir, pe_dir=args.pe_dir, rpe_b0=rpe_b0,
                readout_time=args.readout_time, resume=args.resume,
                n_threads=args.n_threads, eddy_gpu=args.eddy_gpu,
            )
            mask_path = preproc["mask"]
            fit_dwi, fit_bvecs, fit_bvals = (
                preproc["dwi_eddy"], preproc["bvecs"], preproc["bvals"])
            ref_b0 = preproc["mean_b0"]
            ref_b0_brain = preproc["mean_b0_brain"]
            if snr_path is None and preproc.get("noise_map") is not None:
                try:
                    snr_info = nerves_dti.compute_snr_map(
                        ref_b0, preproc["noise_map"], mask_path,
                        preproc_dir / "snr_map.nii.gz")
                    snr_path = snr_info["path"]
                    extra["snr"] = snr_info["stats"]
                    log.info("Median b=0 SNR in mask: %.1f",
                             snr_info["stats"].get("snr_median", float("nan")))
                except ValueError as e:
                    log.warning("SNR map not computed: %s", e)
        else:
            # ----------------------------------------------------------
            # Preprocessed DWI supplied: mean b0, mask, then fit
            # ----------------------------------------------------------
            preproc_dir.mkdir(parents=True, exist_ok=True)
            fit_dwi, fit_bvecs, fit_bvals = args.dwi, args.bvec, args.bval
            ref_b0 = nerves_dti.mean_b0_numpy(
                args.dwi, args.bval, preproc_dir / "mean_b0.nii.gz")
            if mask_path is None:
                if not check_tool("bet"):
                    _fail("--skip_dwi_preproc needs --mask (FSL bet not "
                          "found to create one).")
                brain = preproc_dir / "mean_b0_brain"
                run(["bet", ref_b0, brain, "-m", "-f", "0.3"])
                mask_path = preproc_dir / "mean_b0_brain_mask.nii.gz"
            _validate_shapes({"--dwi (b0 grid)": ref_b0, "--mask": mask_path})
            ref_b0_brain = None
            if registration == "epi_reg" and args.t1 is not None:
                ref_b0_brain = preproc_dir / "mean_b0_masked.nii.gz"
                run(["fslmaths", ref_b0, "-mas", mask_path, ref_b0_brain])

        # --------------------------------------------------------------
        # Step 3: T1 preprocessing and registration (if a T1 is available)
        # --------------------------------------------------------------
        if args.t1 is not None:
            t1_results = preprocess_t1(
                t1=args.t1, dwi_ref=ref_b0, out_dir=out_dir / "t1",
                use_ants=args.use_ants, registration=registration,
                dwi_ref_brain=ref_b0_brain, resume=args.resume,
            )
            pve_csf = t1_results["pve_csf"]
            pve_gm  = t1_results["pve_gm"]
            pve_wm  = t1_results["pve_wm"]
        else:
            log.info("No T1 supplied: tissue-specific QC will be skipped.")

        # --------------------------------------------------------------
        # Optional shell selection, then Step 4: tensor estimation
        # --------------------------------------------------------------
        if args.max_b is not None:
            try:
                sel = nerves_dti.select_shells(
                    fit_dwi, fit_bvecs, fit_bvals, args.max_b,
                    out_dir / "preproc" / "shells")
            except ValueError as e:
                _fail("%s", e)
            fit_dwi, fit_bvecs, fit_bvals = (
                sel["dwi"], sel["bvecs"], sel["bvals"])
            extra["max_b"] = args.max_b
        tensor = estimate_tensor(
            dwi=fit_dwi, bvecs=fit_bvecs, bvals=fit_bvals, mask=mask_path,
            out_dir=out_dir / "tensor", wls=args.wls,
            method=args.tensor_fit, resume=args.resume,
        )
        l1_path = tensor["l1"]
        l2_path = tensor["l2"]
        l3_path = tensor["l3"]
        fa_path = tensor["fa"]
        md_path = tensor["md"]
        extra["tensor_fit"] = {"backend": args.tensor_fit, "wls": args.wls}

    else:
        # Validate that all required pre-computed files were supplied
        missing = [n for n, p in
                   [("--l1", l1_path), ("--l2", l2_path),
                    ("--l3", l3_path), ("--mask", mask_path)]
                   if p is None]
        if missing:
            log.error("--skip_preproc requires: %s", ", ".join(missing))
            sys.exit(1)
        _validate_shapes({"--l1": l1_path, "--l2": l2_path,
                          "--l3": l3_path, "--mask": mask_path,
                          "--fa": args.fa, "--md": args.md,
                          "--snr": args.snr,
                          "--pve_csf": pve_csf, "--pve_gm": pve_gm,
                          "--pve_wm": pve_wm})
        # Derive FA / MD paths if they exist alongside the eigenvalues
        # (assume dtifit naming convention: ..._L1.nii.gz -> ..._FA.nii.gz).
        # --fa/--md, if given, always take priority over auto-detection.
        # Guard on "L1" actually being in the filename first: str.replace()
        # silently no-ops otherwise, which previously made candidate_fa
        # equal to l1_path itself (since that file obviously exists) and
        # fed the eigenvalue map into run_qc() as if it were FA data.
        if args.fa is not None:
            fa_path = args.fa
        elif "L1" in l1_path.name:
            candidate_fa = l1_path.parent / l1_path.name.replace("L1", "FA")
            fa_path = candidate_fa if candidate_fa.exists() else None
        else:
            fa_path = None

        if args.md is not None:
            md_path = args.md
        elif "L1" in l1_path.name:
            candidate_md = l1_path.parent / l1_path.name.replace("L1", "MD")
            md_path = candidate_md if candidate_md.exists() else None
        else:
            md_path = None

        if fa_path is None:
            log.warning("FA map not found alongside L1 (expected dtifit "
                        "naming convention, e.g. ..._L1.nii.gz / "
                        "..._FA.nii.gz); QC scatter plot will be skipped. "
                        "Pass --fa (and --md) explicitly if your files use "
                        "different names.")

    # ------------------------------------------------------------------
    # Step 5: NDE computation (+ companion shape maps)
    # ------------------------------------------------------------------
    nde_dir = out_dir / "nde"
    nde_result = compute_nde(
        l1_path=l1_path,
        l2_path=l2_path,
        l3_path=l3_path,
        mask_path=mask_path,
        out_dir=nde_dir,
        snr_path=snr_path,
        min_snr=args.min_snr,
    )
    if not args.nde_only:
        shape = nerves_dti.compute_shape_metrics(
            l1_path, l2_path, l3_path, mask_path, out_dir / "shape",
            qc_flag_path=nde_result["qc_flag"])
        extra["shape_statistics"] = shape["stats"]
        log.info("Shape maps written: %s",
                 ", ".join(sorted(shape["paths"])))

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
        versions=collect_versions(),
        extra=extra,
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
                    "exceeds 5%. Review preprocessing.")
        sys.exit(2)


if __name__ == "__main__":
    main()
