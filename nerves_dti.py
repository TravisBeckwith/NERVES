#!/usr/bin/env python3
"""
nerves_dti.py
=============
Tool-free (pure NumPy / nibabel) building blocks used by ``nde_pipeline.py``:

* b-value / b-vector reading and validation, shell detection and selection
* a native diffusion-tensor fit (OLS or weighted least squares) so the tensor
  stage -- and therefore a complete, testable end-to-end run -- does not
  require FSL
* tensor-shape metrics that accompany NDE (tensor mode, Westin indices)
* a simple SNR map from a noise-sigma map

Units: b-values in s/mm^2, diffusivities in mm^2/s (as FSL ``dtifit``).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import nibabel as nib
import numpy as np

log = logging.getLogger("nde_pipeline")

#: b-values at or below this are treated as b = 0
B0_THRESHOLD = 50.0
#: DTI assumes Gaussian diffusion; above roughly this b-value kurtosis biases
#: the tensor fit and NDE/FA/MD should be interpreted with care.
MAX_RECOMMENDED_B = 1500.0


# ---------------------------------------------------------------------------
# b-values / b-vectors
# ---------------------------------------------------------------------------

def read_bvals(path: Path) -> np.ndarray:
    """Read an FSL bval file (whitespace separated) -> 1-D float array."""
    return np.atleast_1d(np.loadtxt(str(path), dtype=np.float64)).ravel()


def read_bvecs(path: Path) -> np.ndarray:
    """Read an FSL bvec file -> (N, 3). Accepts 3xN (FSL) or Nx3 layouts."""
    arr = np.atleast_2d(np.loadtxt(str(path), dtype=np.float64))
    if arr.shape[0] == 3:
        return arr.T
    if arr.shape[1] == 3:
        return arr
    raise ValueError(f"bvec file {path} must be 3xN or Nx3, got {arr.shape}")


def shells(bvals: np.ndarray, tol: float = B0_THRESHOLD) -> list:
    """Cluster b-values into shells; returns shell centres (b=0 first)."""
    out = []
    group = []
    for b in np.sort(np.asarray(bvals, dtype=float)):
        if group and b - group[0] > tol:
            out.append(float(np.mean(group)))
            group = []
        group.append(b)
    if group:
        out.append(float(np.mean(group)))
    return out


def volumes_up_to(bvals: np.ndarray, max_b: float,
                  tol: float = B0_THRESHOLD) -> np.ndarray:
    """Indices of volumes with b <= max_b (+ tolerance), including b=0."""
    return np.where(np.asarray(bvals) <= max_b + tol)[0]


def validate_dwi_inputs(dwi_path: Path, bval_path: Path, bvec_path: Path,
                        b0_threshold: float = B0_THRESHOLD) -> dict:
    """Sanity-check DWI + gradient table. Raises ``ValueError`` on problems.

    Checks: file existence, 4-D image, bval/bvec length equals number of
    volumes, at least one b=0 volume, >= 6 diffusion-weighted directions,
    unit-length non-b0 gradient vectors, and warns for b-values outside the
    range where the Gaussian tensor model is appropriate.
    """
    for p in (dwi_path, bval_path, bvec_path):
        if not Path(p).exists():
            raise ValueError(f"Input file not found: {p}")
    shape = nib.load(str(dwi_path)).shape
    if len(shape) != 4:
        raise ValueError(f"DWI must be 4-D, got shape {shape}")
    n_vols = shape[3]
    bvals = read_bvals(bval_path)
    bvecs = read_bvecs(bvec_path)
    if len(bvals) != n_vols or len(bvecs) != n_vols:
        raise ValueError(
            f"Gradient table does not match the image: {n_vols} volumes, "
            f"{len(bvals)} b-values, {len(bvecs)} b-vectors")
    is_b0 = bvals <= b0_threshold
    n_b0 = int(is_b0.sum())
    n_dw = int((~is_b0).sum())
    if n_b0 < 1:
        raise ValueError("No b=0 volumes found (b <= %g)." % b0_threshold)
    if n_dw < 6:
        raise ValueError(f"Need >= 6 diffusion-weighted volumes, found {n_dw}.")
    norms = np.linalg.norm(bvecs[~is_b0], axis=1)
    if np.any(np.abs(norms - 1.0) > 0.1):
        raise ValueError(
            "Some diffusion-weighted b-vectors are not unit length "
            f"(norm range {norms.min():.2f}-{norms.max():.2f}); check the "
            "bvec file (rows/columns, scaling).")
    sh = shells(bvals)
    warnings = check_bvalue_range(bvals)
    for w in warnings:
        log.warning(w)
    return {"n_volumes": n_vols, "n_b0": n_b0, "n_dw": n_dw,
            "shells": sh, "warnings": warnings}


def check_bvalue_range(bvals: np.ndarray) -> list:
    """Warnings about b-values / shells unsuitable for a plain DTI fit."""
    msgs = []
    sh = [s for s in shells(bvals) if s > B0_THRESHOLD]
    if sh and max(sh) > MAX_RECOMMENDED_B:
        msgs.append(
            f"Highest b-value is ~{max(sh):.0f} s/mm^2 (> "
            f"{MAX_RECOMMENDED_B:.0f}): kurtosis biases the DTI fit. "
            "Consider --max_b 1000 (or 1500) to fit a single lower shell.")
    if len(sh) > 1:
        msgs.append(
            f"Multi-shell data (shells ~{', '.join(f'{s:.0f}' for s in sh)}): "
            "all shells are used unless --max_b is given. A single-shell "
            "tensor fit on mixed shells is biased by non-Gaussian diffusion.")
    return msgs


def select_shells(dwi_path: Path, bvecs_path: Path, bvals_path: Path,
                  max_b: float, out_dir: Path) -> dict:
    """Write a copy of the DWI containing only volumes with b <= max_b.

    Returns dict(dwi, bvecs, bvals); the inputs are returned unchanged when
    every volume already satisfies the limit.
    """
    bvals = read_bvals(bvals_path)
    idx = volumes_up_to(bvals, max_b)
    if len(idx) == len(bvals):
        return {"dwi": dwi_path, "bvecs": bvecs_path, "bvals": bvals_path}
    n_dw = int((bvals[idx] > B0_THRESHOLD).sum())
    if n_dw < 6:
        raise ValueError(f"--max_b {max_b:g} leaves only {n_dw} diffusion-"
                         "weighted volumes (need >= 6).")
    out_dir.mkdir(parents=True, exist_ok=True)
    img = nib.load(str(dwi_path))
    data = np.asanyarray(img.dataobj)[..., idx]
    hdr = img.header.copy()
    out_dwi = out_dir / "dwi_selected.nii.gz"
    nib.save(nib.Nifti1Image(data, img.affine, hdr), str(out_dwi))
    bvecs = read_bvecs(bvecs_path)
    out_bval = out_dir / "dwi_selected.bval"
    out_bvec = out_dir / "dwi_selected.bvec"
    np.savetxt(out_bval, bvals[idx][None, :], fmt="%.6g")
    np.savetxt(out_bvec, bvecs[idx].T, fmt="%.8g")
    log.info("Shell selection: kept %d of %d volumes (b <= %g).",
             len(idx), len(bvals), max_b)
    return {"dwi": out_dwi, "bvecs": out_bvec, "bvals": out_bval}


def mean_b0_numpy(dwi_path: Path, bvals_path: Path, out_path: Path,
                  b0_threshold: float = B0_THRESHOLD) -> Path:
    """Mean of the b=0 volumes (no MRtrix needed)."""
    bvals = read_bvals(bvals_path)
    idx = np.where(bvals <= b0_threshold)[0]
    img = nib.load(str(dwi_path))
    mean = np.asanyarray(img.dataobj)[..., idx].astype(np.float64).mean(axis=-1)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    hdr = img.header.copy()
    hdr.set_data_dtype(np.float32)
    nib.save(nib.Nifti1Image(mean.astype(np.float32), img.affine, hdr),
             str(out_path))
    return out_path


# ---------------------------------------------------------------------------
# Native tensor fit
# ---------------------------------------------------------------------------

def design_matrix(bvals: np.ndarray, bvecs: np.ndarray) -> np.ndarray:
    """Log-linear DTI design matrix, columns [lnS0, Dxx, Dyy, Dzz, Dxy, Dxz, Dyz]."""
    b = np.asarray(bvals, dtype=np.float64)
    g = np.asarray(bvecs, dtype=np.float64)
    return np.column_stack([
        np.ones_like(b),
        -b * g[:, 0] ** 2, -b * g[:, 1] ** 2, -b * g[:, 2] ** 2,
        -2 * b * g[:, 0] * g[:, 1],
        -2 * b * g[:, 0] * g[:, 2],
        -2 * b * g[:, 1] * g[:, 2],
    ])


def fit_tensor_signals(signals: np.ndarray, bvals: np.ndarray,
                       bvecs: np.ndarray, method: str = "wls",
                       chunk: int = 50_000) -> np.ndarray:
    """Fit tensors to ``signals`` (N voxels, M volumes).

    ``method``: ``"ols"`` (log-linear least squares, as dtifit by default) or
    ``"wls"`` (weighted LS with weights = fitted signal^2, Salvador et al.
    2005 / Veraart et al. 2013), which is less biased for noisy data.

    Returns the 6 unique tensor elements per voxel as (N, 7) = [lnS0, Dxx,
    Dyy, Dzz, Dxy, Dxz, Dyz].
    """
    if method not in ("ols", "wls"):
        raise ValueError("method must be 'ols' or 'wls'")
    B = design_matrix(bvals, bvecs)
    P = np.linalg.pinv(B)
    out = np.empty((signals.shape[0], 7), dtype=np.float64)
    eps = 1e-12
    for s in range(0, signals.shape[0], chunk):
        S = np.maximum(np.asarray(signals[s:s + chunk], dtype=np.float64), eps)
        y = np.log(S)
        x = y @ P.T
        if method == "wls":
            w = np.exp(2.0 * (x @ B.T))                      # fitted S^2
            A = np.einsum("nm,mi,mj->nij", w, B, B)
            A += 1e-12 * np.eye(7)[None]
            rhs = np.einsum("nm,mi,nm->ni", w, B, y)
            try:
                x = np.linalg.solve(A, rhs[..., None])[..., 0]
            except np.linalg.LinAlgError:
                pass                                          # keep OLS
        out[s:s + chunk] = x
    return out


def tensor_eigenvalues(x: np.ndarray) -> np.ndarray:
    """(N, 7) fit output -> (N, 3) eigenvalues, descending."""
    n = x.shape[0]
    D = np.zeros((n, 3, 3))
    D[:, 0, 0], D[:, 1, 1], D[:, 2, 2] = x[:, 1], x[:, 2], x[:, 3]
    D[:, 0, 1] = D[:, 1, 0] = x[:, 4]
    D[:, 0, 2] = D[:, 2, 0] = x[:, 5]
    D[:, 1, 2] = D[:, 2, 1] = x[:, 6]
    return np.linalg.eigvalsh(D)[:, ::-1]


def fa_md(eig: np.ndarray):
    """FA and MD from (N, 3) eigenvalues."""
    md = eig.mean(axis=1)
    num = 1.5 * ((eig - md[:, None]) ** 2).sum(axis=1)
    den = (eig ** 2).sum(axis=1)
    fa = np.sqrt(np.divide(num, den, out=np.zeros_like(num), where=den > 0))
    return np.clip(fa, 0.0, 1.0), md


def estimate_tensor_numpy(dwi: Path, bvecs: Path, bvals: Path, mask: Path,
                          out_dir: Path, method: str = "ols") -> dict:
    """FSL-free equivalent of ``estimate_tensor``: fits the tensor in NumPy
    and writes dti_L1/L2/L3/FA/MD/RD with the same names as dtifit."""
    out_dir.mkdir(parents=True, exist_ok=True)
    img = nib.load(str(dwi))
    m_img = nib.load(str(mask))
    msk = np.asarray(m_img.dataobj) > 0
    if msk.shape != img.shape[:3]:
        raise ValueError(f"Mask shape {msk.shape} != DWI grid {img.shape[:3]}")
    data = np.asarray(img.dataobj, dtype=np.float32)[msk]       # (N, M)
    log.info("Native tensor fit (%s) on %d voxels, %d volumes",
             method.upper(), data.shape[0], data.shape[1])
    x = fit_tensor_signals(data, read_bvals(bvals), read_bvecs(bvecs), method)
    eig = tensor_eigenvalues(x)
    fa, md = fa_md(eig)
    rd = eig[:, 1:].mean(axis=1)

    hdr = img.header.copy()
    hdr.set_data_dtype(np.float32)
    shape3 = img.shape[:3]
    aff = img.affine

    def _save(name, vals):
        vol = np.zeros(shape3, dtype=np.float32)
        vol[msk] = vals
        p = out_dir / f"dti_{name}.nii.gz"
        h = hdr.copy()
        h.set_data_shape(shape3)
        nib.save(nib.Nifti1Image(vol, aff, h), str(p))
        return p

    paths = {
        "l1": _save("L1", eig[:, 0]), "l2": _save("L2", eig[:, 1]),
        "l3": _save("L3", eig[:, 2]), "fa": _save("FA", fa),
        "md": _save("MD", md), "rd": _save("RD", rd),
    }
    return paths


# ---------------------------------------------------------------------------
# Eigenvalue / tensor consistency checks
# ---------------------------------------------------------------------------

def ordering_violations(l1: np.ndarray, l2: np.ndarray, l3: np.ndarray,
                        mask: np.ndarray) -> dict:
    """Counts of voxels violating L1 >= L2 and L2 >= L3 inside ``mask``."""
    m = np.asarray(mask, dtype=bool)
    return {
        "n_mask": int(m.sum()),
        "l1_lt_l2": int(np.sum(l1[m] < l2[m])),
        "l2_lt_l3": int(np.sum(l2[m] < l3[m])),
    }


def trace_consistency(tensor_path: Path, l1: np.ndarray, l2: np.ndarray,
                      l3: np.ndarray, mask: np.ndarray) -> Optional[float]:
    """Median relative |trace(D) - (L1+L2+L3)| inside the mask.

    ``tensor_path`` is dtifit's ``--save_tensor`` output (Dxx, Dxy, Dxz, Dyy,
    Dyz, Dzz). Returns None if the file is missing or malformed. A large
    value means the eigenvalue and tensor files do not belong together.
    """
    tensor_path = Path(tensor_path)
    if not tensor_path.exists():
        return None
    t = np.asarray(nib.load(str(tensor_path)).dataobj, dtype=np.float64)
    if t.ndim != 4 or t.shape[3] != 6 or t.shape[:3] != l1.shape:
        return None
    m = np.asarray(mask, dtype=bool)
    tr = t[..., 0] + t[..., 3] + t[..., 5]
    esum = l1 + l2 + l3
    denom = np.maximum(np.abs(esum[m]), 1e-12)
    return float(np.median(np.abs(tr[m] - esum[m]) / denom)) if m.any() else None


# ---------------------------------------------------------------------------
# NDE
# ---------------------------------------------------------------------------

def nde_from_eigenvalues(eig: np.ndarray) -> np.ndarray:
    """Normalised Shannon entropy of eigenvalue proportions, (N, 3) -> (N,).

    Single definition used by the pipeline, the validation tools and the
    tests. Eigenvalues must be > 0 (the caller masks invalid voxels).
    0 * ln(0) is taken as 0.
    """
    eig = np.asarray(eig, dtype=np.float64)
    p = eig / eig.sum(axis=-1, keepdims=True)
    log_p = np.where(p > 0.0, np.log(np.where(p > 0.0, p, 1.0)), 0.0)
    return np.clip(-np.sum(p * log_p, axis=-1) / np.log(3.0), 0.0, 1.0)


# ---------------------------------------------------------------------------
# Shape metrics
# ---------------------------------------------------------------------------

def tensor_mode(eig: np.ndarray) -> np.ndarray:
    """Tensor mode (Ennis & Kindlmann 2006), in [-1, 1]; (N, 3) -> (N,).

    +1 = linear/prolate, 0 = orthotropic, -1 = planar/oblate. NaN for
    isotropic tensors (mode undefined).
    """
    eig = np.asarray(eig, dtype=np.float64)
    d = eig - eig.mean(axis=1, keepdims=True)
    nrm = np.linalg.norm(d, axis=1)
    ok = nrm > 1e-12 * np.maximum(np.abs(eig).max(axis=1), 1e-300)
    mode = np.full(eig.shape[0], np.nan)
    dn = d[ok] / nrm[ok, None]
    mode[ok] = np.clip(3.0 * np.sqrt(6.0) * np.prod(dn, axis=1), -1.0, 1.0)
    return mode


def westin_indices(eig: np.ndarray):
    """Westin linear/planar/spherical coefficients (cl, cp, cs), (N,3)->3x(N,)."""
    l1 = np.maximum(eig[:, 0], 1e-300)
    return ((eig[:, 0] - eig[:, 1]) / l1, (eig[:, 1] - eig[:, 2]) / l1,
            eig[:, 2] / l1)


def compute_shape_metrics(l1_path: Path, l2_path: Path, l3_path: Path,
                          mask_path: Path, out_dir: Path,
                          qc_flag_path: Optional[Path] = None) -> dict:
    """Write tensor mode and Westin cl/cp/cs maps (NaN where undefined).

    Voxels flagged in ``qc_flag_path`` (non-physical / low SNR) are NaN, so
    these maps are defined on exactly the same voxels as the NDE map.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    img = nib.load(str(l1_path))
    l1 = np.asarray(img.dataobj, dtype=np.float64)
    l2 = np.asarray(nib.load(str(l2_path)).dataobj, dtype=np.float64)
    l3 = np.asarray(nib.load(str(l3_path)).dataobj, dtype=np.float64)
    mask = np.asarray(nib.load(str(mask_path)).dataobj) > 0
    valid = mask & (l1 > 0) & (l2 > 0) & (l3 > 0)
    if qc_flag_path is not None:
        flag = np.asarray(nib.load(str(qc_flag_path)).dataobj)
        valid &= flag == 0
    eig = np.stack([l1[valid], l2[valid], l3[valid]], axis=1)
    eig = np.sort(eig, axis=1)[:, ::-1]
    mode = tensor_mode(eig)
    cl, cp, cs = westin_indices(eig)

    hdr = img.header.copy()
    hdr.set_data_dtype(np.float32)
    paths, stats = {}, {}
    for name, vals in (("mode", mode), ("cl", cl), ("cp", cp), ("cs", cs)):
        vol = np.zeros(l1.shape, dtype=np.float32)
        vol[mask] = np.nan
        vol[valid] = vals
        p = out_dir / f"{name}.nii.gz"
        nib.save(nib.Nifti1Image(vol, img.affine, hdr), str(p))
        paths[name] = p
        finite = vals[np.isfinite(vals)]
        if finite.size:
            stats[f"{name}_mean"] = float(finite.mean())
            stats[f"{name}_std"] = float(finite.std())
    return {"paths": paths, "stats": stats}


# ---------------------------------------------------------------------------
# SNR
# ---------------------------------------------------------------------------

def compute_snr_map(mean_b0_path: Path, noise_path: Path, mask_path: Path,
                    out_path: Path) -> dict:
    """SNR = mean b=0 signal / noise sigma (e.g. dwidenoise ``-noise`` map).

    Approximate: the sigma map lives in the raw acquisition space while the
    mean b=0 may come from motion-corrected data. Written inside the mask;
    voxels with sigma <= 0 are set to 0.
    """
    b0_img = nib.load(str(mean_b0_path))
    b0 = np.asarray(b0_img.dataobj, dtype=np.float64)
    sig = np.asarray(nib.load(str(noise_path)).dataobj, dtype=np.float64)
    mask = np.asarray(nib.load(str(mask_path)).dataobj) > 0
    if not (b0.shape == sig.shape == mask.shape):
        raise ValueError(f"SNR inputs differ in shape: b0 {b0.shape}, "
                         f"noise {sig.shape}, mask {mask.shape}")
    snr = np.zeros_like(b0)
    ok = mask & (sig > 0)
    snr[ok] = b0[ok] / sig[ok]
    hdr = b0_img.header.copy()
    hdr.set_data_dtype(np.float32)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(snr.astype(np.float32), b0_img.affine, hdr),
             str(out_path))
    vals = snr[ok]
    stats = {}
    if vals.size:
        stats = {"snr_median": float(np.median(vals)),
                 "snr_p5": float(np.percentile(vals, 5)),
                 "snr_p95": float(np.percentile(vals, 95))}
    return {"path": out_path, "stats": stats}
