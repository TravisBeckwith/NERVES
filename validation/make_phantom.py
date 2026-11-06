#!/usr/bin/env python3
"""
validation/make_phantom.py
--------------------------
Synthetic DWI phantom with *known* tensors, for end-to-end testing and for
checking that the pipeline recovers NDE.

Layout (default 16x16x8): an anisotropic "white matter" block with fibres
along x, a "grey matter" slab (mildly anisotropic), a "CSF" corner
(isotropic, 3e-3) and background. Acquisition: n_b0 b=0 volumes + n_dirs
directions at b=1000 s/mm^2, Rician noise at the given b=0 SNR.

Writes (FSL conventions): dwi.nii.gz, dwi.bval, dwi.bvec (3xN), mask.nii.gz,
truth_L1/L2/L3.nii.gz, truth_labels.nii.gz (1=WM, 2=GM, 3=CSF).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from validation.noise_simulation import make_scheme  # noqa: E402

TISSUES = {1: (1.7, 0.35, 0.30), 2: (0.90, 0.80, 0.70), 3: (3.0, 3.0, 3.0)}


def make_phantom(out_dir: Path, snr: float = 30.0, seed: int = 0,
                 shape=(16, 16, 8), n_dirs: int = 30, n_b0: int = 3,
                 b: float = 1000.0) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    X, Y, Z = shape
    labels = np.zeros(shape, dtype=np.uint8)
    labels[2:X - 2, 2:Y - 2, 1:Z - 1] = 2                  # GM slab
    labels[4:X - 4, 4:Y - 4, 2:Z - 2] = 1                  # WM block inside
    labels[:3, :3, :3] = 3                                 # CSF corner
    mask = labels > 0

    bvals, bvecs = make_scheme(n_dirs, n_b0, b)
    # Truth eigenvalues (1e-3 mm^2/s); fibres along x (WM) or random (GM)
    lam = np.zeros(shape + (3,))
    for lab, ev in TISSUES.items():
        lam[labels == lab] = ev
    D = np.zeros(shape + (3, 3))
    idx = np.argwhere(mask)
    for (i, j, k) in idx:
        ev = lam[i, j, k] * 1e-3
        if labels[i, j, k] == 1:
            Q = np.eye(3)                                  # fibres along x
        else:
            Qm, R = np.linalg.qr(rng.standard_normal((3, 3)))
            Q = Qm * np.sign(np.diag(R))
        D[i, j, k] = Q @ np.diag(ev) @ Q.T
    q = np.einsum("mi,xyzij,mj->xyzm", bvecs, D, bvecs)
    sig = np.exp(-bvals * q)
    sig[~mask] = 0.0
    s = 1.0 / snr
    noisy = np.sqrt((sig + s * rng.standard_normal(sig.shape)) ** 2
                    + (s * rng.standard_normal(sig.shape)) ** 2)
    noisy *= 1000.0                                        # arbitrary S0
    aff = np.diag([2.0, 2.0, 2.0, 1.0])

    paths = {}

    def save(name, arr, dtype=np.float32):
        p = out_dir / name
        nib.save(nib.Nifti1Image(arr.astype(dtype), aff), str(p))
        paths[name.split(".")[0]] = p

    save("dwi.nii.gz", noisy)
    save("mask.nii.gz", mask, np.uint8)
    save("truth_labels.nii.gz", labels, np.uint8)
    for n in range(3):
        save(f"truth_L{n + 1}.nii.gz", lam[..., n] * 1e-3)
    np.savetxt(out_dir / "dwi.bval", bvals[None, :], fmt="%g")
    np.savetxt(out_dir / "dwi.bvec", bvecs.T, fmt="%.8g")
    paths["bval"], paths["bvec"] = out_dir / "dwi.bval", out_dir / "dwi.bvec"
    paths["dwi"], paths["mask"] = out_dir / "dwi.nii.gz", out_dir / "mask.nii.gz"
    return paths


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--snr", type=float, default=30.0)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    p = make_phantom(a.out, a.snr, a.seed)
    print("Phantom written to", a.out)
    return p


if __name__ == "__main__":
    main()
