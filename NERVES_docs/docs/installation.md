# Installation Guide

## System Requirements

NERVES has been tested on:
- Ubuntu 20.04 / 22.04
- macOS 12 (Monterey) and later
- CentOS 7 / Rocky Linux 8

Windows is not natively supported. WSL2 (Windows Subsystem for Linux) with
Ubuntu 22.04 is a viable alternative.

---

## 1. Install FSL

Follow the official FSL installation guide at:
https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation

Verify installation:
```bash
flirt -version
# Expected: FLIRT version 6.0 (or higher)
```

Ensure `$FSLDIR/bin` is on your `PATH` and `$FSLDIR` is set. The FSL
installer typically handles this automatically by appending to your
`.bashrc` or `.bash_profile`.

---

## 2. Install MRtrix3

**Option A — Package manager (recommended for Linux):**
```bash
sudo bash -c "$(curl -fsSL https://raw.githubusercontent.com/MRtrix3/mrtrix3/master/install_linux_packages.sh)"
```

**Option B — Conda:**
```bash
conda install -c mrtrix3 mrtrix3
```

**Option C — Build from source:**
```bash
git clone https://github.com/MRtrix3/mrtrix3.git
cd mrtrix3
./configure
./build
export PATH="$(pwd)/bin:$PATH"
```

Verify installation:
```bash
mrconvert --version
# Expected: mrconvert 3.0.x (or higher)
```

---

## 3. Install ANTs (Optional)

ANTs is only required if you use the `--use_ants` flag. Skip this step
if you are using the default FSL FLIRT registration.

**Option A — Pre-compiled binaries (recommended):**
Download from https://github.com/ANTsX/ANTs/releases and add the `bin/`
directory to your `PATH`.

**Option B — Conda:**
```bash
conda install -c aramislab ants
```

Verify installation:
```bash
antsRegistrationSyNQuick.sh --help
```

---

## 4. Install NERVES

```bash
# Clone the repository
git clone https://github.com/yourusername/NERVES.git
cd NERVES

# Create virtual environment
python -m venv nerves_env
source nerves_env/bin/activate

# Install Python dependencies
pip install -r requirements.txt
```

### Python Dependencies

| Package | Minimum version | Purpose |
|---|---|---|
| numpy | 1.21.0 | Core array operations for NDE computation |
| nibabel | 3.2.0 | NIfTI file I/O |
| scipy | 1.7.0 | Statistical functions, KDE |
| matplotlib | 3.4.0 | QC figure generation |
| seaborn | 0.11.0 | Plot styling |
| pandas | 1.3.0 | Report tabulation |

---

## 5. Verify Full Installation

```bash
python nde_pipeline.py --help
```

This will check that all required FSL and MRtrix3 tools are available on
your `PATH` before printing the argument list.

To run a quick sanity check with a synthetic dataset:
```bash
python tests/test_nde_math.py
```

Expected output:
```
Config A: NDE=0.8650  [PASS]
Config B: NDE=0.8173  [PASS]
Config C: NDE=0.9464  [PASS]
Config D: NDE=0.9207  [PASS]
Isotropic: NDE=1.0000  [PASS]
All tests passed.
```

---

## Troubleshooting Installation

**`bet: command not found`**
FSL is not on your PATH. Run `source $FSLDIR/etc/fslconf/fsl.sh` or add
`$FSLDIR/bin` to your PATH manually.

**`mrconvert: command not found`**
MRtrix3 bin directory is not on PATH. Add it with:
`export PATH=/path/to/mrtrix3/bin:$PATH`

**`ImportError: No module named nibabel`**
Your virtual environment is not activated, or `pip install -r requirements.txt`
was not run. Activate the environment and reinstall.

**Python version error**
NERVES requires Python ≥ 3.8. Check with `python --version`. On systems
where Python 3 is invoked as `python3`, use `python3 nde_pipeline.py`.
