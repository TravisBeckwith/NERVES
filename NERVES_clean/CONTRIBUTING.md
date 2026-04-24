# Contributing to NERVES

Thank you for your interest in contributing. This document covers how to
report bugs, suggest features, and submit code.

---

## Reporting Bugs

Before opening an issue, please:
1. Check the [Troubleshooting guide](docs/troubleshooting.md)
2. Search existing [Issues](../../issues) for duplicates

When opening a bug report, use the **Bug report** template and include:
- Full command used
- Contents of `nde_pipeline_report.json` (if generated)
- Terminal output with `--verbose`
- FSL, MRtrix3, and Python versions

---

## Suggesting Features

Open a **Feature request** issue describing:
- The use case this feature addresses
- Whether it relates to preprocessing, NDE computation, QC, or documentation
- Any relevant literature

---

## Contributing Code

### Setup

```bash
git clone https://github.com/yourusername/NERVES.git
cd NERVES
python -m venv nerves_dev
source nerves_dev/bin/activate
pip install -r requirements.txt
pip install pytest flake8
```

### Workflow

1. Fork the repository and create a branch from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   ```
2. Make your changes
3. Run the test suite:
   ```bash
   pytest tests/
   ```
4. Check code style:
   ```bash
   flake8 nde_pipeline.py --max-line-length=100
   ```
5. Open a pull request against `main` with a clear description of
   what changed and why

### Code Style

- Follow PEP 8 with a maximum line length of 100 characters
- All new functions must include a docstring
- Mathematical operations should include a comment referencing the
  equation in [docs/theory.md](docs/theory.md) where applicable
- Use `np.where` rather than masked arrays for NIfTI operations

### Adding Tests

New tests go in `tests/`. The test suite must be runnable without
FSL or MRtrix3 installed — tests that require external tools should
be marked with `@pytest.mark.integration` and skipped by default.

The core NDE math is tested in `tests/test_nde_math.py` and does not
require any external dependencies beyond NumPy.

---

## Documentation

Documentation lives in `docs/`. Edits to the theory, parameters, QC
guide, or troubleshooting pages are always welcome. For substantial
changes to the mathematical content in `docs/theory.md`, please link
to or attach the relevant reference.
