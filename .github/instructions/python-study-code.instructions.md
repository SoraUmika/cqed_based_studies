---
description: "Python code conventions for cQED simulation study scripts. Enforces cqed_sim usage, figure saving, data serialization, and script structure standards."
applyTo: "studies/**/*.py"
---

# Python Study Code Conventions

## Script Structure
- Pin key simulation parameters (Hilbert space dims, time steps, frequencies) as UPPER_CASE named constants at the top of each script.
- Include a module docstring explaining the physical setup and parameters.
- Scripts must be self-contained and runnable from the study folder.

## cqed_sim Usage
- Always import from `cqed_sim` before writing custom simulation logic.
- Consult the API Reference before writing any simulation code.
- Document any gap in cqed_sim coverage in the study README under `## cqed_sim Gap Analysis`.

## Figure Saving
- Save every figure in **both** `.png` (300 dpi) and `.pdf` (vector) formats.
- Use `plt.savefig(fig_path, dpi=300, bbox_inches='tight')` for PNG.
- Use descriptive filenames: `<figure_description>.{png,pdf}`, not numbered names.
- Save to the study's `figures/` directory.
- Use colorblind-friendly palettes (e.g., `tab10`, seaborn `colorblind`).

## Data Serialization
- Save parameters and metadata as **JSON** in `artifacts/`.
- Save large arrays (waveforms, unitaries, sweep grids) as **NPZ** in `artifacts/` or `data/`.
- Save tabular data as **CSV** in `data/`.
- Include metadata fields in JSON artifacts: `study_name`, `date_created`, `description`, `parameters`.

## Import Order
- Standard library imports first, then third-party (numpy, scipy, matplotlib, qutip), then cqed_sim, then local modules.

## Red/Green Validation (Test-First)

Before writing simulation code, create `scripts/test_validation.py` encoding known constraints:
- **Limiting cases**: zero drive → no transfer, infinite detuning → decoupled
- **Conservation laws**: Tr(ρ) = 1, unitarity of U
- **Analytic benchmarks**: perturbation theory predictions, known closed-form results
- **Convergence**: Hilbert space dimension stability, time step stability

Workflow:
1. Copy `tools/templates/test_validation_template.py` to `scripts/test_validation.py`
2. Customize tests for your physics problem
3. Run `python scripts/test_validation.py --red` — confirm all tests FAIL (RED phase)
4. Implement the simulation
5. Run `python scripts/test_validation.py` — all tests must PASS (GREEN phase)
6. Reference test results quantitatively in the Validation section of the report

See the `red-green-validation` skill for detailed guidance.

## Acceleration
- Before starting any simulation expected to run longer than a few minutes, apply parallelization (`joblib.Parallel`, `multiprocessing.Pool`) or vectorization.
- Log wall-clock times for significant simulations in `IMPROVEMENTS.md` under `## Compute & Resource Notes`.
