"""Build the reproducibility notebook for the unconditional displacement study."""

from __future__ import annotations

import json
from pathlib import Path


def markdown_cell(text: str) -> dict:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": text.splitlines(keepends=True),
    }


def code_cell(code: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": code.splitlines(keepends=True),
    }


def main() -> None:
    study_dir = Path(__file__).resolve().parent.parent
    notebook_path = study_dir / "scripts" / "reproducibility_notebook.ipynb"

    cells = [
        markdown_cell(
            """# Unconditional Cavity Displacement in Dispersive cQED - Reproducibility Notebook

This notebook reproduces the main results of the unconditional-displacement study by loading the saved `unconditional_*` artifacts and, if requested, re-running the study scripts.

The workflow is intentionally split into:

1. environment and configuration
2. optional end-to-end rerun
3. artifact loading and headline summaries
4. figure review
5. validation spot checks
6. optimal-control waveform inspection

The main scientific takeaway is that the best simple, physically interpretable protocol is a short two-tone branch-compensated pulse, while the best overall mean fidelity on the broad state-test set comes from the constrained hardware-aware optimal-control waveform.
"""
        ),
        code_cell(
            """import sys
from pathlib import Path

cwd = Path.cwd()
candidates = [
    cwd,
    cwd / "studies" / "waveform_level_gate_realization_dispersive_cqed" / "scripts",
]
SCRIPT_DIR = next((path for path in candidates if (path / "common.py").exists()), None)
if SCRIPT_DIR is None:
    raise FileNotFoundError("Could not locate the study scripts directory from the current working directory.")

STUDY_DIR = SCRIPT_DIR.parent
ARTIFACTS_DIR = STUDY_DIR / "artifacts"
FIGURES_DIR = STUDY_DIR / "figures"
REPORT_DIR = STUDY_DIR / "report"
REPO_ROOT = STUDY_DIR.parent.parent
CQED_SIM_ROOT = REPO_ROOT.parent / "cQED_simulation"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(CQED_SIM_ROOT) not in sys.path:
    sys.path.insert(0, str(CQED_SIM_ROOT))

import common
common.apply_plot_style()

print(f"Study directory: {STUDY_DIR}")
print(f"Artifacts directory: {ARTIFACTS_DIR}")
print(f"Figures directory: {FIGURES_DIR}")
print(f"Report directory: {REPORT_DIR}")
"""
        ),
        markdown_cell(
            """## Tunable Parameters

All tunable parameters used by the notebook live in the next cell. The defaults mirror the saved study configuration. Set the rerun flags to `True` if you want to regenerate the artifacts instead of loading the existing ones.
"""
        ),
        code_cell(
            """# Study-scale knobs
N_TR = 3
N_CAV = 15
OPTIMAL_N_CAV = 12
DEFAULT_DT_NS = 0.5

# Physical parameters
CHI_MHZ = -2.84
CHI_PRIME_KHZ = -21.0
KERR_KHZ = -28.0

# Sweep grids
ALPHA_TARGETS = [0.5, 1.0, 1.5, 2.0]
SINGLE_TONE_DURATIONS_NS = [5.0, 10.0, 20.0, 40.0, 80.0, 160.0]
ECHO_TOTAL_DURATIONS_NS = [60.0, 80.0, 120.0, 160.0]
OPTIMAL_DURATIONS_NS = [40.0, 80.0]
CHI_SCALE_FACTORS = [0.5, 1.0, 1.5, 2.0]

# Validation states
QUBIT_STATES = ["g", "e", "plus_x", "plus_y"]
CAVITY_STATES = ["vacuum", "fock1", "fock2", "fock3", "coherent"]

# Notebook control flags
RERUN_MAIN_DRIVER = False
RERUN_VALIDATION = False
RERUN_OPTIMAL_WAVEFORM_PLOT = False

print("Configuration loaded.")
print(f"n_tr={N_TR}, n_cav={N_CAV}, optimal_n_cav={OPTIMAL_N_CAV}, dt={DEFAULT_DT_NS} ns")
print(f"chi/2pi={CHI_MHZ} MHz, chi'/2pi={CHI_PRIME_KHZ} kHz, K/2pi={KERR_KHZ} kHz")
"""
        ),
        markdown_cell(
            """## Optional End-to-End Rerun

Turn on any of the rerun flags in the configuration cell and execute the next cell if you want to regenerate artifacts or figures. The default path is load-only so the notebook stays lightweight.
"""
        ),
        code_cell(
            """if RERUN_MAIN_DRIVER:
    import unconditional_displacement_study
    unconditional_displacement_study.main()

if RERUN_VALIDATION:
    import unconditional_validation
    unconditional_validation.main()

if RERUN_OPTIMAL_WAVEFORM_PLOT:
    import plot_unconditional_optimal_waveform
    plot_unconditional_optimal_waveform.main()

print("Optional rerun cell finished.")
"""
        ),
        markdown_cell(
            """## Load Saved Artifacts

The helper `common.load_json` returns the artifact payload directly, so the notebook can work with the saved JSON without extra schema handling.
"""
        ),
        code_cell(
            """summary = common.load_json(ARTIFACTS_DIR / "unconditional_displacement_summary.json")
protocols = common.load_json(ARTIFACTS_DIR / "unconditional_protocol_comparison.json")
single_tone = common.load_json(ARTIFACTS_DIR / "unconditional_single_tone_summary.json")
two_tone = common.load_json(ARTIFACTS_DIR / "unconditional_two_tone_summary.json")
echo_data = common.load_json(ARTIFACTS_DIR / "unconditional_echo_summary.json")
optimal = common.load_json(ARTIFACTS_DIR / "unconditional_optimal_control_summary.json")
validation = common.load_json(ARTIFACTS_DIR / "unconditional_validation_spotcheck.json")

print("Loaded unconditional-displacement artifacts:")
for name in [
    "unconditional_displacement_summary.json",
    "unconditional_protocol_comparison.json",
    "unconditional_single_tone_summary.json",
    "unconditional_two_tone_summary.json",
    "unconditional_echo_summary.json",
    "unconditional_optimal_control_summary.json",
    "unconditional_validation_spotcheck.json",
]:
    print(" -", name)
"""
        ),
        markdown_cell(
            """## Headline Numerical Results

This cell prints the main benchmark points used throughout the report.
"""
        ),
        code_cell(
            """baseline = summary["baseline_square_80ns"]
fast = summary["best_fast_gaussian"]
twotone_best = summary["best_two_tone"]
echo_best = summary["best_echo"]
optimal_best = summary["best_optimal"]

print("Baseline naive square (80 ns, alpha=1):")
print(f"  delta_alpha = {baseline['delta_alpha']:.6f}")
print(f"  F_g = {baseline['g_fidelity']:.6f}, F_e = {baseline['e_fidelity']:.6f}, F_+x = {baseline['plus_x_fidelity']:.6f}")
print(f"  entanglement = {baseline['plus_x_entanglement_bits']:.6f} bits")

print("\\nBest fast Gaussian:")
print(f"  duration = {fast['duration_ns']:.1f} ns, bandwidth = {fast['bandwidth_mhz']:.1f} MHz")
print(f"  delta_alpha = {fast['delta_alpha']:.6f}, F_+x = {fast['plus_x_fidelity']:.6f}")

print("\\nBest two-tone:")
print(f"  duration = {twotone_best['duration_ns']:.1f} ns")
print(f"  delta_alpha = {twotone_best['delta_alpha']:.6e}, F_+x = {twotone_best['plus_x_fidelity']:.6f}")

print("\\nBest echoed case:")
print(f"  total duration = {echo_best['total_duration_ns']:.1f} ns")
print(f"  delta_alpha = {echo_best['delta_alpha']:.6f}, F_+x = {echo_best['plus_x_fidelity']:.6f}")

print("\\nBest constrained optimal-control case:")
print(f"  duration = {optimal_best['duration_ns']:.1f} ns")
print(f"  mean state fidelity = {optimal_best['full_metrics']['state_test_mean_fidelity']:.6f}")
print(f"  minimum state fidelity = {optimal_best['full_metrics']['state_test_min_fidelity']:.6f}")
print(f"  vacuum delta_alpha = {optimal_best['vacuum_metrics']['delta_alpha']:.6f}")
"""
        ),
        markdown_cell(
            """## Protocol Ranking Table

The next cell shows the common state-test-set ranking used in the report.
"""
        ),
        code_cell(
            """rows = []
for item in protocols["protocols"]:
    rows.append({
        "label": item["label"],
        "duration_ns": item["duration_ns"],
        "mean_fidelity": round(item["full_metrics"]["state_test_mean_fidelity"], 6),
        "min_fidelity": round(item["full_metrics"]["state_test_min_fidelity"], 6),
        "vacuum_delta_alpha": round(item["vacuum_metrics"]["delta_alpha"], 6),
        "vacuum_plus_x_fidelity": round(item["vacuum_metrics"]["plus_x_fidelity"], 6),
        "vacuum_plus_x_ent_bits": round(item["vacuum_metrics"]["plus_x_entanglement_bits"], 6),
        "complexity": item["complexity"],
    })

rows
"""
        ),
        markdown_cell(
            """## Figure Review

The report relies on the saved figure files. Displaying them here makes it easy to inspect the branch mismatch, entanglement, filter tradeoff, protocol summary, Wigner comparison, and optimal waveform directly from the notebook.
"""
        ),
        code_cell(
            """from IPython.display import Image, display

for name in [
    "unconditional_branch_mismatch.png",
    "unconditional_superposition_entanglement.png",
    "unconditional_filter_tradeoff.png",
    "unconditional_protocol_summary.png",
    "unconditional_wigner_comparison.png",
    "unconditional_optimal_waveform.png",
]:
    path = FIGURES_DIR / name
    print(f"Displaying {name}")
    display(Image(filename=str(path)))
"""
        ),
        markdown_cell(
            """## Validation Spot Checks

The study includes a direct unconditional-displacement validation artifact in addition to the broader shared-model convergence evidence from the earlier waveform-level study.
"""
        ),
        code_cell(
            """print("Ideal limit:")
ideal = validation["ideal_limit"]
print(f"  delta_alpha = {ideal['delta_alpha']:.6e}")
print(f"  F_+x = {ideal['plus_x_fidelity']:.6f}")

print("\\nn_cav sweep:")
for row in validation["n_cav_sweep"]:
    print(f"  n_cav={row['n_cav']}: delta_alpha={row['delta_alpha']:.6f}, F_+x={row['plus_x_fidelity']:.6f}")

print("\\ndt sweep:")
for row in validation["dt_sweep"]:
    print(f"  dt={row['dt_ns']:.2f} ns: delta_alpha={row['delta_alpha']:.6f}, F_+x={row['plus_x_fidelity']:.6f}")

print("\\nn_tr sweep:")
for row in validation["n_tr_sweep"]:
    print(f"  n_tr={row['n_tr']}: delta_alpha={row['delta_alpha']:.6f}, F_+x={row['plus_x_fidelity']:.6f}")
"""
        ),
        markdown_cell(
            """## Optimal-Control Waveform Details

The next cell exposes the actual optimized schedule used in the appendix, together with the main hardware metrics. This keeps the optimization result transparent rather than treating it as a black box.
"""
        ),
        code_cell(
            """best_case = max(optimal["cases"], key=lambda item: item["full_metrics"]["state_test_mean_fidelity"])
schedule = best_case["schedule_values"]
hardware = best_case["hardware_metrics"]

print("Held-sample schedule values (rad/s):")
print(schedule)

print("\\nHardware metrics:")
for key in [
    "physical_max_abs_amplitude",
    "physical_rms_amplitude",
    "physical_max_slew",
    "hardware_map_count",
]:
    print(f"  {key}: {hardware[key]}")
"""
        ),
        markdown_cell(
            """## Final Interpretation

The saved data support three practical conclusions:

1. a long naive cavity pulse should not be modeled as unconditional in this dispersive regime,
2. short two-tone compensation is the best simple experimental strategy here,
3. constrained optimal control gives the best broad-set mean fidelity but at the cost of interpretability.

If you change the configuration or rerun flags above, re-execute the notebook from the configuration cell onward so the plots and summaries stay consistent with the modified workflow.
"""
        ),
    ]

    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {
                "name": "python",
                "version": "3.12",
            },
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }

    notebook_path.write_text(json.dumps(notebook, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
