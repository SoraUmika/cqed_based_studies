from __future__ import annotations

from pathlib import Path

import nbformat as nbf


NOTEBOOK_PATH = Path(__file__).resolve().parent / "reproducibility_notebook.ipynb"


def build_notebook() -> None:
    nb = nbf.v4.new_notebook()
    cells = [
        nbf.v4.new_markdown_cell(
            """# Reproducibility Notebook: Microwave Component Thermalization in cQED

This notebook reproduces the main results for the study of hot microwave components in a dispersive cQED platform. It is written to support two workflows:

1. Full reproduction: rerun the four simulation thrusts and the validation checks with tunable parameters exposed below.
2. Fast inspection: load the saved JSON artifacts and figures that were already generated in the study directory.

The quantum model resolves effective bath occupation, qubit heating, dephasing, multimode back-action, and intrinsic Lindblad response times. It does **not** model macroscopic thermal transport, boundary resistance, or component heat-flow delays; those must be coupled in as a separate thermal layer.
"""
        ),
        nbf.v4.new_code_cell(
            """from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def locate_scripts_dir() -> Path:
    candidates = [
        Path.cwd(),
        Path.cwd() / "scripts",
        Path.cwd() / "studies" / "microwave_component_thermalization_cqed" / "scripts",
    ]
    for candidate in candidates:
        if (candidate / "run_study.py").exists():
            return candidate.resolve()
    raise FileNotFoundError("Could not locate the study scripts directory from the current working directory.")


SCRIPTS_DIR = locate_scripts_dir()
STUDY_ROOT = SCRIPTS_DIR.parent
REPO_ROOT = STUDY_ROOT.parent.parent
DATA_DIR = STUDY_ROOT / "data"
ARTIFACTS_DIR = STUDY_ROOT / "artifacts"
FIGURES_DIR = STUDY_ROOT / "figures"

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from common import MultimodeStudyConfig, SingleModeStudyConfig  # noqa: E402
from run_study import run_full_study  # noqa: E402
from validate_results import run_validation  # noqa: E402


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


print(f"Scripts directory: {SCRIPTS_DIR}")
print(f"Study root:        {STUDY_ROOT}")
"""
        ),
        nbf.v4.new_code_cell(
            """# Tunable parameters for end-to-end reruns.
# Set RUN_FULL_REPRODUCTION = False if you only want to inspect the saved artifacts.

RUN_FULL_REPRODUCTION = True
SAVE_OUTPUTS = True

single_config = replace(
    SingleModeStudyConfig(),
    n_cav_steady=30,
    n_cav_dynamic=12,
    n_tr=3,
    temperature_grid=(0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.50, 1.0, 2.0),
    dephasing_temperatures=(0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.50),
    cold_temperature=0.05,
    transient_temperatures=(0.20, 0.50),
    coherence_target=20.0e-6,
)

multimode_config = replace(
    MultimodeStudyConfig(),
    n_readout=6,
    n_storage=6,
    hot_storage_temperature=0.35,
    cold_readout_temperature=0.05,
    detuning_grid_mhz=(-240.0, -180.0, -120.0, -60.0, 0.0, 60.0, 120.0, 180.0, 240.0),
    coupling_grid_mhz=(0.0, 1.5, 3.0, 4.5, 6.0, 7.5, 9.0, 10.5, 12.0),
)

print("Single-mode config:", single_config)
print("Multimode config:", multimode_config)
"""
        ),
        nbf.v4.new_code_cell(
            """if RUN_FULL_REPRODUCTION:
    study_results = run_full_study(single_config=single_config, multimode_config=multimode_config, save_outputs=SAVE_OUTPUTS)
    validation_results = run_validation(single=single_config, multi=multimode_config, save_output=SAVE_OUTPUTS)
else:
    study_results = None
    validation_results = None

summary = load_json(DATA_DIR / "study_summary.json")
thermometer = load_json(ARTIFACTS_DIR / "thermometer_summary.json")
dephasing = load_json(ARTIFACTS_DIR / "dephasing_summary.json")
multimode = load_json(ARTIFACTS_DIR / "multimode_summary.json")
transient = load_json(ARTIFACTS_DIR / "transient_summary.json")
validation = load_json(ARTIFACTS_DIR / "validation_summary.json")

print(json.dumps(summary, indent=2))
print("All validation checks passed:", validation["passed_all"])
"""
        ),
        nbf.v4.new_code_cell(
            """selected_temperatures = [0.05, 0.10, 0.20, 0.50, 2.0]
temperature_to_index = {round(temp, 2): index for index, temp in enumerate(thermometer["temperatures_K"])}

print("Thermometer calibration points")
for temperature in selected_temperatures:
    index = temperature_to_index[round(temperature, 2)]
    print(
        f"T = {temperature:>4.2f} K | "
        f"n_th = {thermometer['nth_values'][index]:>7.4f} | "
        f"n_cav = {thermometer['dispersive_cavity_occupation'][index]:>7.4f} | "
        f"P_e = {thermometer['dressed_qubit_excited_population'][index]:>7.4f} | "
        f"width = {thermometer['spectroscopy_width_MHz'][index]:>6.3f} MHz"
    )

print()
print("Dephasing checkpoints")
for temperature in [0.03, 0.05, 0.07, 0.10, 0.15, 0.20]:
    index = dephasing["temperatures_K"].index(temperature)
    print(
        f"T = {temperature:>4.2f} K | "
        f"gamma_phi = {dephasing['gamma_pure_per_s'][index]:>10.1f} 1/s | "
        f"gamma_up = {dephasing['gamma_up_per_s'][index]:>10.1f} 1/s | "
        f"T2 = {1.0e6 * dephasing['t2_total_s'][index]:>7.2f} us"
    )

print()
danger = summary["multimode"]["most_dangerous_point"]
print("Multimode safe fraction:", summary["multimode"]["safe_fraction"])
print("Most dangerous multimode point:", danger)

print()
for key, value in summary["transient"].items():
    print(
        f"{key}: cavity tau = {value['cavity_response_tau_us']:.3f} us, "
        f"qubit tau = {value['qubit_response_tau_us']:.3f} us"
    )
"""
        ),
        nbf.v4.new_code_cell(
            """print("Validation audit")
for check in validation["checks"]:
    status = "PASS" if check["passed"] else "FAIL"
    print(f"[{status}] {check['name']}: {check['criterion']}")
"""
        ),
        nbf.v4.new_code_cell(
            """figure_stems = summary["figure_stems"]
figure, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

for axis, stem in zip(axes.ravel(), figure_stems):
    image = plt.imread(FIGURES_DIR / f"{stem}.png")
    axis.imshow(image)
    axis.set_title(stem.replace("_", " "))
    axis.axis("off")

plt.show()
"""
        ),
        nbf.v4.new_markdown_cell(
            """## How to extend this notebook

- Change the truncations, temperature grids, or safe-regime thresholds in the parameter cell and rerun the notebook.
- Replace the representative multimode frequencies, linewidths, or couplings with measured device values once those are available.
- Couple a separate thermal RC or finite-element model into the notebook by generating a time-dependent effective bath occupation and feeding that into the existing cQED layer.
"""
        ),
    ]

    nb["cells"] = cells
    nb["metadata"] = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {
            "name": "python",
            "version": "3.12",
        },
    }
    NOTEBOOK_PATH.write_text(nbf.writes(nb), encoding="utf-8")


if __name__ == "__main__":
    build_notebook()
