"""Build the reproducibility notebook for the gf-sideband cooling study."""

from __future__ import annotations

from pathlib import Path

import nbformat as nbf


SCRIPT_DIR = Path(__file__).resolve().parent
NOTEBOOK_PATH = SCRIPT_DIR / "reproducibility_notebook.ipynb"


def md(text: str):
    return nbf.v4.new_markdown_cell(text)


def code(text: str):
    return nbf.v4.new_code_cell(text)


def main() -> None:
    nb = nbf.v4.new_notebook()
    nb["metadata"]["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    nb["metadata"]["language_info"] = {"name": "python", "version": "3.12"}

    cells = [
        md(
            "# Sequential Active Cooling Through the `|q,n_r,n_s>` Ladder\n\n"
            "This notebook is the interactive reproducibility companion for the storage-cooling study. "
            "It follows the same state convention used in the report, `|q,n_r,n_s>`, even though the "
            "internal simulator tensor ordering is `(transmon, storage, readout)`. The notebook is "
            "load-first by default and includes commented rerun cells for the expensive calculations."
        ),
        md(
            "## 1. Title and Overview\n\n"
            "This notebook reproduces the main claims of the study:\n"
            "1. The dressed Step A ladder `|g,0_r,n_s> <-> |f,0_r,n_s-1>` and Step B ladder "
            "`|f,0_r,n_s-1> <-> |g,1_r,n_s-1>` are resolved through `n_s = 4`.\n"
            "2. Paper-guided smooth ramps matter for Step A: the `bump` family wins for three of the four manifolds.\n"
            "3. The recommended pulses show no resolved detuning shift within the `+-2 MHz` effective-model scan.\n"
            "4. The full open-system primitive cools both basis states and non-Fock inputs."
        ),
        md(
            "## 2. Environment Setup\n\n"
            "This cell locates the study directories, loads standard analysis libraries, and defines a few "
            "helpers for reading saved artifacts and displaying figures. It also chooses the study root "
            "robustly so the notebook can be executed from either the `scripts/` directory or the study root."
        ),
        code(
            "from __future__ import annotations\n\n"
            "import csv\n"
            "import json\n"
            "from pathlib import Path\n\n"
            "import matplotlib.image as mpimg\n"
            "import matplotlib.pyplot as plt\n"
            "from IPython.display import Markdown, display\n\n"
            "def find_study_dir() -> Path:\n"
            "    candidates = [Path.cwd().resolve(), Path.cwd().resolve().parent]\n"
            "    for candidate in candidates:\n"
            "        if (candidate / 'data').exists() and (candidate / 'figures').exists() and (candidate / 'artifacts').exists():\n"
            "            return candidate\n"
            "    raise RuntimeError('Could not locate the study directory from the current working directory.')\n\n"
            "study_dir = find_study_dir()\n"
            "scripts_dir = study_dir / 'scripts'\n"
            "data_dir = study_dir / 'data'\n"
            "artifacts_dir = study_dir / 'artifacts'\n"
            "figures_dir = study_dir / 'figures'\n\n"
            "def read_json(path: Path):\n"
            "    return json.loads(path.read_text(encoding='utf-8'))\n\n"
            "def read_csv(path: Path):\n"
            "    with path.open(newline='', encoding='utf-8') as handle:\n"
            "        return list(csv.DictReader(handle))\n\n"
            "def show_png(stem: str, title: str | None = None, figsize=(6.5, 4.5)):\n"
            "    fig, ax = plt.subplots(figsize=figsize)\n"
            "    image = mpimg.imread(figures_dir / f'{stem}.png')\n"
            "    ax.imshow(image)\n"
            "    ax.axis('off')\n"
            "    if title:\n"
            "        ax.set_title(title)\n"
            "    plt.show()\n\n"
            "display(Markdown(f'Using study directory: `{study_dir}`'))"
        ),
        md(
            "## 3. User-Tunable Parameters\n\n"
            "All adjustable knobs are collected in one place here. Re-running this cell and the next one is enough "
            "to propagate a new truncation, timestep, or calibration window into the optional rerun cells."
        ),
        code(
            "params = {\n"
            "    'n_tr': 4,\n"
            "    'n_storage': 7,\n"
            "    'n_readout': 3,\n"
            "    'default_dt_ns': 0.25,\n"
            "    'ringdown_multiple': 4.0,\n"
            "    'selected_n': 4,\n"
            "    'detuning_scan_offsets_MHz': [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0],\n"
            "    'robustness_detuning_window_MHz': 1.0,\n"
            "    'robustness_amplitude_window_fraction': 0.10,\n"
            "    'coherent_alpha': 1.1,\n"
            "    'thermal_probs': [0.44, 0.28, 0.16, 0.08, 0.04],\n"
            "    'floquet_n_tr': 6,\n"
            "    'floquet_n_storage': 8,\n"
            "    'floquet_n_readout': 4,\n"
            "}\n\n"
            "param_effects = {\n"
            "    'n_tr': 'Baseline transmon truncation used by the study.',\n"
            "    'n_storage': 'Storage cutoff used for the production run.',\n"
            "    'n_readout': 'Readout cutoff used for the production run.',\n"
            "    'default_dt_ns': 'Time resolution for pulse simulation.',\n"
            "    'ringdown_multiple': 'Passive wait after Step B in units of 1/kappa_r.',\n"
            "    'selected_n': 'Which basis-state primitive to inspect in detail.',\n"
            "    'detuning_scan_offsets_MHz': 'Offsets used in the saved frequency-calibration scan.',\n"
            "    'robustness_detuning_window_MHz': 'Neighborhood used for the local sensitivity summary.',\n"
            "    'robustness_amplitude_window_fraction': 'Amplitude-error neighborhood used for the local sensitivity summary.',\n"
            "    'coherent_alpha': 'Amplitude of the coherent-state repeated-cycle test.',\n"
            "    'thermal_probs': 'Initial thermal-like storage population used in the mixed-state test.',\n"
            "    'floquet_n_tr': 'Transmon cutoff used in the qualitative Floquet check.',\n"
            "    'floquet_n_storage': 'Storage cutoff used in the qualitative Floquet check.',\n"
            "    'floquet_n_readout': 'Readout cutoff used in the qualitative Floquet check.',\n"
            "}\n\n"
            "for key, value in params.items():\n"
            "    print(f'{key:>33}: {value}')\n"
            "    print(f'    effect: {param_effects[key]}')"
        ),
        md(
            "## 4. Derived Objects\n\n"
            "This cell builds the local model, rotating frame, and noise model from the parameter cell above. "
            "It also imports the helper functions needed by the optional rerun cells."
        ),
        code(
            "import sys\n"
            "sys.path.insert(0, str(scripts_dir))\n\n"
            "from common import build_frame, build_model, build_noise\n\n"
            "model = build_model(\n"
            "    n_tr=params['n_tr'],\n"
            "    n_storage=params['n_storage'],\n"
            "    n_readout=params['n_readout'],\n"
            ")\n"
            "frame = build_frame(model)\n"
            "noise = build_noise()\n\n"
            "display(Markdown(\n"
            "    f\"Model dimensions: transmon={model.n_tr}, readout={model.n_readout}, storage={model.n_storage}.\"\n"
            "))"
        ),
        md(
            "## 5. Step-by-Step Reproduction\n\n"
            "Each subsection below follows the same pattern: a markdown explanation, a load-saved-results cell, "
            "and a commented rerun cell."
        ),
        md(
            "### 5.1 Dressed Frequencies and Labels\n\n"
            "This section reproduces the spectroscopy table that guides experiment. It shows the user-facing "
            "labels `|q,n_r,n_s>`, the Step A and Step B frequencies, and the nearest-line spacing."
        ),
        code(
            "# --- Load saved results (default) ---\n"
            "frequency_rows = read_csv(data_dir / 'frequency_table.csv')\n"
            "for row in frequency_rows:\n"
            "    print(\n"
            "        f\"n={row['n']}: {row['initial_state']} -> {row['step_a_target_state']} at {row['storage_sideband_lab_GHz']} GHz; \"\n"
            "        f\"{row['step_a_target_state']} -> {row['step_b_target_state']} at {row['readout_dump_lab_GHz']} GHz; \"\n"
            "        f\"nearest spacing = {row['nearest_storage_sideband_detuning_MHz']} MHz\"\n"
            "    )\n"
            "show_png('transition_map', title='Transition Map')"
        ),
        code(
            "# --- Re-run with current parameters ---\n"
            "# from run_study import spectrum_and_frequency_table\n"
            "# fresh_frequency_data = spectrum_and_frequency_table(model, frame)\n"
            "# fresh_frequency_data['frequency_rows']"
        ),
        md(
            "### 5.2 Pulse-Family Comparison\n\n"
            "This section reloads the best pulse for each family and each storage manifold. It is the clearest "
            "place to see the paper-guided pulse-design outcome: `bump` wins Step A for three of the four `n_s` values."
        ),
        code(
            "# --- Load saved results (default) ---\n"
            "study_results = read_json(data_dir / 'study_results.json')\n"
            "print('Best Step A pulse by storage manifold:')\n"
            "for n in range(1, 5):\n"
            "    row = study_results['best_storage'][str(n)]\n"
            "    print(f\"n={n}: {row['family']}, amp={row['amplitude_MHz']} MHz, dur={row['duration_ns']} ns, transfer={row['target_probability']:.6f}\")\n"
            "print('\\nBest Step B pulse by storage manifold:')\n"
            "for n in range(1, 5):\n"
            "    row = study_results['best_dump'][str(n)]\n"
            "    print(f\"n={n}: {row['family']}, amp={row['amplitude_MHz']} MHz, dur={row['duration_ns']} ns, transfer={row['target_probability']:.6f}\")\n"
            "show_png('storage_family_comparison', title='Step A Family Comparison')"
        ),
        code(
            "# --- Re-run with current parameters ---\n"
            "# import run_study\n"
            "# fresh_pulse_data = run_study.optimize_pulses(model, frame)\n"
            "# fresh_pulse_data['best_storage_by_family'][str(params['selected_n'])]"
        ),
        md(
            "### 5.3 Frequency Calibration and Floquet Check\n\n"
            "This section reproduces the paper-motivated strong-drive diagnostics. The calibration scan checks "
            "for an effective detuning shift around each recommended pulse, and the Floquet summary checks "
            "whether the dominant target doublet follows the expected bosonic splitting."
        ),
        code(
            "# --- Load saved results (default) ---\n"
            "calibration_rows = study_results['frequency_calibration']\n"
            "floquet_rows = study_results['floquet_summary']\n"
            "print('Calibration summary:')\n"
            "for row in calibration_rows:\n"
            "    print(f\"mode={row['mode']}, n={row['n']}, family={row['family']}, best detuning={row['optimal_detuning_MHz']} MHz\")\n"
            "print('\\nFloquet summary:')\n"
            "for row in floquet_rows:\n"
            "    print(\n"
            "        f\"mode={row['mode']}, n={row['n']}, split={row['quasienergy_split_MHz']} MHz, \"\n"
            "        f\"expected={row['expected_bosonic_split_MHz']} MHz, warning={row['warnings'][0]}\"\n"
            "    )\n"
            "show_png('frequency_calibration_curves', title='Frequency Calibration Curves')\n"
            "show_png('floquet_doublet_splitting', title='Floquet Doublet Splitting')"
        ),
        code(
            "# --- Re-run with current parameters ---\n"
            "# from run_study import calibration_frequency_scan, floquet_sideband_summary\n"
            "# fresh_calibration = calibration_frequency_scan(model, frame, study_results['best_storage'], study_results['best_dump'])\n"
            "# fresh_floquet = floquet_sideband_summary(study_results['best_storage'], study_results['best_dump'])\n"
            "# fresh_calibration['summary'], fresh_floquet"
        ),
        md(
            "### 5.4 Full Cooling Primitive\n\n"
            "This section inspects one selected basis-state primitive in detail. The final metrics tell us "
            "whether the sequence cooled the storage mode, left residual transmon excitation, or stranded "
            "population in the readout."
        ),
        code(
            "# --- Load saved results (default) ---\n"
            "selected_n = str(params['selected_n'])\n"
            "primitive = study_results['cooling_primitives'][selected_n]\n"
            "for key in ['success_probability', 'residual_same_n_probability', 'final_mean_storage_n', 'final_transmon_excited_population', 'final_readout_n']:\n"
            "    print(f'{key:>35}: {primitive[key]}')\n"
            "print('\\nDominant final basis states:')\n"
            "for label, prob in primitive['final_stage_top']:\n"
            "    print(f'  {label}: {prob:.6f}')\n"
            "show_png('cooling_primitive_trajectories', title=f'Single-Cycle Primitive for n={selected_n}')"
        ),
        code(
            "# --- Re-run with current parameters ---\n"
            "# from run_study import simulate_cooling_primitive\n"
            "# fresh_primitive = simulate_cooling_primitive(\n"
            "#     model,\n"
            "#     frame,\n"
            "#     noise,\n"
            "#     n=params['selected_n'],\n"
            "#     best_storage_case=study_results['best_storage'][str(params['selected_n'])],\n"
            "#     best_dump_case=study_results['best_dump'][str(params['selected_n'])],\n"
            "# )\n"
            "# fresh_primitive['success_probability']"
        ),
        md(
            "### 5.5 Repeated-Cycle Cooling\n\n"
            "This section reloads the repeated-cycle cooling results for a basis state, a coherent state, "
            "and a thermal-like mixture. It answers whether the protocol is genuinely useful as a cooling scheme."
        ),
        code(
            "# --- Load saved results (default) ---\n"
            "print('Basis-state ladder:')\n"
            "for row in study_results['ladder_basis_g4']:\n"
            "    print(row)\n"
            "print('\\nCoherent-state ladder:')\n"
            "for row in study_results['ladder_coherent']:\n"
            "    print(row)\n"
            "print('\\nThermal-like ladder:')\n"
            "for row in study_results['ladder_thermal']:\n"
            "    print(row)\n"
            "show_png('cooling_per_cycle', title='Repeated-Cycle Cooling')"
        ),
        code(
            "# --- Re-run with current parameters ---\n"
            "# import qutip as qt\n"
            "# from run_study import simulate_ladder_protocol\n"
            "# initial_basis = model.basis_state(0, 4, 0)\n"
            "# fresh_ladder = simulate_ladder_protocol(\n"
            "#     model,\n"
            "#     frame,\n"
            "#     noise,\n"
            "#     study_results['best_storage'],\n"
            "#     study_results['best_dump'],\n"
            "#     initial_basis,\n"
            "# )\n"
            "# fresh_ladder['cycles']"
        ),
        md(
            "## 6. Validation\n\n"
            "This section reloads the analytic resonance checks, the convergence pass, and the local sensitivity plot. "
            "It is the fastest way to confirm that the study is numerically stable and not just visually plausible."
        ),
        code(
            "# --- Load saved results (default) ---\n"
            "validation_rows = read_json(artifacts_dir / 'validation_summary.json')\n"
            "convergence_payload = read_json(artifacts_dir / 'convergence_checks.json')\n"
            "print('Analytic resonance checks:')\n"
            "for row in validation_rows:\n"
            "    print(row)\n"
            "print('\\nConvergence checks:')\n"
            "for row in convergence_payload['rows']:\n"
            "    print(row['label'], row['delta_success_probability'], row['delta_final_mean_storage_n'])\n"
            "show_png('sensitivity_heatmaps', title='Local Robustness Heatmaps', figsize=(9, 6.5))"
        ),
        code(
            "# --- Re-run with current parameters ---\n"
            "# from run_study import validation_summary\n"
            "# fresh_validation = validation_summary(model, frame, study_results['best_storage'])\n"
            "# fresh_validation"
        ),
        md(
            "## 7. Key Figures\n\n"
            "These cells re-display the main figures used in the report so the study can be skimmed quickly from the notebook."
        ),
        code(
            "show_png('transition_map', title='Transition Map')\n"
            "show_png('storage_family_comparison', title='Step A Family Comparison')\n"
            "show_png('frequency_calibration_curves', title='Frequency Calibration Curves')\n"
            "show_png('floquet_doublet_splitting', title='Floquet Doublet Splitting')\n"
            "show_png('cooling_per_cycle', title='Repeated-Cycle Cooling')\n"
            "show_png('coherent_final_wigner', title='Final Coherent-State Wigner Function', figsize=(5.5, 4.8))"
        ),
        md(
            "## 8. Summary\n\n"
            "The notebook reproduces the main conclusions of the study from the saved artifacts. The protocol is "
            "viable in the effective model through `n_s = 4`, the `bump` family is the best Step A choice for most "
            "manifolds, the Step B dump is best served by a cosine-squared pulse, and repeated application cools the "
            "storage mode strongly for both basis and non-Fock inputs."
        ),
        code(
            "summary_rows = [\n"
            "    ('step A n=1 frequency (GHz)', frequency_rows[0]['storage_sideband_lab_GHz']),\n"
            "    ('step A n=4 frequency (GHz)', frequency_rows[3]['storage_sideband_lab_GHz']),\n"
            "    ('best step A family at n=4', study_results['best_storage']['4']['family']),\n"
            "    ('single-cycle n=4 success', study_results['cooling_primitives']['4']['success_probability']),\n"
            "    ('basis ladder final mean n', study_results['ladder_basis_g4'][-1]['final_mean_storage_n']),\n"
            "    ('coherent ladder final mean n', study_results['ladder_coherent'][-1]['final_mean_storage_n']),\n"
            "    ('thermal ladder final mean n', study_results['ladder_thermal'][-1]['final_mean_storage_n']),\n"
            "]\n"
            "for label, value in summary_rows:\n"
            "    print(f'{label:>34}: {value}')"
        ),
    ]

    nb["cells"] = cells
    NOTEBOOK_PATH.write_text(nbf.writes(nb), encoding="utf-8")
    print(f"Wrote {NOTEBOOK_PATH}")


if __name__ == "__main__":
    main()
