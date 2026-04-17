"""Build the reproducibility notebook for the strict no-detuning multitone-SQR study."""

from __future__ import annotations

import json
from pathlib import Path


STUDY_DIR = Path(__file__).resolve().parent.parent
NOTEBOOK_PATH = STUDY_DIR / "scripts" / "reproducibility_notebook.ipynb"


def markdown_cell(text: str) -> dict:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": text.splitlines(keepends=True),
    }


def code_cell(text: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": text.splitlines(keepends=True),
    }


def main() -> None:
    cells = [
        markdown_cell(
            "# Reproducibility Notebook\n\n"
            "This notebook reproduces the main results of the strict no-detuning multitone-SQR study. "
            "It is designed for both users and future agents: the early parameter cell collects the main knobs, "
            "the rerun cell can regenerate the saved artifacts, and the remaining cells load the archived results "
            "to reproduce the headline tables and figures."
        ),
        code_cell(
            "from pathlib import Path\n"
            "import json\n"
            "import subprocess\n"
            "import sys\n"
            "\n"
            "# User-tunable parameters\n"
            "STUDY_DIR = Path.cwd().resolve().parent if Path.cwd().name == 'scripts' else Path(r'" + str(STUDY_DIR).replace("\\", "\\\\") + "')\n"
            "DATA_DIR = STUDY_DIR / 'data'\n"
            "FIGURES_DIR = STUDY_DIR / 'figures'\n"
            "SCRIPTS_DIR = STUDY_DIR / 'scripts'\n"
            "RERUN_MAIN = False\n"
            "RERUN_VALIDATION = False\n"
            "RUN_STUDY_ARGS = ['--n-starts', '1', '--maxiter', '5']\n"
            "SHOW_FIGURES = True\n"
        ),
        markdown_cell(
            "## Optional rerun\n\n"
            "Set `RERUN_MAIN` and/or `RERUN_VALIDATION` to `True` in the parameter cell above if you want to "
            "regenerate the archived results instead of only loading them."
        ),
        code_cell(
            "if RERUN_MAIN:\n"
            "    subprocess.run([sys.executable, str(SCRIPTS_DIR / 'run_study.py'), *RUN_STUDY_ARGS], check=True, cwd=STUDY_DIR)\n"
            "if RERUN_VALIDATION:\n"
            "    subprocess.run([sys.executable, str(SCRIPTS_DIR / 'validate_results.py')], check=True, cwd=STUDY_DIR)\n"
            "print('Rerun step complete.')\n"
        ),
        markdown_cell("## Load archived artifacts"),
        code_cell(
            "def load_json(path: Path):\n"
            "    return json.loads(path.read_text(encoding='utf-8'))\n"
            "\n"
            "results = load_json(DATA_DIR / 'study_results.json')\n"
            "summary = load_json(DATA_DIR / 'study_summary.json')\n"
            "analytic = load_json(DATA_DIR / 'analytic_summary.json')\n"
            "validation = load_json(DATA_DIR / 'validation_summary.json')\n"
            "validation_details = load_json(DATA_DIR / 'validation_details.json')\n"
            "audit = load_json(DATA_DIR / 'prior_audit.json')\n"
            "print('Loaded rows:', len(results['case_rows']))\n"
            "print('Runtime [s]:', results['runtime_s'])\n"
        ),
        markdown_cell("## Headline findings"),
        code_cell(
            "headline = {\n"
            "    'strict_shared_line_mean_fidelity': summary['strict_shared_line_mean_fidelity'],\n"
            "    'strict_shared_line_best_fidelity': summary['strict_shared_line_best_fidelity'],\n"
            "    'strict_shared_line_mean_max_residual_z_error_rad': summary['strict_shared_line_mean_max_residual_z_error_rad'],\n"
            "    'decoupled_block_min_fidelity': summary['decoupled_block_min_fidelity'],\n"
            "    'ideal_echo_mean_fidelity': summary['ideal_echo_mean_fidelity'],\n"
            "    'finite_echo_mean_fidelity': summary['finite_echo_mean_fidelity'],\n"
            "}\n"
            "headline\n"
        ),
        markdown_cell("## Tabulate the saved case rows"),
        code_cell(
            "import pandas as pd\n"
            "\n"
            "df = pd.DataFrame(results['case_rows'])\n"
            "full = df[df['construction'] == 'full_shared_line']\n"
            "display(full.groupby('family')['restricted_average_gate_fidelity'].agg(['mean', 'min', 'max']).round(6))\n"
            "display(full.groupby('n_active')['restricted_average_gate_fidelity'].agg(['mean', 'min', 'max']).round(6))\n"
            "display(df.groupby('construction')['restricted_average_gate_fidelity'].agg(['mean', 'min', 'max']).round(6))\n"
        ),
        markdown_cell("## Validation details"),
        code_cell(
            "validation_details\n"
        ),
        markdown_cell("## Prior-work audit"),
        code_cell(
            "audit\n"
        ),
        markdown_cell("## Main figures"),
        code_cell(
            "if SHOW_FIGURES:\n"
            "    from IPython.display import Image, display\n"
            "    for name in [\n"
            "        'duration_fidelity_tradeoff.png',\n"
            "        'blockwise_residual_z_vs_duration.png',\n"
            "        'addressed_subspace_scaling.png',\n"
            "        'plain_vs_echo_comparison.png',\n"
            "    ]:\n"
            "        path = FIGURES_DIR / name\n"
            "        print(path.name)\n"
            "        display(Image(filename=str(path)))\n"
        ),
        markdown_cell(
            "## Next steps\n\n"
            "If you want to extend the study, the most natural follow-ups are:\n"
            "1. Replace the strict shared-line square multitone ansatz with a richer segmented or sampled envelope while still forbidding artificial per-tone detuning.\n"
            "2. Test whether selective or manifold-adaptive refocusing pulses can outperform the strict echoed construction studied here.\n"
            "3. Extend the same no-go analysis to a model with explicit open-system noise or additional hardware filtering."
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
                "version": "3.12.10",
            },
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    NOTEBOOK_PATH.write_text(json.dumps(notebook, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
