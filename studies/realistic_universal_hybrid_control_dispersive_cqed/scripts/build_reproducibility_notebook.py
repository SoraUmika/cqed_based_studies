from __future__ import annotations

from pathlib import Path

import nbformat as nbf

from common import DATA_DIR, FIGURES_DIR, STUDY_ROOT


def main() -> None:
    nb = nbf.v4.new_notebook()
    cells = []

    cells.append(
        nbf.v4.new_markdown_cell(
            "# Realistic Universal Hybrid Control in Dispersive cQED\n\n"
            "This notebook reproduces the synthesis study that asks whether the ideal primitive "
            "universal-control claim for a transmon-cavity system survives the realistic dispersive "
            "Hamiltonian once `chi`, `chi'`, Kerr, and finite pulse durations are enforced."
        )
    )
    cells.append(
        nbf.v4.new_markdown_cell(
            "## Environment Setup\n\n"
            "The notebook only loads saved results by default. It can optionally rerun the synthesis "
            "builder script, which itself aggregates validated repository artifacts rather than rerunning "
            "every upstream study."
        )
    )
    cells.append(
        nbf.v4.new_code_cell(
            "from pathlib import Path\n"
            "import json\n"
            "import subprocess\n"
            "import matplotlib.pyplot as plt\n"
            "from matplotlib.image import imread\n\n"
            f"study_root = Path(r'{STUDY_ROOT}')\n"
            f"data_dir = Path(r'{DATA_DIR}')\n"
            f"figures_dir = Path(r'{FIGURES_DIR}')\n"
            "print('Study root:', study_root)\n"
            "print('Data dir:', data_dir)\n"
            "print('Figures dir:', figures_dir)"
        )
    )

    cells.append(
        nbf.v4.new_markdown_cell(
            "## User-Tunable Parameters\n\n"
            "All user-facing knobs are collected here. The default workflow loads saved artifacts, but "
            "you can flip `rerun_builder` to regenerate the synthesis summary and figures."
        )
    )
    cells.append(
        nbf.v4.new_code_cell(
            "rerun_builder = False\n"
            "show_phase_budget = True\n"
            "show_timescale_hierarchy = True\n"
            "max_rows_to_print = 7\n\n"
            "print('rerun_builder =', rerun_builder)\n"
            "print('show_phase_budget =', show_phase_budget)\n"
            "print('show_timescale_hierarchy =', show_timescale_hierarchy)\n"
            "print('max_rows_to_print =', max_rows_to_print)"
        )
    )

    cells.append(
        nbf.v4.new_markdown_cell(
            "## Derived Objects\n\n"
            "This cell derives the saved-file paths from the tunable parameters so that rerunning "
            "Sections 3 and 4 is enough to propagate any parameter change."
        )
    )
    cells.append(
        nbf.v4.new_code_cell(
            "summary_path = data_dir / 'synthesis_summary.json'\n"
            "timescale_png = figures_dir / 'timescale_hierarchy.png'\n"
            "phase_budget_png = figures_dir / 'phase_budget.png'\n"
            "builder_script = study_root / 'scripts' / 'build_synthesis_dataset.py'\n"
            "summary_path, timescale_png, phase_budget_png, builder_script"
        )
    )

    cells.append(
        nbf.v4.new_markdown_cell(
            "## Step-by-Step Reproduction\n\n"
            "The default path is to load the saved synthesis outputs. A commented rerun cell is included "
            "immediately after it."
        )
    )
    cells.append(
        nbf.v4.new_code_cell(
            "# --- Load saved results (default) ---\n"
            "summary = json.loads(summary_path.read_text(encoding='utf-8'))\n"
            "summary['top_level_verdict']"
        )
    )
    cells.append(
        nbf.v4.new_code_cell(
            "# --- Re-run with current parameters ---\n"
            "# if rerun_builder:\n"
            "#     subprocess.run(['python', str(builder_script)], check=True)\n"
            "#     summary = json.loads(summary_path.read_text(encoding='utf-8'))\n"
            "# summary['top_level_verdict']"
        )
    )

    cells.append(
        nbf.v4.new_markdown_cell(
            "## Primitive Verdict Table\n\n"
            "The table below reproduces the unified primitive-level verdicts used in the report."
        )
    )
    cells.append(
        nbf.v4.new_code_cell(
            "rows = summary['primitive_verdicts'][:max_rows_to_print]\n"
            "for row in rows:\n"
            "    print(row['primitive'])\n"
            "    print('  family   :', row['best_structured_family'])\n"
            "    print('  metric   :', row['best_metric_label'], row['best_metric_value'])\n"
            "    print('  chiT/2pi :', row['chi_t_over_2pi'])\n"
            "    print('  verdict  :', row['verdict'])\n"
            "    print('  obstacle :', row['main_obstruction'])\n"
            "    print()"
        )
    )

    cells.append(
        nbf.v4.new_markdown_cell(
            "## Validation\n\n"
            "This section reproduces the study's key consistency checks: the derived timing hierarchy and "
            "the phase-budget trends."
        )
    )
    cells.append(
        nbf.v4.new_code_cell(
            "summary['derived_timescales']"
        )
    )
    cells.append(
        nbf.v4.new_code_cell(
            "phase_rows = summary['phase_budget_rows'][:12]\n"
            "phase_rows"
        )
    )

    cells.append(
        nbf.v4.new_markdown_cell(
            "## Key Figures\n\n"
            "These are the report figures generated by the synthesis script."
        )
    )
    cells.append(
        nbf.v4.new_code_cell(
            "if show_timescale_hierarchy:\n"
            "    img = imread(timescale_png)\n"
            "    plt.figure(figsize=(12, 6))\n"
            "    plt.imshow(img)\n"
            "    plt.axis('off')\n"
            "    plt.show()"
        )
    )
    cells.append(
        nbf.v4.new_code_cell(
            "if show_phase_budget:\n"
            "    img = imread(phase_budget_png)\n"
            "    plt.figure(figsize=(12, 6))\n"
            "    plt.imshow(img)\n"
            "    plt.axis('off')\n"
            "    plt.show()"
        )
    )

    cells.append(
        nbf.v4.new_markdown_cell(
            "## Summary\n\n"
            "The saved outputs support a clear conclusion: the strict ideal primitive gate set does not survive "
            "the realistic dispersive Hamiltonian as a literal statement. A weaker phase-aware constructive "
            "library does survive, but it is low-dimensional, gauge-relaxed, and not yet a fully demonstrated "
            "non-GRAPE universal stack.\n\n"
            "| Tunable parameter | Default value | Effect |\n"
            "|---|---:|---|\n"
            "| `rerun_builder` | `False` | Regenerates the unified summary and figures from the aggregation script. |\n"
            "| `show_phase_budget` | `True` | Displays the phase-budget figure. |\n"
            "| `show_timescale_hierarchy` | `True` | Displays the timing-hierarchy figure. |\n"
            "| `max_rows_to_print` | `7` | Limits how many primitive verdict rows are printed in the notebook. |"
        )
    )

    nb["cells"] = cells
    nb["metadata"] = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.12"},
    }

    out_path = STUDY_ROOT / "scripts" / "reproducibility_notebook.ipynb"
    out_path.write_text(nbf.writes(nb), encoding="utf-8")


if __name__ == "__main__":
    main()
