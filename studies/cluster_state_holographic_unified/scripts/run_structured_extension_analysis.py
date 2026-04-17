from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path
from typing import Any, Sequence

import matplotlib.pyplot as plt
import numpy as np
import runtime_compat  # noqa: F401
import qutip as qt


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common as c
import run_design_space_study as ds


STYLE_PATH = SCRIPT_DIR.parents[2] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))


DEFAULT_N_CAV = 12
DEFAULT_LOGICAL_MAXITER = 8
DEFAULT_AUGMENTED_MAXITER = 24
DEFAULT_ORDERING_MAXITER = 24
WIGNER_XVEC = np.linspace(-3.0, 3.0, 151)

BEST_ARTIFACTS = {
    "drsqr": ("corrected_best_sqr.json", "Best D + R + SQR"),
    "drcpsqr": ("corrected_best_cpsqr.json", "Best D + R + CPSQR"),
}


def _write_summary_checkpoint(
    *,
    n_cav: int,
    logical_maxiter: int,
    augmented_maxiter: int,
    ordering_maxiter: int,
    parameter_tables: dict[str, Any],
    basis_extension: dict[str, Any],
    ordering_study: dict[str, Any],
) -> None:
    summary = {
        "study_name": "cluster_state_holographic_unified",
        "date_created": time.strftime("%Y-%m-%d"),
        "n_cav": int(n_cav),
        "logical_maxiter": int(logical_maxiter),
        "augmented_maxiter": int(augmented_maxiter),
        "ordering_maxiter": int(ordering_maxiter),
        "parameter_tables": parameter_tables,
        "basis_extension": {
            family_key: _slim_basis_payload(payload, n_cav=int(n_cav))
            for family_key, payload in basis_extension.items()
        },
        "ordering_study": ordering_study,
    }
    c.save_json(c.DATA_DIR / "structured_extension_summary.json", summary)


def _load_best_record(family_key: str) -> dict[str, Any]:
    artifact_name, _label = BEST_ARTIFACTS[str(family_key)]
    payload = c.load_json(c.ARTIFACT_DIR / artifact_name)
    return dict(payload["record"])


def _sequence_from_record(record: dict[str, Any], *, n_cav: int) -> c.GateSequence:
    return ds.apply_solution_to_case(record, n_cav=int(n_cav))


def _warm_start_from_record(record: dict[str, Any], *, n_cav: int) -> dict[str, Any]:
    sequence = _sequence_from_record(record, n_cav=int(n_cav))
    return {
        "parameter_vector": sequence.get_parameter_vector().tolist(),
        "time_raw_vector": sequence.get_time_raw_vector(active_only=True).tolist(),
    }


def _rebuild_case(base_record: dict[str, Any], *, order_tokens: Sequence[str] | None = None) -> dict[str, Any]:
    case = ds.case_from_record(base_record)
    if order_tokens is not None:
        case["order_tokens"] = tuple(str(token) for token in order_tokens)
        case["order_label"] = ds.compact_order_label(order_tokens)
        builder_kwargs = dict(case.get("builder_kwargs", {}))
        builder_kwargs["order"] = tuple(str(token) for token in order_tokens)
        case["builder_kwargs"] = builder_kwargs
        case["variant_key"] = f"blk{int(case['blocks'])}_a{int(case['max_tones'])}_ord{case['order_label']}"
        case["variant_label"] = f"{int(case['blocks'])} blocks / {case['order_label']} / n_active={int(case['max_tones'])}"
    return case


def _record_from_fit(
    base_record: dict[str, Any],
    fit: dict[str, Any],
    *,
    n_cav: int,
    search_phase: str,
    seed: int,
    subspace_label: str,
    order_tokens: Sequence[str] | None = None,
) -> dict[str, Any]:
    case = _rebuild_case(base_record, order_tokens=order_tokens)
    sequence = fit["result"].sequence
    record = {
        "case_id": f"{ds.case_id(case)}_{subspace_label}",
        "family_key": str(case["family_key"]),
        "family_label": str(case["family_label"]),
        "variant_key": str(case["variant_key"]),
        "variant_label": f"{case['variant_label']} / {subspace_label}",
        "builder_name": str(case["builder_name"]),
        "builder_kwargs": dict(case.get("builder_kwargs", {})),
        "order_tokens": [str(token) for token in case.get("order_tokens", ())],
        "order_label": str(case.get("order_label", "")),
        "levels": None if case.get("levels") is None else [int(level) for level in case["levels"]],
        "max_tones": None if case.get("max_tones") is None else int(case["max_tones"]),
        "blocks": None if case.get("blocks") is None else int(case["blocks"]),
        "search_phase": str(search_phase),
        "optimization_n_cav": int(n_cav),
        "seed": int(seed),
        "init_guess": "heuristic",
        "maxiter": int(fit.get("maxiter", 0)) if isinstance(fit.get("maxiter"), int) else None,
        "multistart": 1,
        "fidelity": float(fit["fidelity"]),
        "objective": float(fit["objective"]),
        "success": bool(fit["success"]),
        "message": str(fit["message"]),
        "summary": dict(fit["summary"]),
        "metrics": dict(fit["metrics"]),
        "sequence": fit["sequence_payload"],
        "parameter_vector": sequence.get_parameter_vector().tolist(),
        "time_vector": sequence.get_time_vector(active_only=False).tolist(),
    }
    record["physical"] = ds.evaluate_physical_candidate(record)
    return record


def _augmented_spectator_levels(record: dict[str, Any], *, n_cav: int) -> list[int]:
    active = record["physical"]["by_n_cav"][str(int(n_cav))]["active_subspace"]["candidate_active_levels"]
    return [int(level) for level in active if int(level) not in (0, 1)]


def _augmented_target_for_spectators(spectator_levels: Sequence[int], *, n_cav: int) -> tuple[c.Subspace, np.ndarray]:
    levels = [0, 1] + [int(level) for level in spectator_levels]
    subspace = c.level_subspace(levels, n_cav=int(n_cav))
    spectator_dim = 2 * len(tuple(int(level) for level in spectator_levels))
    target = np.eye(4 + spectator_dim, dtype=np.complex128)
    target[:4, :4] = c.TARGET_UNITARY
    return subspace, target


def _wigner_metrics_for_record(record: dict[str, Any], *, n_cav: int) -> dict[str, Any]:
    by_input: dict[str, Any] = {}
    rms_values: list[float] = []
    l1_values: list[float] = []
    state_fidelities: list[float] = []
    for label in c.LOGICAL_LABELS:
        target_state = ds.target_output_state(label, n_cav=int(n_cav))
        candidate_state = ds.candidate_output_state(record, n_cav=int(n_cav), input_label=label)
        target_rho = ds.reduced_cavity_density(target_state, n_cav=int(n_cav))
        candidate_rho = ds.reduced_cavity_density(candidate_state, n_cav=int(n_cav))
        target_wigner = qt.wigner(target_rho, WIGNER_XVEC, WIGNER_XVEC)
        candidate_wigner = qt.wigner(candidate_rho, WIGNER_XVEC, WIGNER_XVEC)
        delta = candidate_wigner - target_wigner
        rms = float(np.sqrt(np.mean(np.square(delta))))
        l1 = float(np.mean(np.abs(delta)))
        state_fidelity = float(qt.fidelity(target_rho, candidate_rho))
        by_input[str(label)] = {
            "reduced_state_fidelity": state_fidelity,
            "wigner_rms": rms,
            "wigner_l1": l1,
        }
        rms_values.append(rms)
        l1_values.append(l1)
        state_fidelities.append(state_fidelity)
    return {
        "by_input": by_input,
        "mean_reduced_state_fidelity": float(np.mean(state_fidelities)),
        "min_reduced_state_fidelity": float(np.min(state_fidelities)),
        "mean_wigner_rms": float(np.mean(rms_values)),
        "max_wigner_rms": float(np.max(rms_values)),
        "mean_wigner_l1": float(np.mean(l1_values)),
    }


def _plot_wigner_triptych(
    *,
    title: str,
    baseline_record: dict[str, Any],
    augmented_record: dict[str, Any],
    stem: str,
    n_cav: int,
) -> None:
    selected_rows: list[tuple[str, dict[str, Any] | None]] = [
        ("Target", None),
        ("Logical-basis refit", baseline_record),
        ("Augmented-basis refit", augmented_record),
    ]

    fig = plt.figure(figsize=(8.8, 9.8), constrained_layout=True)
    grid = fig.add_gridspec(
        len(c.LOGICAL_LABELS),
        len(selected_rows) + 1,
        width_ratios=[1.0] * len(selected_rows) + [0.08],
        wspace=0.08,
        hspace=0.10,
    )
    axes = np.empty((len(c.LOGICAL_LABELS), len(selected_rows)), dtype=object)
    vmax = 0.35
    mesh = None
    shared_axis = None
    for row_index, input_label in enumerate(c.LOGICAL_LABELS):
        for col_index, (panel_title, record) in enumerate(selected_rows):
            if shared_axis is None:
                ax = fig.add_subplot(grid[row_index, col_index])
                shared_axis = ax
            else:
                ax = fig.add_subplot(grid[row_index, col_index], sharex=shared_axis, sharey=shared_axis)
            axes[row_index, col_index] = ax
            if record is None:
                state = ds.target_output_state(input_label, n_cav=int(n_cav))
            else:
                state = ds.candidate_output_state(record, n_cav=int(n_cav), input_label=input_label)
            rho_cav = ds.reduced_cavity_density(state, n_cav=int(n_cav))
            wig = qt.wigner(rho_cav, WIGNER_XVEC, WIGNER_XVEC)
            mesh = ax.pcolormesh(WIGNER_XVEC, WIGNER_XVEC, wig, cmap="RdBu_r", shading="auto", vmin=-vmax, vmax=vmax, rasterized=True)
            ax.set_aspect("equal")
            if row_index == 0:
                ax.set_title(panel_title)
            if col_index == 0:
                ax.set_ylabel(f"{input_label}\nIm(alpha)")
            if row_index == len(c.LOGICAL_LABELS) - 1:
                ax.set_xlabel("Re(alpha)")
    for ax in axes.ravel():
        ax.label_outer()
    cax = fig.add_subplot(grid[:, -1])
    colorbar = fig.colorbar(mesh, cax=cax)
    colorbar.set_label("W(alpha)")
    fig.suptitle(title)
    fig.savefig(c.FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(c.FIG_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _format_vector(values: Sequence[float], *, precision: int = 3) -> str:
    fmt = f"{{:.{int(precision)}f}}"
    return "[" + ", ".join(fmt.format(float(value)) for value in values) + "]"


def _gate_parameter_rows(record: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for gate in record["sequence"]:
        params = [float(value) for value in gate.get("parameters", [])]
        metadata = gate.get("metadata", {}) if isinstance(gate.get("metadata"), dict) else {}
        note = ""
        if gate["type"] == "Displacement":
            parameter_text = f"Re(alpha)={params[0]:.3f}, Im(alpha)={params[1]:.3f}"
        elif gate["type"] == "QubitRotation":
            parameter_text = f"theta={params[0]:.3f}, phi={params[1]:.3f}"
        elif gate["type"] == "PrimitiveGate" and metadata.get("ideal_kind") == "MaskedSQR":
            tone_count = len(metadata.get("levels", []))
            theta = params[:tone_count]
            phi = params[tone_count:]
            parameter_text = f"theta={_format_vector(theta)}, phi={_format_vector(phi)}"
            note = "levels=" + ",".join(str(level) for level in metadata.get("levels", []))
        elif gate["type"] == "PrimitiveGate" and metadata.get("ideal_kind") == "MaskedCPSQR":
            parameter_text = f"phases={_format_vector(params)}"
            note = "levels=" + ",".join(str(level) for level in metadata.get("levels", []))
        else:
            parameter_text = _format_vector(params) if params else "-"
        rows.append(
            {
                "gate": str(gate["name"]),
                "kind": str(gate["type"] if gate["type"] != "PrimitiveGate" else metadata.get("ideal_kind", gate["type"])),
                "duration_ns": float(gate["duration"]) * 1.0e9,
                "parameters": parameter_text,
                "note": note,
            }
        )
    return rows


def _write_parameter_table_files(family_key: str, record: dict[str, Any]) -> dict[str, str]:
    rows = _gate_parameter_rows(record)
    csv_path = c.DATA_DIR / f"{family_key}_parameter_table.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gate", "kind", "duration_ns", "parameters", "note"])
        writer.writeheader()
        writer.writerows(rows)

    generated_dir = c.REPORT_DIR / "generated"
    generated_dir.mkdir(parents=True, exist_ok=True)
    tex_path = generated_dir / f"{family_key}_parameter_rows.tex"
    tex_lines = []
    for row in rows:
        tex_lines.append(
            f"{row['gate']} & {row['kind']} & {row['duration_ns']:.1f} & {row['parameters']} & {row['note'] or '-'} \\\\" 
        )
    tex_path.write_text("\n".join(tex_lines) + "\n", encoding="utf-8")
    return {"csv": str(csv_path), "tex_rows": str(tex_path)}


def _run_basis_extension(
    base_record: dict[str, Any],
    *,
    family_key: str,
    n_cav: int,
    logical_maxiter: int,
    augmented_maxiter: int,
    use_fast_path: bool,
) -> dict[str, Any]:
    print(f"[structured-extension] {family_key}: logical refit start", flush=True)
    logical_warm_start = _warm_start_from_record(base_record, n_cav=int(n_cav))
    logical_sequence = _sequence_from_record(base_record, n_cav=int(n_cav))
    logical_fit = c.fit_sequence(
        logical_sequence,
        n_cav=int(n_cav),
        seed=17,
        init_guess="heuristic",
        multistart=1,
        maxiter=int(logical_maxiter),
        use_fast_path=bool(use_fast_path),
        warm_start=logical_warm_start,
    )
    logical_fit["maxiter"] = int(logical_maxiter)
    logical_record = _record_from_fit(
        base_record,
        logical_fit,
        n_cav=int(n_cav),
        search_phase="logical_refit",
        seed=17,
        subspace_label="logical_refit",
    )
    logical_wigner = _wigner_metrics_for_record(logical_record, n_cav=int(n_cav))
    print(
        f"[structured-extension] {family_key}: logical refit complete "
        f"F12={logical_record['physical']['by_n_cav'][str(int(n_cav))]['fidelity']:.6f}",
        flush=True,
    )

    spectator_levels = _augmented_spectator_levels(logical_record, n_cav=int(n_cav))
    print(
        f"[structured-extension] {family_key}: augmented refit start spectators={spectator_levels}",
        flush=True,
    )
    augmented_subspace, augmented_target = _augmented_target_for_spectators(spectator_levels, n_cav=int(n_cav))
    augmented_sequence = _sequence_from_record(base_record, n_cav=int(n_cav))
    augmented_fit = c.fit_sequence(
        augmented_sequence,
        n_cav=int(n_cav),
        seed=17,
        init_guess="heuristic",
        multistart=1,
        maxiter=int(augmented_maxiter),
        use_fast_path=bool(use_fast_path),
        warm_start=logical_warm_start,
        target_unitary=augmented_target,
        subspace=augmented_subspace,
    )
    augmented_fit["maxiter"] = int(augmented_maxiter)
    augmented_record = _record_from_fit(
        base_record,
        augmented_fit,
        n_cav=int(n_cav),
        search_phase="augmented_refit",
        seed=17,
        subspace_label="augmented_refit",
    )
    augmented_wigner = _wigner_metrics_for_record(augmented_record, n_cav=int(n_cav))
    augmented_target_eval = c.evaluate_sequence(
        _sequence_from_record(augmented_record, n_cav=int(n_cav)),
        n_cav=int(n_cav),
        target_unitary=augmented_target,
        subspace=augmented_subspace,
    )
    print(
        f"[structured-extension] {family_key}: augmented refit complete "
        f"F12={augmented_record['physical']['by_n_cav'][str(int(n_cav))]['fidelity']:.6f}",
        flush=True,
    )

    interpretation = (
        "supports"
        if augmented_wigner["mean_wigner_rms"] < logical_wigner["mean_wigner_rms"]
        and augmented_record["physical"]["by_n_cav"][str(int(n_cav))]["fidelity"] >= logical_record["physical"]["by_n_cav"][str(int(n_cav))]["fidelity"]
        else "does_not_support"
    )
    stem = f"{family_key}_wigner_basis_extension"
    title = f"{BEST_ARTIFACTS[family_key][1]}: logical vs augmented spectator-constrained refit"
    print(f"[structured-extension] {family_key}: plotting Wigner triptych", flush=True)
    _plot_wigner_triptych(
        title=title,
        baseline_record=logical_record,
        augmented_record=augmented_record,
        stem=stem,
        n_cav=int(n_cav),
    )
    return {
        "base_case_id": str(base_record["case_id"]),
        "spectator_levels": [int(level) for level in spectator_levels],
        "logical_refit": {
            "record": logical_record,
            "wigner": logical_wigner,
        },
        "augmented_refit": {
            "record": augmented_record,
            "wigner": augmented_wigner,
            "augmented_target_fidelity": float(augmented_target_eval["fidelity"]),
        },
        "interpretation": interpretation,
        "figure_stem": stem,
    }


def _slim_basis_payload(payload: dict[str, Any], *, n_cav: int) -> dict[str, Any]:
    logical_record = payload["logical_refit"]["record"]
    augmented_record = payload["augmented_refit"]["record"]
    logical_n = logical_record["physical"]["by_n_cav"][str(int(n_cav))]
    augmented_n = augmented_record["physical"]["by_n_cav"][str(int(n_cav))]
    return {
        "base_case_id": payload["base_case_id"],
        "spectator_levels": payload["spectator_levels"],
        "interpretation": payload["interpretation"],
        "figure_stem": payload["figure_stem"],
        "logical_refit": {
            "case_id": logical_record["case_id"],
            "optimized_fidelity": float(logical_record["fidelity"]),
            "physical_fidelity_n12": float(logical_n["fidelity"]),
            "leakage_worst_n12": float(logical_n["leakage_worst"]),
            "wigner": payload["logical_refit"]["wigner"],
        },
        "augmented_refit": {
            "case_id": augmented_record["case_id"],
            "optimized_fidelity": float(augmented_record["fidelity"]),
            "physical_fidelity_n12": float(augmented_n["fidelity"]),
            "leakage_worst_n12": float(augmented_n["leakage_worst"]),
            "augmented_target_fidelity": float(payload["augmented_refit"]["augmented_target_fidelity"]),
            "wigner": payload["augmented_refit"]["wigner"],
        },
    }


def _run_ordering_study(
    base_record: dict[str, Any],
    *,
    n_cav: int,
    maxiter: int,
    use_fast_path: bool,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for order_tokens in ds.SQR_ORDERINGS:
        label = ds.compact_order_label(order_tokens)
        case = _rebuild_case(base_record, order_tokens=order_tokens)
        sequence = ds.build_sequence_from_case(case, n_cav=int(n_cav))
        started = time.perf_counter()
        fit = c.fit_sequence(
            sequence,
            n_cav=int(n_cav),
            seed=17,
            init_guess="heuristic",
            multistart=1,
            maxiter=int(maxiter),
            use_fast_path=bool(use_fast_path),
        )
        fit["maxiter"] = int(maxiter)
        record = _record_from_fit(
            base_record,
            fit,
            n_cav=int(n_cav),
            search_phase="ordering_refine",
            seed=17,
            subspace_label=f"ordering_{label}",
            order_tokens=order_tokens,
        )
        n12 = record["physical"]["by_n_cav"][str(int(n_cav))]
        rows.append(
            {
                "order_label": label,
                "order_pretty": ds.pretty_order_label(order_tokens),
                "optimized_fidelity": float(record["fidelity"]),
                "physical_fidelity_n12": float(n12["fidelity"]),
                "leakage_worst_n12": float(n12["leakage_worst"]),
                "total_duration_ns": float(record["summary"]["total_duration_ns"]),
                "elapsed_s": float(time.perf_counter() - started),
                "record": record,
            }
        )
    ordered = sorted(rows, key=lambda row: row["physical_fidelity_n12"], reverse=True)
    csv_path = c.DATA_DIR / "sqr_ordering_fixed_budget.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["order_label", "order_pretty", "optimized_fidelity", "physical_fidelity_n12", "leakage_worst_n12", "total_duration_ns", "elapsed_s"],
        )
        writer.writeheader()
        for row in ordered:
            writer.writerow({key: row[key] for key in writer.fieldnames})

    labels = [row["order_label"] for row in ordered]
    fidelities = [row["physical_fidelity_n12"] for row in ordered]
    leakages = [row["leakage_worst_n12"] for row in ordered]
    fig, ax1 = plt.subplots(figsize=(7.2, 4.2))
    ax1.bar(labels, fidelities, color="#1b6ca8", alpha=0.9)
    ax1.set_ylabel(r"$F_{12}$")
    ax1.set_ylim(min(fidelities) - 0.01, max(fidelities) + 0.01)
    ax1.set_xlabel("Ordering")
    ax2 = ax1.twinx()
    ax2.plot(labels, leakages, color="#c84c09", marker="o", linewidth=1.5)
    ax2.set_ylabel("Worst leakage")
    fig.tight_layout()
    fig.savefig(c.FIG_DIR / "sqr_ordering_fixed_budget.pdf", bbox_inches="tight")
    fig.savefig(c.FIG_DIR / "sqr_ordering_fixed_budget.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    spread = float(max(fidelities) - min(fidelities)) if fidelities else float("nan")
    return {
        "rows": [
            {key: value for key, value in row.items() if key != "record"}
            for row in ordered
        ],
        "best_order": ordered[0]["order_label"],
        "worst_order": ordered[-1]["order_label"],
        "physical_fidelity_spread": spread,
        "csv_path": str(csv_path),
        "figure_stem": "sqr_ordering_fixed_budget",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extend the unified structured cluster-state study.")
    parser.add_argument("--n-cav", type=int, default=DEFAULT_N_CAV)
    parser.add_argument("--logical-maxiter", type=int, default=DEFAULT_LOGICAL_MAXITER)
    parser.add_argument("--augmented-maxiter", type=int, default=DEFAULT_AUGMENTED_MAXITER)
    parser.add_argument("--ordering-maxiter", type=int, default=DEFAULT_ORDERING_MAXITER)
    parser.add_argument("--skip-basis", action="store_true")
    parser.add_argument("--skip-ordering", action="store_true")
    parser.add_argument("--disable-fast-path", action="store_true")
    args = parser.parse_args()

    best_records = {family_key: _load_best_record(family_key) for family_key in BEST_ARTIFACTS}
    parameter_tables = {
        family_key: _write_parameter_table_files(family_key, record)
        for family_key, record in best_records.items()
    }
    print("[structured-extension] parameter tables exported", flush=True)

    basis_extension = {}
    ordering_study: dict[str, Any] = {}
    _write_summary_checkpoint(
        n_cav=int(args.n_cav),
        logical_maxiter=int(args.logical_maxiter),
        augmented_maxiter=int(args.augmented_maxiter),
        ordering_maxiter=int(args.ordering_maxiter),
        parameter_tables=parameter_tables,
        basis_extension=basis_extension,
        ordering_study=ordering_study,
    )
    print("[structured-extension] wrote initial checkpoint", flush=True)

    if not args.skip_basis:
        for family_key, record in best_records.items():
            print(f"[structured-extension] basis extension start for {family_key}", flush=True)
            payload = _run_basis_extension(
                record,
                family_key=family_key,
                n_cav=int(args.n_cav),
                logical_maxiter=int(args.logical_maxiter),
                augmented_maxiter=int(args.augmented_maxiter),
                use_fast_path=not bool(args.disable_fast_path),
            )
            basis_extension[family_key] = payload
            c.save_json(c.DATA_DIR / f"{family_key}_basis_extension_detail.json", payload)
            _write_summary_checkpoint(
                n_cav=int(args.n_cav),
                logical_maxiter=int(args.logical_maxiter),
                augmented_maxiter=int(args.augmented_maxiter),
                ordering_maxiter=int(args.ordering_maxiter),
                parameter_tables=parameter_tables,
                basis_extension=basis_extension,
                ordering_study=ordering_study,
            )
            print(f"[structured-extension] basis extension complete for {family_key}", flush=True)

    if not args.skip_ordering:
        print("[structured-extension] ordering study start", flush=True)
        ordering_study = _run_ordering_study(
            best_records["drsqr"],
            n_cav=int(args.n_cav),
            maxiter=int(args.ordering_maxiter),
            use_fast_path=not bool(args.disable_fast_path),
        )
        _write_summary_checkpoint(
            n_cav=int(args.n_cav),
            logical_maxiter=int(args.logical_maxiter),
            augmented_maxiter=int(args.augmented_maxiter),
            ordering_maxiter=int(args.ordering_maxiter),
            parameter_tables=parameter_tables,
            basis_extension=basis_extension,
            ordering_study=ordering_study,
        )

    _write_summary_checkpoint(
        n_cav=int(args.n_cav),
        logical_maxiter=int(args.logical_maxiter),
        augmented_maxiter=int(args.augmented_maxiter),
        ordering_maxiter=int(args.ordering_maxiter),
        parameter_tables=parameter_tables,
        basis_extension=basis_extension,
        ordering_study=ordering_study,
    )
    print("[structured-extension] wrote data/structured_extension_summary.json", flush=True)
    for family_key, payload in basis_extension.items():
        logical = payload["logical_refit"]
        augmented = payload["augmented_refit"]
        logical_f12 = logical["record"]["physical"]["by_n_cav"][str(int(args.n_cav))]["fidelity"]
        augmented_f12 = augmented["record"]["physical"]["by_n_cav"][str(int(args.n_cav))]["fidelity"]
        print(
            f"[structured-extension] {family_key}: F12 logical={logical_f12:.6f}, F12 augmented={augmented_f12:.6f}, "
            f"Wigner RMS logical={logical['wigner']['mean_wigner_rms']:.5f}, augmented={augmented['wigner']['mean_wigner_rms']:.5f}",
            flush=True,
        )
    if ordering_study:
        print(
            f"[structured-extension] best fixed-budget ordering={ordering_study['best_order']} "
            f"spread={ordering_study['physical_fidelity_spread']:.6f}",
            flush=True,
        )


if __name__ == "__main__":
    main()