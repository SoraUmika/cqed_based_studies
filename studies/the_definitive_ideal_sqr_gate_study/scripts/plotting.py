from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DIRECT_COLOR = "#0072B2"
SYMMETRIC_ECHO_COLOR = "#E69F00"
INDEPENDENT_ECHO_COLOR = "#009E73"
ASYMMETRIC_ECHO_COLOR = "#D55E00"
HYBRID_COLOR = "#CC79A7"
OTHER_COLOR = "#666666"


def construction_color(construction_family: str) -> str:
    family = str(construction_family).lower()
    if family == "direct":
        return DIRECT_COLOR
    if family == "symmetric echo":
        return SYMMETRIC_ECHO_COLOR
    if family == "independent echo":
        return INDEPENDENT_ECHO_COLOR
    if family == "asymmetric echo":
        return ASYMMETRIC_ECHO_COLOR
    if family == "hybrid":
        return HYBRID_COLOR
    return OTHER_COLOR


def _construction_family_proxy(name: str) -> str:
    lowered = str(name).lower()
    if "independent" in lowered:
        return "independent echo"
    if "symmetric" in lowered or lowered == "echoed":
        return "symmetric echo"
    if "asymmetric" in lowered:
        return "asymmetric echo"
    if "hybrid" in lowered:
        return "hybrid"
    return "direct"


def apply_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "font.size": 9,
            "font.family": "serif",
            "mathtext.fontset": "cm",
            "axes.grid": True,
            "grid.alpha": 0.25,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.prop_cycle": plt.cycler(
                color=[DIRECT_COLOR, SYMMETRIC_ECHO_COLOR, INDEPENDENT_ECHO_COLOR, ASYMMETRIC_ECHO_COLOR, HYBRID_COLOR]
            ),
        }
    )


def save_figure_bundle(fig: Any, stem: str, figure_data: dict[str, Any], *, figure_root: Path) -> None:
    pdf_dir = figure_root / "pdf"
    png_dir = figure_root / "png"
    data_dir = figure_root / "data"
    pdf_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(pdf_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(png_dir / f"{stem}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    (data_dir / f"{stem}.json").write_text(json.dumps(figure_data, indent=2), encoding="utf-8")


def unified_prior_work_comparison(df: pd.DataFrame, *, figure_root: Path) -> str:
    subset = df.dropna(subset=["strict_process_fidelity"]).copy()
    grouped = (
        subset.groupby(["study_short", "construction_display", "construction_family"], as_index=False)["strict_process_fidelity"]
        .max()
        .sort_values(["study_short", "strict_process_fidelity"], ascending=[True, False])
    )
    studies = list(dict.fromkeys(grouped["study_short"].tolist()))
    constructions = list(dict.fromkeys(grouped["construction_display"].tolist()))
    x = np.arange(len(studies))
    width = 0.82 / max(1, len(constructions))
    fig, ax = plt.subplots(figsize=(7.0, 4.1))
    for idx, construction in enumerate(constructions):
        rows = grouped[grouped["construction_display"] == construction]
        values = []
        colors = []
        for study in studies:
            match = rows[rows["study_short"] == study]
            if match.empty:
                values.append(np.nan)
                colors.append(OTHER_COLOR)
            else:
                values.append(float(match.iloc[0]["strict_process_fidelity"]))
                colors.append(construction_color(str(match.iloc[0]["construction_family"])))
        ax.bar(
            x + (idx - 0.5 * (len(constructions) - 1)) * width,
            values,
            width=width,
            label=construction,
            color=colors,
        )
    ax.set_xticks(x, studies, rotation=15, ha="right")
    ax.set_ylabel("Best strict process fidelity")
    ax.set_ylim(0.0, 1.05)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    save_figure_bundle(fig, "unified_prior_work_comparison", {"rows": grouped.to_dict(orient="records")}, figure_root=figure_root)
    return "unified_prior_work_comparison"


def per_manifold_error_budget(error_df: pd.DataFrame, *, figure_root: Path) -> str:
    best_cases = (
        error_df.groupby(["construction", "case_id"], as_index=False)["strict_process_fidelity"]
        .max()
        .sort_values("strict_process_fidelity", ascending=False)
        .drop_duplicates(subset=["construction"])
        .head(4)
    )
    chosen = error_df.merge(best_cases[["construction", "case_id"]], on=["construction", "case_id"], how="inner")
    constructions = list(dict.fromkeys(chosen["construction"].tolist()))
    fig, axes = plt.subplots(1, len(constructions), figsize=(3.4 * max(1, len(constructions)), 3.0), sharey=True)
    if len(constructions) == 1:
        axes = [axes]
    for ax, construction in zip(axes, constructions, strict=True):
        rows = chosen[chosen["construction"] == construction].sort_values("n")
        x = np.arange(len(rows))
        width = 0.24
        ax.bar(x - width, np.abs(rows["eps_x"]), width=width, label=r"$|\epsilon_x|$", color="#56B4E9")
        ax.bar(x, np.abs(rows["eps_y"]), width=width, label=r"$|\epsilon_y|$", color="#009E73")
        ax.bar(x + width, np.abs(rows["eps_z"]), width=width, label=r"$|\epsilon_z|$", color="#D55E00")
        ax.set_title(str(construction).replace("_", " "))
        ax.set_xticks(x, [int(v) for v in rows["n"]])
        ax.set_xlabel("Manifold $n$")
    axes[0].set_ylabel("Error-generator component (rad)")
    axes[0].legend(frameon=False, fontsize=7)
    save_figure_bundle(fig, "per_manifold_error_budget", {"rows": chosen.to_dict(orient="records")}, figure_root=figure_root)
    return "per_manifold_error_budget"


def direct_vs_echoed_scatter(df: pd.DataFrame, *, figure_root: Path) -> str:
    ideal = df[(df["target_kind"] == "ideal_sqr") & df["strict_process_fidelity"].notna()].copy()
    rows = []
    for case_id, group in ideal.groupby("case_id"):
        direct = group[group["construction_family"] == "direct"]["strict_process_fidelity"].max()
        non_direct_group = group[group["construction_family"] != "direct"]
        if pd.isna(direct) or non_direct_group.empty:
            continue
        best_row = non_direct_group.sort_values("strict_process_fidelity", ascending=False).iloc[0]
        rows.append(
            {
                "case_id": case_id,
                "direct_strict_process_fidelity": float(direct),
                "non_direct_strict_process_fidelity": float(best_row["strict_process_fidelity"]),
                "construction_family": str(best_row["construction_family"]),
            }
        )
    plot_df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    for family, group in plot_df.groupby("construction_family"):
        ax.scatter(
            group["direct_strict_process_fidelity"],
            group["non_direct_strict_process_fidelity"],
            label=family,
            color=construction_color(str(family)),
            s=28,
        )
    ax.plot([0.0, 1.0], [0.0, 1.0], linestyle="--", color="black", linewidth=1.0)
    ax.set_xlabel("Best direct strict fidelity")
    ax.set_ylabel("Best non-direct strict fidelity")
    ax.legend(frameon=False, fontsize=7)
    save_figure_bundle(fig, "direct_vs_echoed_fidelity_scatter", {"rows": rows}, figure_root=figure_root)
    return "direct_vs_echoed_fidelity_scatter"


def fidelity_vs_duration(df: pd.DataFrame, *, figure_root: Path) -> str:
    native = df[(df["study"] == "native_rich_multitone_sqr_cpsqr_feasibility") & df["strict_process_fidelity"].notna()].copy()
    native["chi_t_over_2pi"] = native["case_id"].str.extract(r"chiT([0-9p]+)")[0].str.replace("p", ".", regex=False).astype(float)
    fig, ax = plt.subplots(figsize=(6.7, 4.1))
    for family_name, group in native.groupby("construction"):
        stats = group.groupby("chi_t_over_2pi")["strict_process_fidelity"].agg(["mean", "min", "max"]).reset_index()
        color = construction_color(_construction_family_proxy(family_name))
        ax.plot(stats["chi_t_over_2pi"], stats["mean"], marker="o", label=family_name, color=color)
        ax.fill_between(stats["chi_t_over_2pi"], stats["min"], stats["max"], alpha=0.12, color=color)
    ax.axhline(0.99, linestyle="--", color="black", linewidth=1.0)
    ax.set_xlabel(r"$|\chi|T/2\pi$")
    ax.set_ylabel("Strict process fidelity")
    ax.legend(frameon=False, fontsize=6, ncol=2)
    save_figure_bundle(fig, "fidelity_vs_duration_scaling", {"rows": native.to_dict(orient="records")}, figure_root=figure_root)
    return "fidelity_vs_duration_scaling"


def ansatz_comparison_fixed_duration(df: pd.DataFrame, *, figure_root: Path) -> str:
    native = df[(df["study"] == "native_rich_multitone_sqr_cpsqr_feasibility") & df["strict_process_fidelity"].notna()].copy()
    native["chi_t_over_2pi"] = native["case_id"].str.extract(r"chiT([0-9p]+)")[0].str.replace("p", ".", regex=False).astype(float)
    longest = float(native["chi_t_over_2pi"].max())
    longest_df = native[np.isclose(native["chi_t_over_2pi"], longest)]
    grouped = longest_df.groupby("construction")["strict_process_fidelity"].agg(["mean", "std"]).reset_index()
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.bar(
        grouped["construction"],
        grouped["mean"],
        yerr=grouped["std"].fillna(0.0),
        color=[construction_color(_construction_family_proxy(v)) for v in grouped["construction"]],
    )
    ax.set_ylabel("Strict process fidelity")
    ax.set_title(rf"Longest duration slice: $|\chi|T/2\pi = {longest:g}$")
    ax.tick_params(axis="x", rotation=25)
    save_figure_bundle(fig, "ansatz_comparison_fixed_duration", {"rows": grouped.to_dict(orient="records")}, figure_root=figure_root)
    return "ansatz_comparison_fixed_duration"


def parameter_count_vs_fidelity(df: pd.DataFrame, *, figure_root: Path) -> str:
    subset = df[df["strict_process_fidelity"].notna()].copy()
    fig, ax = plt.subplots(figsize=(5.9, 4.1))
    for family, group in subset.groupby("construction_family"):
        ax.scatter(group["parameter_count"], group["strict_process_fidelity"], label=family, s=24, color=construction_color(str(family)))
    ax.set_xlabel("Approximate ansatz parameter count")
    ax.set_ylabel("Strict process fidelity")
    ax.legend(frameon=False, fontsize=7)
    save_figure_bundle(fig, "parameter_count_vs_fidelity", {"rows": subset.to_dict(orient="records")}, figure_root=figure_root)
    return "parameter_count_vs_fidelity"


def full_case_heatmap(df: pd.DataFrame, *, figure_root: Path) -> str:
    ideal = df[(df["target_kind"] == "ideal_sqr") & df["strict_process_fidelity"].notna()].copy()
    heat = (
        ideal.groupby(["case_id", "construction"], as_index=False)["strict_process_fidelity"]
        .max()
        .pivot(index="case_id", columns="construction", values="strict_process_fidelity")
        .sort_index()
    )
    fig, ax = plt.subplots(figsize=(7.0, max(4.0, 0.2 * len(heat))))
    image = ax.imshow(heat.fillna(np.nan).values, aspect="auto", vmin=0.0, vmax=1.0, cmap="viridis")
    ax.set_xticks(np.arange(len(heat.columns)), heat.columns, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(heat.index)), heat.index)
    ax.set_xlabel("Construction")
    ax.set_ylabel("Case")
    fig.colorbar(image, ax=ax, label="Strict process fidelity")
    save_figure_bundle(fig, "full_case_construction_heatmap", {"rows": heat.reset_index().to_dict(orient="records")}, figure_root=figure_root)
    return "full_case_construction_heatmap"


def strict_vs_cpsqr_joint_comparison(df: pd.DataFrame, *, figure_root: Path) -> str:
    native = df[
        (df["study"] == "native_rich_multitone_sqr_cpsqr_feasibility")
        & df["strict_process_fidelity"].notna()
        & df["cpsqr_process_fidelity"].notna()
    ].copy()
    grouped = native.groupby(["construction", "construction_family"], as_index=False)[["strict_process_fidelity", "cpsqr_process_fidelity"]].mean()
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    for _, row in grouped.iterrows():
        ax.scatter(
            row["strict_process_fidelity"],
            row["cpsqr_process_fidelity"],
            color=construction_color(str(row["construction_family"])),
            s=36,
        )
        ax.text(row["strict_process_fidelity"], row["cpsqr_process_fidelity"], str(row["construction"]), fontsize=7)
    ax.plot([0.0, 1.0], [0.0, 1.0], linestyle="--", color="black", linewidth=1.0)
    ax.set_xlabel("Mean strict joint fidelity")
    ax.set_ylabel("Mean CPSQR joint fidelity")
    save_figure_bundle(fig, "strict_vs_cpsqr_joint_comparison", {"rows": grouped.to_dict(orient="records")}, figure_root=figure_root)
    return "strict_vs_cpsqr_joint_comparison"


def one_state_vs_quartet_validation(df: pd.DataFrame, *, figure_root: Path) -> str:
    native = df[(df["study"] == "native_rich_multitone_sqr_cpsqr_feasibility") & df["strict_process_fidelity"].notna()].copy()
    grouped = (
        native.groupby("construction", as_index=False)[["reduced_quartet_fidelity", "full_quartet_fidelity", "strict_process_fidelity"]]
        .mean()
        .sort_values("strict_process_fidelity", ascending=False)
    )
    x = np.arange(len(grouped))
    width = 0.25
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.bar(x - width, grouped["reduced_quartet_fidelity"], width=width, label="Reduced quartet")
    ax.bar(x, grouped["full_quartet_fidelity"], width=width, label="Full quartet")
    ax.bar(x + width, grouped["strict_process_fidelity"], width=width, label="Joint process")
    ax.set_xticks(x, grouped["construction"], rotation=25, ha="right")
    ax.set_ylabel("Mean fidelity")
    ax.legend(frameon=False, fontsize=7)
    save_figure_bundle(fig, "one_state_vs_quartet_validation", {"rows": grouped.to_dict(orient="records")}, figure_root=figure_root)
    return "one_state_vs_quartet_validation"


def n_active_scaling(scaling_df: pd.DataFrame, scaling_payload: dict[str, Any], *, figure_root: Path) -> str:
    grouped = scaling_df.copy()
    if grouped.empty:
        return ""
    fig, ax = plt.subplots(figsize=(6.3, 4.1))
    for (construction, model_variant), group in grouped.groupby(["construction", "model_variant"]):
        ordered = group.sort_values("n_active")
        linestyle = "-" if str(model_variant) == "chi_only" else "--"
        color = construction_color(_construction_family_proxy(construction))
        label = f"{construction} ({model_variant})"
        ax.plot(ordered["n_active"], ordered["strict_process_fidelity"], marker="o", linestyle=linestyle, color=color, label=label)
    for fit in scaling_payload.get("fits", []):
        x = np.arange(2, 9, dtype=float)
        predicted = 1.0 - np.exp(float(fit["intercept"]) + float(fit["slope"]) * x)
        ax.plot(x, predicted, linestyle=":", linewidth=1.0, color=construction_color(_construction_family_proxy(str(fit["construction"]))))
    ax.set_xlabel(r"Addressed manifold count $N_{\mathrm{active}}$")
    ax.set_ylabel("Best strict process fidelity")
    ax.set_ylim(0.0, 1.02)
    ax.legend(frameon=False, fontsize=6, ncol=2)
    save_figure_bundle(
        fig,
        "n_active_scaling",
        {"rows": grouped.to_dict(orient="records"), "fits": scaling_payload.get("fits", [])},
        figure_root=figure_root,
    )
    return "n_active_scaling"


def spectral_crowding_diagram(spectral_payload: dict[str, Any], *, figure_root: Path) -> str:
    rows = pd.DataFrame(spectral_payload.get("rows", []))
    plot_rows = pd.DataFrame(spectral_payload.get("plot_rows", []))
    if rows.empty or plot_rows.empty:
        return ""
    focus_duration = 5.0
    focus_rows = rows[np.isclose(rows["chi_t_over_2pi"], focus_duration)]
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.2), sharex=False)
    for ax, model_variant in zip(axes, ["chi_only", "chi_plus_chiprime"], strict=True):
        freq_rows = plot_rows[plot_rows["model_variant"] == model_variant].sort_values("level")
        model_rows = focus_rows[focus_rows["model_variant"] == model_variant].sort_values("lower_level")
        for _, row in freq_rows.iterrows():
            freq_mhz = float(row["transition_frequency_hz"]) / 1.0e6
            ax.axvline(freq_mhz, color="black", linewidth=1.0)
            ax.text(freq_mhz, 1.02, f"$n={int(row['level'])}$", rotation=90, va="bottom", ha="center", fontsize=7)
        for _, row in model_rows.head(5).iterrows():
            center = float(freq_rows.iloc[int(row["lower_level"])]["transition_frequency_hz"]) / 1.0e6
            half_band = 0.5 * float(row["tone_bandwidth_hz"]) / 1.0e6
            ax.axvspan(center - half_band, center + half_band, color=DIRECT_COLOR, alpha=0.12)
        ax.set_ylabel(str(model_variant).replace("_", " "))
        ax.set_yticks([])
    axes[-1].set_xlabel("Transition frequency (MHz, rotating-frame offset)")
    save_figure_bundle(
        fig,
        "spectral_crowding_diagram",
        {"rows": rows.to_dict(orient="records"), "plot_rows": plot_rows.to_dict(orient="records")},
        figure_root=figure_root,
    )
    return "spectral_crowding_diagram"


def refocusing_pulse_manifold_dependence(xpi_df: pd.DataFrame, *, figure_root: Path) -> str:
    if xpi_df.empty:
        return ""
    plot_df = xpi_df[xpi_df["variant"].isin({"baseline", "robust"})].copy()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.2), sharex="col")
    variant_styles = {
        "baseline": {"label": "baseline Gaussian $X_\\pi$", "linestyle": "-", "marker": "o"},
        "robust": {"label": "compromise robust $X_\\pi$", "linestyle": "--", "marker": "s"},
    }
    for row_idx, model_variant in enumerate(["chi_only", "chi_plus_chiprime"]):
        model_rows = plot_df[plot_df["model_variant"] == model_variant].sort_values(["variant", "level"])
        ax_theta = axes[row_idx, 0]
        ax_fid = axes[row_idx, 1]
        for variant, group in model_rows.groupby("variant"):
            style = variant_styles[str(variant)]
            ax_theta.plot(group["level"], group["theta_rad"] / np.pi, color=DIRECT_COLOR if variant == "baseline" else SYMMETRIC_ECHO_COLOR, **style)
            ax_fid.plot(group["level"], group["process_fidelity"], color=DIRECT_COLOR if variant == "baseline" else SYMMETRIC_ECHO_COLOR, **style)
        ax_theta.set_ylabel(f"{model_variant}\n$\\theta_n/\\pi$")
        ax_fid.set_ylabel("Channel fidelity")
        ax_theta.set_ylim(0.0, 1.2)
        ax_fid.set_ylim(0.0, 1.02)
    axes[0, 0].legend(frameon=False, fontsize=7)
    axes[1, 0].set_xlabel("Fock manifold $n$")
    axes[1, 1].set_xlabel("Fock manifold $n$")
    save_figure_bundle(fig, "refocusing_pulse_manifold_dependence", {"rows": plot_df.to_dict(orient="records")}, figure_root=figure_root)
    return "refocusing_pulse_manifold_dependence"


def generate_all_figures(
    df: pd.DataFrame,
    error_df: pd.DataFrame,
    *,
    figure_root: Path,
    xpi_df: pd.DataFrame | None = None,
    spectral_payload: dict[str, Any] | None = None,
    scaling_payload: dict[str, Any] | None = None,
) -> list[str]:
    apply_style()
    stems: list[str] = [
        unified_prior_work_comparison(df, figure_root=figure_root),
        direct_vs_echoed_scatter(df, figure_root=figure_root),
        fidelity_vs_duration(df, figure_root=figure_root),
        ansatz_comparison_fixed_duration(df, figure_root=figure_root),
        parameter_count_vs_fidelity(df, figure_root=figure_root),
        full_case_heatmap(df, figure_root=figure_root),
        strict_vs_cpsqr_joint_comparison(df, figure_root=figure_root),
        one_state_vs_quartet_validation(df, figure_root=figure_root),
    ]
    if not error_df.empty:
        stems.append(per_manifold_error_budget(error_df, figure_root=figure_root))
    if xpi_df is not None and not xpi_df.empty:
        stem = refocusing_pulse_manifold_dependence(xpi_df, figure_root=figure_root)
        if stem:
            stems.append(stem)
    if spectral_payload:
        stem = spectral_crowding_diagram(spectral_payload, figure_root=figure_root)
        if stem:
            stems.append(stem)
    if scaling_payload:
        scaling_df = pd.DataFrame(scaling_payload.get("rows", []))
        if not scaling_df.empty:
            stem = n_active_scaling(scaling_df, scaling_payload, figure_root=figure_root)
            if stem:
                stems.append(stem)
    return stems
