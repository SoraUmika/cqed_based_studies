"""Generate publication-quality figures for the open-system SQR deep-dive study."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_DIR = STUDY_DIR / "data"
FIGURES_DIR = STUDY_DIR / "figures"
STYLE_PATH = SCRIPT_DIR.parents[2] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"

NOISE_PATH = DATA_DIR / "noise_channel_sweeps.npz"
PURCELL_PATH = DATA_DIR / "purcell_and_backaction.npz"
THREE_MODE_PATH = DATA_DIR / "three_mode_readout_effects.npz"
GRAPE_PATH = DATA_DIR / "grape_noisy_replay.npz"
CONVERGENCE_PATH = DATA_DIR / "convergence" / "a2_convergence.npz"
CONVERGENCE_TARGET = 5.0e-4

REPRESENTATIVE_T1_US = np.array([30.0, 10.0], dtype=float)
REPRESENTATIVE_CHI_T = 2.0
THREE_MODE_REPRESENTATIVE_CHI_T = 3.0
REPRESENTATIVE_THERMAL_OCCUPATION = 0.02

FAMILY_LABELS = {
    "single_tone_gaussian": "Single-tone Gaussian",
    "square": "Square",
    "cosine_squared": "Cosine-squared",
    "multitone_one_segment": "Multitone",
}

FIGURES_DIR.mkdir(parents=True, exist_ok=True)
plt.style.use(str(STYLE_PATH))


def load_payload(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"{path} not found.")
    return np.load(path, allow_pickle=True)


def family_color_map(family_names: list[str]) -> dict[str, str]:
    palette = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]
    return {family_name: palette[idx % len(palette)] for idx, family_name in enumerate(family_names)}


def family_label(name: str) -> str:
    return FAMILY_LABELS.get(name, name.replace("_", " ").title())


def match_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(np.asarray(values, dtype=float) - float(target))))


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {stem}.png/.pdf")


def representative_family_index(noise_payload) -> int:
    chi_t_values = np.asarray(noise_payload["chi_t_values"], dtype=float)
    chi_mask = (chi_t_values >= 1.0) & (chi_t_values <= 5.0)
    t1_cases = np.asarray(noise_payload["multilevel_t1_cases_us"], dtype=float)
    case_index = int(np.argmin(np.sum(np.abs(t1_cases - REPRESENTATIVE_T1_US), axis=1)))
    multilevel_target = np.asarray(noise_payload["multilevel_fidelity_to_target"], dtype=float)
    scores = np.mean(multilevel_target[case_index][:, chi_mask], axis=1)
    return int(np.argmax(scores))


def plot_legacy_and_multilevel(noise_payload) -> None:
    family_names = [str(name) for name in noise_payload["family_names"]]
    chi_t_values = np.asarray(noise_payload["chi_t_values"], dtype=float)
    legacy_target = np.asarray(noise_payload["legacy_fidelity_to_target"], dtype=float)
    multilevel_target = np.asarray(noise_payload["multilevel_fidelity_to_target"], dtype=float)
    t1_cases = np.asarray(noise_payload["multilevel_t1_cases_us"], dtype=float)
    t1_case_index = int(np.argmin(np.sum(np.abs(t1_cases - REPRESENTATIVE_T1_US), axis=1)))
    colors = family_color_map(family_names)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    ax = axes[0]
    for family_index, family_name in enumerate(family_names):
        ax.plot(
            chi_t_values,
            legacy_target[family_index],
            marker="o",
            label=family_label(family_name),
            color=colors[family_name],
        )
    ax.set_xlabel(r"$|\chi_s| T / 2\pi$")
    ax.set_ylabel("Mean fidelity to target")
    ax.set_ylim(0.5, 1.0)
    ax.set_title("Legacy open-system SQR fidelity window")
    ax.legend(ncol=2)

    ax = axes[1]
    for family_index, family_name in enumerate(family_names):
        penalty = multilevel_target[t1_case_index, family_index] - legacy_target[family_index]
        ax.plot(
            chi_t_values,
            penalty,
            marker="o",
            label=family_label(family_name),
            color=colors[family_name],
        )
    ax.axhline(0.0, color="#666666", linewidth=1.0, linestyle="--")
    ax.set_xlabel(r"$|\chi_s| T / 2\pi$")
    ax.set_ylabel("Fidelity shift vs legacy model")
    ax.set_title(r"Explicit multilevel $T_1$ penalty at $(30, 10)\,\mu$s")

    fig.tight_layout()
    save_figure(fig, "fig1_legacy_and_multilevel")


def plot_thermal_and_branch(noise_payload) -> None:
    family_names = [str(name) for name in noise_payload["family_names"]]
    chi_t_values = np.asarray(noise_payload["chi_t_values"], dtype=float)
    thermal_target = np.asarray(noise_payload["thermal_fidelity_to_target"], dtype=float)
    thermal_occupations = np.asarray(noise_payload["thermal_occupations"], dtype=float)
    legacy_branch = np.asarray(noise_payload["legacy_branch_to_target"], dtype=float)
    rep_family_idx = representative_family_index(noise_payload)
    rep_family_name = family_names[rep_family_idx]
    rep_chi_idx = match_index(chi_t_values, REPRESENTATIVE_CHI_T)
    colors = family_color_map(family_names)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    ax = axes[0]
    for family_index, family_name in enumerate(family_names):
        ax.plot(
            thermal_occupations,
            thermal_target[:, family_index, rep_chi_idx],
            marker="o",
            label=family_label(family_name),
            color=colors[family_name],
        )
    ax.set_xlabel(r"Storage thermal occupation $n_{\mathrm{th},s}$")
    ax.set_ylabel("Mean fidelity to target")
    ax.set_ylim(0.5, 1.0)
    ax.set_title(rf"Thermal sensitivity at $|\chi_s|T/2\pi = {chi_t_values[rep_chi_idx]:.1f}$")
    ax.legend(ncol=2)

    ax = axes[1]
    branch_indices = np.arange(legacy_branch.shape[-1], dtype=int)
    chi_indices = [match_index(chi_t_values, value) for value in (1.0, 2.0, 3.0, 5.0)]
    for line_index, chi_index in enumerate(chi_indices):
        ax.plot(
            branch_indices,
            legacy_branch[rep_family_idx, chi_index],
            marker="o",
            linewidth=1.8,
            label=rf"$|\chi_s|T/2\pi = {chi_t_values[chi_index]:.1f}$",
            color=["#4477AA", "#EE6677", "#228833", "#AA3377"][line_index],
        )
    ax.axvline(1, color="#666666", linewidth=1.0, linestyle="--")
    ax.set_xlabel("Storage branch")
    ax.set_ylabel("Branch-averaged target fidelity")
    ax.set_ylim(0.5, 1.0)
    ax.set_xticks(branch_indices)
    ax.set_title(f"Selectivity profile for {family_label(rep_family_name)}")
    ax.legend()

    fig.tight_layout()
    save_figure(fig, "fig2_thermal_and_branch")


def plot_purcell_and_backaction(purcell_payload) -> None:
    detuning = np.asarray(purcell_payload["detuning_ghz"], dtype=float)
    durations = np.asarray(purcell_payload["backaction_durations_ns"], dtype=float)
    epsilon_values = np.asarray(purcell_payload["readout_epsilon_mhz"], dtype=float)
    t1_no_filter = np.asarray(purcell_payload["purcell_t1_no_filter_s"], dtype=float) * 1.0e6
    t1_with_filter = np.asarray(purcell_payload["purcell_t1_with_filter_s"], dtype=float) * 1.0e6
    coherence = np.asarray(purcell_payload["coherence_after_backaction"], dtype=float)
    excited = np.asarray(purcell_payload["excited_population_after_backaction"], dtype=float)
    palette = ["#4477AA", "#EE6677", "#228833", "#AA3377"]

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.0))

    ax = axes[0]
    ax.plot(detuning, t1_no_filter, marker="o", label="No filter", color="#EE6677")
    ax.plot(detuning, t1_with_filter, marker="o", label="With Purcell filter", color="#4477AA")
    ax.set_xlabel(r"$|\omega_q - \omega_r| / 2\pi$ (GHz)")
    ax.set_ylabel(r"Purcell-limited $T_1$ ($\mu$s)")
    ax.set_yscale("log")
    ax.set_title("Purcell protection from the readout chain")
    ax.legend()

    ax = axes[1]
    for idx, epsilon_mhz in enumerate(epsilon_values):
        ax.plot(durations, coherence[idx], marker="o", label=rf"$\epsilon/2\pi = {epsilon_mhz:.1f}$ MHz", color=palette[idx % len(palette)])
    ax.set_xlabel("Interaction duration (ns)")
    ax.set_ylabel(r"$|\rho_{ge}|$ after backaction")
    ax.set_title("Measurement-induced dephasing")
    ax.legend(fontsize=9)

    ax = axes[2]
    for idx, epsilon_mhz in enumerate(epsilon_values):
        ax.plot(durations, excited[idx], marker="o", label=rf"$\epsilon/2\pi = {epsilon_mhz:.1f}$ MHz", color=palette[idx % len(palette)])
    ax.set_xlabel("Interaction duration (ns)")
    ax.set_ylabel(r"Excited-state population $P_e$")
    ax.set_title("Readout-induced relaxation")

    fig.tight_layout()
    save_figure(fig, "fig3_purcell_and_backaction")


def plot_three_mode(three_mode_payload, noise_payload) -> None:
    family_names = [str(name) for name in three_mode_payload["family_names"]]
    chi_t_values = np.asarray(three_mode_payload["chi_t_values"], dtype=float)
    readout_amplitudes = np.asarray(three_mode_payload["readout_amplitudes_mhz"], dtype=float)
    reduced = np.asarray(three_mode_payload["reduced_target_fidelity"], dtype=float)
    coherence_ratio = np.asarray(three_mode_payload["storage_coherence_ratio"], dtype=float)
    rep_family_idx = representative_family_index(noise_payload)
    rep_family_name = [str(name) for name in noise_payload["family_names"]][rep_family_idx]
    rep_three_mode_family_idx = family_names.index(rep_family_name) if rep_family_name in family_names else 0
    rep_chi_idx = match_index(chi_t_values, THREE_MODE_REPRESENTATIVE_CHI_T)
    colors = family_color_map(family_names)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    ax = axes[0]
    for family_index, family_name in enumerate(family_names):
        ax.plot(
            readout_amplitudes,
            reduced[family_index, rep_chi_idx],
            marker="o",
            label=family_label(family_name),
            color=colors[family_name],
        )
    ax.set_xlabel(r"Readout drive amplitude $\epsilon / 2\pi$ (MHz)")
    ax.set_ylabel("Reduced qubit-storage fidelity")
    ax.set_ylim(0.0, 1.0)
    ax.set_title(rf"Three-mode spectator fidelity at $|\chi_s|T/2\pi = {chi_t_values[rep_chi_idx]:.1f}$")
    ax.legend()

    ax = axes[1]
    for chi_index, chi_t_value in enumerate(chi_t_values):
        ax.plot(
            readout_amplitudes,
            coherence_ratio[rep_three_mode_family_idx, chi_index],
            marker="o",
            label=rf"$|\chi_s|T/2\pi = {chi_t_value:.1f}$",
        )
    ax.set_xlabel(r"Readout drive amplitude $\epsilon / 2\pi$ (MHz)")
    ax.set_ylabel("Storage coherence ratio")
    ax.set_ylim(0.0, 1.05)
    ax.set_title(f"Readout-induced storage dephasing for {family_label(family_names[rep_three_mode_family_idx])}")
    ax.legend(fontsize=9)

    fig.tight_layout()
    save_figure(fig, "fig4_three_mode_readout")


def plot_grape_vs_parametric(noise_payload, grape_payload) -> None:
    family_names = [str(name) for name in noise_payload["family_names"]]
    chi_t_values = np.asarray(noise_payload["chi_t_values"], dtype=float)
    thermal_target = np.asarray(noise_payload["thermal_fidelity_to_target"], dtype=float)
    thermal_occupations = np.asarray(noise_payload["thermal_occupations"], dtype=float)
    rep_family_idx = representative_family_index(noise_payload)
    rep_family_name = family_names[rep_family_idx]
    nth_index = match_index(thermal_occupations, REPRESENTATIVE_THERMAL_OCCUPATION)
    parametric_trace = thermal_target[nth_index, rep_family_idx]

    grape_chi_t = np.asarray(grape_payload["chi_t_values"], dtype=float)
    grape_objective = np.asarray(grape_payload["objective_fidelity"], dtype=float)
    grape_noisy = np.asarray(grape_payload["noisy_fidelity_to_target"], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    ax = axes[0]
    ax.plot(
        chi_t_values,
        parametric_trace,
        marker="o",
        linewidth=2.0,
        label=rf"{family_label(rep_family_name)} with $n_{{\mathrm{{th}},s}} = 0.02$",
        color="#4477AA",
    )
    ax.plot(
        grape_chi_t,
        grape_noisy,
        marker="s",
        linewidth=2.0,
        label="GRAPE noisy replay",
        color="#EE6677",
    )
    ax.set_xlabel(r"$|\chi_s| T / 2\pi$")
    ax.set_ylabel("Mean fidelity to target")
    ax.set_ylim(0.5, 1.0)
    ax.set_title("Parametric baseline versus GRAPE under realistic noise")
    ax.legend()

    ax = axes[1]
    ax.plot(grape_chi_t, grape_objective, marker="o", linewidth=2.0, label="GRAPE objective fidelity", color="#228833")
    ax.plot(grape_chi_t, grape_noisy, marker="s", linewidth=2.0, label="GRAPE noisy target fidelity", color="#EE6677")
    parametric_interp = np.interp(grape_chi_t, chi_t_values, parametric_trace)
    ax.plot(grape_chi_t, grape_noisy - parametric_interp, marker="^", linewidth=2.0, label="GRAPE minus parametric", color="#AA3377")
    ax.axhline(0.0, color="#666666", linewidth=1.0, linestyle="--")
    ax.set_xlabel(r"$|\chi_s| T / 2\pi$")
    ax.set_ylabel("Fidelity / fidelity gap")
    ax.set_title("Closed-system advantage versus noisy realized gain")
    ax.legend(fontsize=9)

    fig.tight_layout()
    save_figure(fig, "fig5_grape_vs_parametric")


def plot_convergence(convergence_payload) -> None:
    labels = [r"Two-mode $\Delta t$", "Two-mode trunc.", r"Three-mode $\Delta t$", "Three-mode trunc.", "GRAPE objective", "GRAPE noisy"]
    values = np.array(
        [
            float(convergence_payload["two_mode_dt_delta"]),
            float(convergence_payload["two_mode_dims_delta"]),
            float(convergence_payload["three_mode_dt_delta"]),
            float(convergence_payload["three_mode_dims_delta"]),
            float(convergence_payload["grape_objective_delta"]),
            float(convergence_payload["grape_noisy_delta"]),
        ],
        dtype=float,
    )
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.bar(np.arange(values.size), values, color=["#4477AA", "#66CCEE", "#EE6677", "#CCBB44", "#228833", "#AA3377"])
    ax.axhline(CONVERGENCE_TARGET, color="#666666", linestyle="--", linewidth=1.0, label=r"$5\times 10^{-4}$ target")
    ax.set_yscale("log")
    ax.set_xticks(np.arange(values.size))
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("Absolute metric change")
    ax.set_title("Representative convergence checks")
    ax.legend()
    fig.tight_layout()
    save_figure(fig, "fig6_convergence")


def main() -> None:
    noise = load_payload(NOISE_PATH)
    purcell = load_payload(PURCELL_PATH)
    three_mode = load_payload(THREE_MODE_PATH)
    grape = load_payload(GRAPE_PATH)

    plot_legacy_and_multilevel(noise)
    plot_thermal_and_branch(noise)
    plot_purcell_and_backaction(purcell)
    plot_three_mode(three_mode, noise)
    plot_grape_vs_parametric(noise, grape)
    if CONVERGENCE_PATH.exists():
        convergence = load_payload(CONVERGENCE_PATH)
        plot_convergence(convergence)


if __name__ == "__main__":
    main()