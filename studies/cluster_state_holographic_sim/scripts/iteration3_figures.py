"""
Iteration 3 figures for cluster-state holographic simulation study.

Generates:
  1. hilbert_space_validity.{png,pdf}  – N_cav=2 vs N_cav=8 fidelity comparison
  2. hilbert_convergence.{png,pdf}     – Fidelity vs Hilbert space dimension
  3. combined_fidelity_budget.{png,pdf} – GRAPE fidelity × coherence across all durations
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
FIG_DIR = BASE / "figures"
DATA_DIR = BASE / "data"
FIG_DIR.mkdir(exist_ok=True)

# ── Style ────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "legend.fontsize": 9,
    "figure.dpi": 150,
})

# ── Data from iteration 3 evaluation (terminal output, confirmed) ───────
# Strategy B: D+SQR+CP, 2 blocks
eval_B = {
    "ideal_ncav2": 1.000000,
    "embed_ncav8_nodrift": 0.111825,
    "model_ncav8_drift": 0.093660,
    "model_ncav12_drift": 0.074658,
    "model_ncav15_drift": 0.074458,
    "leak_ncav8_nodrift": 0.919084,
    "leak_ncav8_drift": 0.933435,
}

# Strategy D: D+R+FE, 2 blocks
eval_D = {
    "ideal_ncav2": 0.999869,
    "embed_ncav8": 0.187811,
    "embed_ncav12": 0.274515,
    "leak_ncav8": 0.744690,
}

# FE wait-time analysis
TAU_CZ_NS = 176.1
FE0_NS = 177.0
FE1_NS = 177.0

# GRAPE sweep (from results_combined.json)
grape_data = [
    (50,  0.6337), (100, 0.9494), (150, 0.9561),
    (200, 0.9966), (300, 0.9957), (400, 0.9990),
    (500, 0.99999), (600, 1.0000), (800, 1.0000),
]

# Coherence budget data
with open(DATA_DIR / "iteration3_coherence_budget.json") as f:
    coh_data = json.load(f)
budget = coh_data["budget"]
T1_q = coh_data["parameters"]["T1_qubit_us"]
T2_q = coh_data["parameters"]["T2_qubit_us"]

# ═══════════════════════════════════════════════════════════════════════════
# Figure 1: Hilbert Space Validity – N_cav=2 vs N_cav=8
# ═══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Left panel: fidelity comparison
strategies = ["D+SQR+CP\n(2 blocks)", "D+R+FE\n(2 blocks)", "GRAPE\n400 ns"]
fid_ncav2 = [eval_B["ideal_ncav2"], eval_D["ideal_ncav2"], 0.999]
fid_ncav8 = [eval_B["model_ncav8_drift"], eval_D["embed_ncav8"], 0.999]  # GRAPE is model-based already

x = np.arange(len(strategies))
w = 0.35
bars1 = axes[0].bar(x - w/2, fid_ncav2, w, label=r"$N_\mathrm{cav}=2$ (ideal)", color="#2196F3", alpha=0.85)
bars2 = axes[0].bar(x + w/2, fid_ncav8, w, label=r"$N_\mathrm{cav}=8$ (model)", color="#F44336", alpha=0.85)

axes[0].set_ylabel("Subspace Fidelity")
axes[0].set_title("Fidelity: Ideal Mode vs. Full Model")
axes[0].set_xticks(x)
axes[0].set_xticklabels(strategies)
axes[0].set_ylim(0, 1.15)
axes[0].legend(loc="upper left")
axes[0].axhline(0.99, color="gray", ls="--", lw=0.8, alpha=0.5)
axes[0].text(2.4, 0.995, r"$\mathcal{F}=0.99$", fontsize=8, color="gray")

# Add value labels
for bar in bars1:
    h = bar.get_height()
    axes[0].text(bar.get_x() + bar.get_width()/2, h + 0.02,
                 f"{h:.3f}", ha="center", va="bottom", fontsize=8)
for bar in bars2:
    h = bar.get_height()
    axes[0].text(bar.get_x() + bar.get_width()/2, h + 0.02,
                 f"{h:.3f}", ha="center", va="bottom", fontsize=8)

# Right panel: leakage comparison
strategies_leak = ["D+SQR+CP\n(no drift)", "D+SQR+CP\n(full drift)", "D+R+FE"]
leakages = [eval_B["leak_ncav8_nodrift"], eval_B["leak_ncav8_drift"], eval_D["leak_ncav8"]]
colors = ["#FF9800", "#F44336", "#9C27B0"]
bars = axes[1].bar(strategies_leak, leakages, color=colors, alpha=0.85, width=0.5)
axes[1].set_ylabel("Average Leakage")
axes[1].set_title(r"Leakage at $N_\mathrm{cav}=8$")
axes[1].set_ylim(0, 1.1)

for bar in bars:
    h = bar.get_height()
    axes[1].text(bar.get_x() + bar.get_width()/2, h + 0.02,
                 f"{h:.1%}", ha="center", va="bottom", fontsize=9, fontweight="bold")

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"hilbert_space_validity.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  [1/3] hilbert_space_validity.{png,pdf}")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2: Hilbert Space Convergence – Fidelity vs N_cav
# ═══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(6, 4))

ncav_B = [2, 8, 12, 15]
fid_B = [eval_B["ideal_ncav2"], eval_B["model_ncav8_drift"],
         eval_B["model_ncav12_drift"], eval_B["model_ncav15_drift"]]

ncav_D = [2, 8, 12]
fid_D = [eval_D["ideal_ncav2"], eval_D["embed_ncav8"], eval_D["embed_ncav12"]]

ax.plot(ncav_B, fid_B, "o-", color="#2196F3", lw=2, ms=8, label="D+SQR+CP (2 blocks)")
ax.plot(ncav_D, fid_D, "s--", color="#9C27B0", lw=2, ms=8, label="D+R+FE (2 blocks)")

# GRAPE reference line (model-based at N_cav=8, no truncation issue)
ax.axhline(0.999, color="#4CAF50", ls="-.", lw=1.5, alpha=0.7, label="GRAPE 400 ns (model-based)")

ax.set_xlabel(r"Cavity Hilbert Space Dimension $N_\mathrm{cav}$")
ax.set_ylabel("Subspace Fidelity")
ax.set_title("Hilbert Space Convergence of Parametric Decompositions")
ax.set_xticks([2, 4, 6, 8, 10, 12, 15])
ax.set_ylim(-0.05, 1.1)
ax.legend(loc="center right")
ax.axhline(0.5, color="gray", ls=":", lw=0.8, alpha=0.4)
ax.text(14.5, 0.52, "random", fontsize=7, color="gray")

# Annotate the crash
ax.annotate(r"$>90\%$ leakage", xy=(8, eval_B["model_ncav8_drift"]),
            xytext=(10, 0.4), fontsize=9, color="#F44336",
            arrowprops=dict(arrowstyle="->", color="#F44336", lw=1.5))

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"hilbert_convergence.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  [2/3] hilbert_convergence.{png,pdf}")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 3: Combined Fidelity Budget (unitary × coherence)
# ═══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Left: GRAPE unitary fidelity + coherence limit + combined
durations = [g[0] for g in grape_data]
fid_unitary = [g[1] for g in grape_data]

# Coherence limit: F_coh = exp(-t/T1 - t/T2) / 2 (simplified; use budget data where available)
# More precisely from the analytical formula
T1_s = T1_q * 1e-6
T2_s = T2_q * 1e-6
fid_coherence = [np.exp(-t*1e-9 * (1/(2*T1_s) + 1/(2*T2_s))) for t in durations]

fid_combined = [fu * fc for fu, fc in zip(fid_unitary, fid_coherence)]

ax = axes[0]
ax.semilogy(durations, [1-f for f in fid_unitary], "o-", color="#2196F3", lw=2, ms=6,
            label=r"Unitary infidelity $1-\mathcal{F}_U$")
ax.semilogy(durations, [1-f for f in fid_coherence], "s--", color="#FF9800", lw=2, ms=6,
            label=r"Decoherence limit $1-\mathcal{F}_\mathrm{coh}$")
ax.semilogy(durations, [max(1-f, 1e-10) for f in fid_combined], "D-", color="#4CAF50", lw=2, ms=6,
            label=r"Combined $1-\mathcal{F}_U \cdot \mathcal{F}_\mathrm{coh}$")

ax.set_xlabel("Gate Duration (ns)")
ax.set_ylabel("Infidelity")
ax.set_title("GRAPE: Unitary vs. Coherence Tradeoff")
ax.legend(fontsize=8)
ax.set_ylim(1e-8, 1)
ax.grid(True, alpha=0.3)

# Mark optimal region
ax.axvspan(150, 450, alpha=0.08, color="green")
ax.text(300, 3e-8, "Sweet spot", ha="center", fontsize=8, color="green")

# Right: All strategies combined fidelity comparison
ax = axes[1]
labels = []
combined_vals = []
colors_bar = []

# Parametric strategies (at N_cav=2, ideal — as disclaimer)
labels.append("D+SQR+CP\n(ideal, N=2)")
combined_vals.append(budget["B_D+SQR+CP_2blk_ideal"]["F_combined"])
colors_bar.append("#BBDEFB")  # light blue to indicate "ideal only"

labels.append("D+R+FE\n(ideal, N=2)")
combined_vals.append(budget["D_D+R+FE_2blk_ideal"]["F_combined"])
colors_bar.append("#E1BEE7")  # light purple

# GRAPE strategies (model-based, validated)
for dns in [200, 300, 400]:
    key = f"GRAPE_{dns}ns"
    labels.append(f"GRAPE\n{dns} ns")
    combined_vals.append(budget[key]["F_combined"])
    colors_bar.append("#4CAF50")

bars = ax.bar(labels, combined_vals, color=colors_bar, alpha=0.85, width=0.6,
              edgecolor="gray", linewidth=0.5)
ax.set_ylabel(r"Combined Fidelity $\mathcal{F}_U \cdot \mathcal{F}_\mathrm{coh}$")
ax.set_title("Strategy Comparison with Decoherence")
ax.set_ylim(0.8, 1.01)

for bar in bars:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, h + 0.002,
            f"{h:.3f}", ha="center", va="bottom", fontsize=8)

# Add asterisks for ideal-mode disclaimer
ax.text(0, 0.82, "*ideal mode only;\nfails at N_cav=8", fontsize=7, color="#F44336",
        ha="center")
ax.text(1, 0.82, "*ideal mode only;\nfails at N_cav=8", fontsize=7, color="#F44336",
        ha="center")

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"combined_fidelity_budget.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  [3/3] combined_fidelity_budget.{png,pdf}")

print("\nAll iteration 3 figures generated.")
