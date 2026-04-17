"""
generate_appendix_figures.py

Regenerates the appendix figures for the corrected unified study.

Run from the scripts/ directory:
    python generate_appendix_figures.py
"""

from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import common as c
from run_design_space_study import plot_wigner_panels

matplotlib.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 150,
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
STUDY_ROOT = Path(__file__).resolve().parents[1]
FIG_DIR     = STUDY_ROOT / "figures"
FIG_DIR.mkdir(exist_ok=True)

CORRECTED_SUMMARY = STUDY_ROOT / "data" / "corrected_scope_summary.json"


# ---------------------------------------------------------------------------
# Fig A1 — Cavity Wigner functions from the corrected-scope study
# ---------------------------------------------------------------------------
def make_wigner_figure() -> None:
    if not CORRECTED_SUMMARY.exists():
        raise FileNotFoundError(
            "corrected_scope_summary.json is missing. Run run_design_space_study.py first."
        )
    summary = c.load_json(CORRECTED_SUMMARY)
    plot_wigner_panels(summary)
    print(f"Saved {FIG_DIR / 'appendix_wigner.pdf'}")


# ---------------------------------------------------------------------------
# Fig A2 — Qubit Bloch-vector components and cavity Fock populations (target)
# ---------------------------------------------------------------------------
def make_state_observables_figure() -> None:
    # -------------------------------------------------------------------
    # Compute analytically from U_target = SWAP · CZ · (H ⊗ I).
    #
    # Basis: {|g,0>, |g,1>, |e,0>, |e,1>} ↔ indices 0,1,2,3
    # U_target columns (= output states for each input):
    #   |g,0> → (1/√2)(|g,0> + |g,1>)      qubit=g,   cavity=(|0>+|1>)/√2
    #   |g,1> → (1/√2)(|e,0> - |e,1>)      qubit=e,   cavity=(|0>-|1>)/√2
    #   |e,0> → (1/√2)(|g,0> - |g,1>)      qubit=g,   cavity=(|0>-|1>)/√2
    #   |e,1> → (1/√2)(|e,0> + |e,1>)      qubit=e,   cavity=(|0>+|1>)/√2
    # -------------------------------------------------------------------
    U = np.array([
        [1,  0,  1,  0],
        [1,  0, -1,  0],
        [0,  1,  0,  1],
        [0, -1,  0,  1],
    ], dtype=complex) / np.sqrt(2)

    inputs = [
        (r"$|g,0\rangle$", np.array([1, 0, 0, 0], dtype=complex)),
        (r"$|g,1\rangle$", np.array([0, 1, 0, 0], dtype=complex)),
        (r"$|e,0\rangle$", np.array([0, 0, 1, 0], dtype=complex)),
        (r"$|e,1\rangle$", np.array([0, 0, 0, 1], dtype=complex)),
    ]

    # Pauli matrices (qubit subspace)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)

    def qubit_dm(psi_4):
        """Reduced density matrix for qubit, tracing out cavity."""
        # Basis: qubit ⊗ cavity → indices: |g,n> = n, |e,n> = n+2 (for N_cav=2)
        # psi shape (4,), ordering [|g,0>, |g,1>, |e,0>, |e,1>]
        # Qubit density matrix: rho_q[i,j] = sum_n <q_i,n|psi><psi|q_j,n>
        rho = np.zeros((2, 2), dtype=complex)
        # g component: coeffs 0,1 (cavity 0,1 with qubit g)
        # e component: coeffs 2,3
        amp_g = psi_4[:2]   # |g,0>, |g,1>
        amp_e = psi_4[2:]   # |e,0>, |e,1>
        rho[0, 0] = np.dot(amp_g.conj(), amp_g).real
        rho[1, 1] = np.dot(amp_e.conj(), amp_e).real
        rho[0, 1] = np.dot(amp_g.conj(), amp_e)
        rho[1, 0] = rho[0, 1].conj()
        return rho

    def cavity_populations(psi_4):
        """Fock-level populations of the cavity, tracing over qubit."""
        p0 = abs(psi_4[0])**2 + abs(psi_4[2])**2   # n=0
        p1 = abs(psi_4[1])**2 + abs(psi_4[3])**2   # n=1
        return [p0, p1]

    records = []
    for label, psi_in in inputs:
        psi_out = U @ psi_in
        rho_q   = qubit_dm(psi_out)
        expX = float(np.real(np.trace(sx @ rho_q)))
        expY = float(np.real(np.trace(sy @ rho_q)))
        expZ = float(np.real(np.trace(sz @ rho_q)))
        pops = cavity_populations(psi_out)
        records.append(dict(label=label, expX=expX, expY=expY, expZ=expZ,
                            pop0=pops[0], pop1=pops[1], psi_out=psi_out))

    # ---- Plot ----
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4))

    # Left: Bloch-vector components as grouped bar chart
    ax = axes[0]
    n = len(records)
    x = np.arange(n)
    w = 0.25
    xlabels = [r["label"] for r in records]
    ax.bar(x - w,  [r["expX"] for r in records], w, label=r"$\langle X\rangle$", color="C0")
    ax.bar(x,      [r["expY"] for r in records], w, label=r"$\langle Y\rangle$", color="C1")
    ax.bar(x + w,  [r["expZ"] for r in records], w, label=r"$\langle Z\rangle$", color="C2")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels)
    ax.set_ylabel("Expectation value")
    ax.set_title("Transmon Bloch-vector components\n"
                 r"(output of $U_{\rm target}$, exact)")
    ax.legend(fontsize=8)
    ax.set_ylim(-1.2, 1.4)
    for xi, r in zip(x, records):
        ax.text(xi + w, r["expZ"] + 0.05, f"{r['expZ']:+.2f}", ha="center",
                va="bottom", fontsize=7)

    # Right: cavity Fock populations per sector
    ax2 = axes[1]
    w2 = 0.35
    ax2.bar(x - w2/2, [r["pop0"] for r in records], w2,
            label=r"$P(n=0)$", color="C3", alpha=0.85)
    ax2.bar(x + w2/2, [r["pop1"] for r in records], w2,
            label=r"$P(n=1)$", color="C4", alpha=0.85)
    ax2.set_xticks(x)
    ax2.set_xticklabels(xlabels)
    ax2.set_ylabel("Population")
    ax2.set_title("Cavity Fock populations\n"
                  r"(output of $U_{\rm target}$, exact)")
    ax2.legend(fontsize=8)
    ax2.set_ylim(0, 0.65)
    for xi, r in zip(x, records):
        ax2.text(xi - w2/2, r["pop0"] + 0.01, f"{r['pop0']:.2f}", ha="center",
                 va="bottom", fontsize=7, color="C3")
        ax2.text(xi + w2/2, r["pop1"] + 0.01, f"{r['pop1']:.2f}", ha="center",
                 va="bottom", fontsize=7, color="C4")

    fig.suptitle(
        "Appendix A2 — Target output state observables for each logical input\n"
        r"$U_{\rm target}=\mathrm{SWAP}\cdot\mathrm{CZ}\cdot(H\otimes I)$, "
        r"basis $\{|g,0\rangle,|g,1\rangle,|e,0\rangle,|e,1\rangle\}$",
        fontsize=10,
    )
    fig.tight_layout()
    out = FIG_DIR / "appendix_state_observables.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    print(f"Saved {out}")
    plt.close(fig)


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Generating appendix figures ...")
    make_wigner_figure()
    make_state_observables_figure()
    print("Done.")
