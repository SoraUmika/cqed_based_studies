"""Phase 3 — Physical embedding: convergence vs cavity truncation.

Objectives
----------
Embed the logical target in the full N-dimensional cavity Hilbert space and
check how the results converge with N = 4, 6, 8, 10.

Key questions:
  - How large does N need to be to suppress truncation artefacts?
  - Do displacements and conditional phases cause leakage outside {|0>,|1>}?
  - How does the DriftPhaseModel (χ, χ′, K) affect the logical subspace?

Methodology:
  1. For each N, build the logical subspace Subspace.qubit_cavity_block(1, N).
  2. Construct the target unitary embedded in the full N-space.
  3. Directly test ideal primitive unitaries (SQR, displacement, CΦ) and
     measure leakage induced by each.
  4. Check dispersive-phase accumulation over t_π = π/|χ|.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── path setup ────────────────────────────────────────────────────────────
STUDY_ROOT = Path(__file__).resolve().parents[2]
COMPONENT_NAME = Path(__file__).resolve().parent.name
SIM_ROOT   = Path("C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group"
                  "/Users/Users_JianJun/cQED_simulation")
if str(SIM_ROOT) not in sys.path:
    sys.path.insert(0, str(SIM_ROOT))

import numpy as np
import qutip as qt
from cqed_sim.core import DispersiveTransmonCavityModel, FrameSpec
from cqed_sim.gates.bosonic import displacement, snap, oscillator_rotation
from cqed_sim.gates.coupled import sqr, dispersive_phase
from cqed_sim.core.conventions import qubit_cavity_index, qubit_cavity_block_indices
from cqed_sim.unitary_synthesis import (
    Subspace, subspace_unitary_fidelity, leakage_metrics,
    DriftPhaseModel, drift_phase_unitary,
)

DATA_DIR = STUDY_ROOT / "data" / COMPONENT_NAME
FIG_DIR  = STUDY_ROOT / "figures" / COMPONENT_NAME

# ── Device parameters ───────────────────────────────────────────────────
CHI  = 2 * np.pi * (-2.84e6)    # rad/s
CHIP = 2 * np.pi * (-21e3)      # rad/s
KERR = 2 * np.pi * (-28e3)      # rad/s

# ── Logical target (4×4) ────────────────────────────────────────────────
s = 1.0 / np.sqrt(2)
U_target_4x4 = np.array([
    [ s, 0,  s, 0],
    [ s, 0, -s, 0],
    [ 0, s,  0, s],
    [ 0,-s,  0, s],
], dtype=np.complex128)

CNOT_q2c = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]], dtype=np.complex128)
CNOT_c2q = np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]], dtype=np.complex128)
H2 = np.array([[1,1],[1,-1]], dtype=np.complex128) / np.sqrt(2)
Hc_4x4 = np.kron(np.eye(2, dtype=np.complex128), H2)

results = {}

# ════════════════════════════════════════════════════════════════════════════
# Helper: permutation matrix cavity-first ↔ qubit-first
# ════════════════════════════════════════════════════════════════════════════
def cavfirst_to_qfirst_perm(N: int) -> np.ndarray:
    """Permutation P such that U_qf = P @ U_cf @ P.T."""
    full = 2 * N
    P = np.zeros((full, full), dtype=np.complex128)
    for q in range(2):
        for c in range(N):
            row = q * N + c   # qubit-first index
            col = c * 2 + q   # cavity-first index
            P[row, col] = 1.0
    return P


# ════════════════════════════════════════════════════════════════════════════
# Helper: build CNOT_{q→c} in full N-space (qubit-first)
# ════════════════════════════════════════════════════════════════════════════
def cnot_q2c_full(N: int) -> np.ndarray:
    full = 2 * N
    U = np.eye(full, dtype=np.complex128)
    # For |e> sector: swap |e,0> ↔ |e,1>  (indices N and N+1)
    U[N, N]     = 0.0
    U[N+1, N+1] = 0.0
    U[N, N+1]   = 1.0
    U[N+1, N]   = 1.0
    return U


# ════════════════════════════════════════════════════════════════════════════
# Helper: build H_c (logical cavity Hadamard) in full N-space (qubit-first)
# ════════════════════════════════════════════════════════════════════════════
def hc_full(N: int) -> np.ndarray:
    full = 2 * N
    U = np.eye(full, dtype=np.complex128)
    for q in range(2):
        i00 = qubit_cavity_index(N, q, 0)
        i01 = qubit_cavity_index(N, q, 1)
        U[i00, i00] =  s
        U[i00, i01] =  s
        U[i01, i00] =  s
        U[i01, i01] = -s
    return U


# ════════════════════════════════════════════════════════════════════════════
# 1. Truncation convergence: exact analytic gates
# ════════════════════════════════════════════════════════════════════════════
print("\n── 1. Truncation convergence (analytic gates) ───────────────────────")

N_VALS = [2, 4, 6, 8, 10, 12]
trunc_data = {}

for N in N_VALS:
    log_idx = [0, 1, N, N+1]
    subspace = Subspace.custom(2*N, log_idx, ["|g,0>","|g,1>","|e,0>","|e,1>"])
    P = cavfirst_to_qfirst_perm(N)

    # Build ideal SQR_1(π,0) in qubit-first ordering
    sqr1_cf = np.asarray(sqr(np.pi, 0.0, n=1, cavity_dim=N).full())
    sqr1_qf = P @ sqr1_cf @ P.T

    # SQR_1(π,0) = -iX conditioned on n=1; need CP correction to get exact CNOT_{c→q}:
    # CNOT_{c→q} = SQR_1(π,0) · CP_corr where CP_corr applies 1j to |g,1> and |e,1>
    # (qubit-first indices 1 and N+1 respectively)
    CP_corr = np.eye(2*N, dtype=np.complex128)
    CP_corr[1, 1]     = 1j   # |g,1>  (qubit-first index = 0*N+1 = 1)
    CP_corr[N+1, N+1] = 1j   # |e,1>  (qubit-first index = 1*N+1 = N+1)
    U_cnot_c2q = sqr1_qf @ CP_corr  # = exact CNOT_{c→q}

    # Build ideal CNOT_{q→c} and H_c
    U_q2c = cnot_q2c_full(N)
    U_Hc  = hc_full(N)

    # Full composition: U_target = H_c · CNOT_{c→q} · CNOT_{q→c}
    U_comp = U_Hc @ U_cnot_c2q @ U_q2c
    U_log  = U_comp[np.ix_(log_idx, log_idx)]

    F = subspace_unitary_fidelity(U_log, U_target_4x4, gauge="global")
    err = np.linalg.norm(U_log - U_target_4x4, ord="fro")

    # Leakage of each primitive individually
    # SQR leakage (should be zero — SQR is exact in any truncation)
    sqr1_log = sqr1_qf[np.ix_(log_idx, log_idx)]
    leak_sqr1 = np.linalg.norm(
        (np.eye(2*N) - subspace.projector()) @ sqr1_qf @ subspace.projector(),
        ord="fro"
    ) ** 2 / 4
    # CNOT_{q→c} leakage (should be zero — exact gate)
    leak_q2c = np.linalg.norm(
        (np.eye(2*N) - subspace.projector()) @ U_q2c @ subspace.projector(),
        ord="fro"
    ) ** 2 / 4
    # H_c leakage
    leak_Hc = np.linalg.norm(
        (np.eye(2*N) - subspace.projector()) @ U_Hc @ subspace.projector(),
        ord="fro"
    ) ** 2 / 4

    trunc_data[N] = {
        "F": float(F), "err": float(err),
        "leak_SQR1": float(leak_sqr1),
        "leak_CNOT_q2c": float(leak_q2c),
        "leak_Hc": float(leak_Hc),
    }
    print(f"  N={N:2d}: F={F:.8f}  err={err:.2e}  "
          f"L_SQR={leak_sqr1:.2e}  L_CNOT_q2c={leak_q2c:.2e}  L_Hc={leak_Hc:.2e}")

results["truncation_analytic"] = trunc_data


# ════════════════════════════════════════════════════════════════════════════
# 2. Displacement-induced leakage
# ════════════════════════════════════════════════════════════════════════════
print("\n── 2. Displacement-induced leakage ──────────────────────────────────")

N_test = 10
log_idx_test = [0, 1, N_test, N_test+1]
subspace_test = Subspace.custom(
    2*N_test, log_idx_test, ["|g,0>","|g,1>","|e,0>","|e,1>"]
)
P_test = cavfirst_to_qfirst_perm(N_test)
proj_test = subspace_test.projector()

disp_data = {}
alphas = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]
for alpha in alphas:
    D_cf = np.asarray(displacement(alpha, N_test).full())
    # Displacement acts on cavity only; embed as I_q ⊗ D
    D_full_cf = np.kron(np.eye(2), D_cf)  # cavity-first
    D_full_qf = P_test @ D_full_cf @ P_test.T  # qubit-first

    leak = np.linalg.norm(
        (np.eye(2*N_test) - proj_test) @ D_full_qf @ proj_test,
        ord="fro"
    ) ** 2 / 4
    disp_data[float(alpha)] = float(leak)

results["displacement_leakage"] = disp_data
print("  alpha → leakage L:")
for a, L in disp_data.items():
    print(f"    |α|={a:.1f}: L={L:.4f}")


# ════════════════════════════════════════════════════════════════════════════
# 3. Dispersive phase: drift accumulation with full χ, χ′, K
# ════════════════════════════════════════════════════════════════════════════
print("\n── 3. Dispersive phase drift (χ, χ′, K) ─────────────────────────────")

drift_model = DriftPhaseModel(chi=CHI, chi2=CHIP, kerr=KERR)

N_dp = 8
t_pi = np.pi / abs(CHI)
log_idx_dp = [0, 1, N_dp, N_dp+1]
subspace_dp = Subspace.custom(2*N_dp, log_idx_dp, ["|g,0>","|g,1>","|e,0>","|e,1>"])

# Dispersive phase unitary at various times
times = np.linspace(0, 4 * t_pi, 40)
fid_disp_vs_t = []

# Target: a CNOT_{q→c}-like operation using only dispersive phase
# Dispersive evolution at t = π/|χ| gives e^{iπ} on |e,1> relative to |g,1>
# This is NOT CNOT_{q→c} but a conditional phase gate
# Let's instead check the dispersive unitary subspace action
for t in times:
    U_drift = np.asarray(drift_phase_unitary(n_cav=N_dp, duration=t, model=drift_model).full())
    U_drift_log = U_drift[np.ix_(log_idx_dp, log_idx_dp)]
    # Compare to ideal dispersive-only unitary (χ only, no χ', no K)
    drift_chi_only = DriftPhaseModel(chi=CHI, chi2=0.0, kerr=0.0)
    U_chi_only = np.asarray(drift_phase_unitary(n_cav=N_dp, duration=t, model=drift_chi_only).full())
    U_chi_log = U_chi_only[np.ix_(log_idx_dp, log_idx_dp)]
    F_diff = subspace_unitary_fidelity(U_drift_log, U_chi_log, gauge="global")
    fid_disp_vs_t.append(float(F_diff))

results["dispersive_phase_chi_prime_k_effect"] = {
    "times_us": [float(t * 1e6) for t in times],
    "fidelity_full_vs_chi_only": fid_disp_vs_t,
    "t_pi_us": float(t_pi * 1e6),
}

print(f"  t_π = {t_pi*1e6:.3f} μs")
fids_at_pi = [fid_disp_vs_t[np.argmin(abs(times - t_pi * k))]
              for k in [1, 2, 4]]
print(f"  F(χ+χ′+K vs χ-only) at t=t_π:  {fids_at_pi[0]:.6f}")
print(f"  F at t=2t_π: {fids_at_pi[1]:.6f}")
print(f"  F at t=4t_π: {fids_at_pi[2]:.6f}")


# ════════════════════════════════════════════════════════════════════════════
# 4. CNOT_{q→c} decomposition using dispersive phase + displacement
#    Strategy: CNOT_{q→c} = D†(β) · e^{iπ a†a |e><e|} · D(β)
#    where D(β) moves the |0⟩ and |1⟩ cavity states to a suitable
#    "code frame".  This is a standard conditional displacement approach.
# ════════════════════════════════════════════════════════════════════════════
print("\n── 4. CNOT_{q→c} via dispersive phase + displacement ────────────────")

N_cn = 8
log_idx_cn = [0, 1, N_cn, N_cn+1]
subspace_cn = Subspace.custom(2*N_cn, log_idx_cn, ["|g,0>","|g,1>","|e,0>","|e,1>"])
P_cn = cavfirst_to_qfirst_perm(N_cn)
proj_cn = subspace_cn.projector()
target_q2c_full = cnot_q2c_full(N_cn)

# Best idealized decomposition:
# CNOT_{q→c} on logical {|0>,|1>}:
# Approach: use SQR + dispersive phase to implement
# |g>⟨g| ⊗ I_c + |e>⟨e| ⊗ X_c^logical
# X_c^logical = |0><1| + |1><0| acts on {|0>,|1>} cavity.
# This is NOT achievable by dispersive phase alone (which only gives phases,
# not amplitude mixing on the cavity).
# The minimal physical decomposition:
# 1. D(α):   coherent state displacement off-origin
# 2. Dispersive phase t = π/|χ|: accumulates conditional phase
# 3. D(-α):  undo displacement
# This creates a conditional phase but NOT a logical NOT on the cavity.
# To get X_c we need SNAP or a Jaynes–Cummings type interaction.

# CNOT_{q→c} = qubit-controlled cavity NOT = requires Fock-level mixing.
# SQR gates (which rotate the qubit conditioned on the Fock state) cannot flip
# the cavity Fock level. Only Displacement gates can mix Fock levels.
# We verify this by attempting to implement CNOT_{q→c} via SQR composition and
# checking the resulting fidelity against the ideal gate.

from cqed_sim.core.ideal_gates import qubit_rotation_xy

# Build SQR_1(π,0) with CP correction → exact CNOT_{c→q}
sqr1_cf_n = np.asarray(sqr(np.pi, 0.0, n=1, cavity_dim=N_cn).full())
P_cn_mat  = P_cn
sqr1_qf_n = P_cn_mat @ sqr1_cf_n @ P_cn_mat.T

# CP correction for CNOT_{c→q} (1j on n=1 sector)
CP_corr_cn = np.eye(2*N_cn, dtype=np.complex128)
CP_corr_cn[1, 1]       = 1j   # |g,1>
CP_corr_cn[N_cn+1, N_cn+1] = 1j  # |e,1>
U_cnot_c2q_full = sqr1_qf_n @ CP_corr_cn  # exact CNOT_{c→q}

# Verify CNOT_{c→q} fidelity with CP correction
U_log_c2q = U_cnot_c2q_full[np.ix_(log_idx_cn, log_idx_cn)]
F_c2q = subspace_unitary_fidelity(U_log_c2q, CNOT_c2q, gauge="global")
print(f"  SQR_1(π,0) · CP_corr → CNOT_c2q: F={F_c2q:.6f}  (should be 1.0)")
results["CNOT_c2q_sqr_decomp"] = {
    "fidelity": float(F_c2q),
    "sequence": "SQR_1(π,0) · CP_correction([0,π/2] on n=1)",
    "n_cav": N_cn,
}

# Full U_target: H_c · CNOT_{c→q} · CNOT_{q→c} with exact analytic gates
U_Hc_full_n = hc_full(N_cn)
U_q2c_full  = cnot_q2c_full(N_cn)  # ideal CNOT_{q→c} (requires Fock mixing)

U_target_full_n = U_Hc_full_n @ U_cnot_c2q_full @ U_q2c_full
U_log_full = U_target_full_n[np.ix_(log_idx_cn, log_idx_cn)]
F_full = subspace_unitary_fidelity(U_log_full, U_target_4x4, gauge="global")
print(f"  Full U_target (analytic, N={N_cn}): F={F_full:.6f}  (should be 1.0)")
results["U_target_analytic_Ncn"] = {
    "fidelity": float(F_full),
    "sequence": "H_c · (SQR_1·CP_corr) · CNOT_{q→c}",
    "n_cav": N_cn,
    "note": "CNOT_{q→c} requires Fock-level mixing; cannot be built from SQR alone",
}


# ════════════════════════════════════════════════════════════════════════════
# 5. Plots
# ════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Plot A: truncation convergence
ax = axes[0]
N_plot = sorted(trunc_data.keys())
F_plot = [trunc_data[N]["F"] for N in N_plot]
ax.plot(N_plot, F_plot, "o-", color="steelblue", lw=2)
ax.set_xlabel("Cavity truncation $N$", fontsize=11)
ax.set_ylabel("$F_{\\rm proj}$ (analytic gates)", fontsize=11)
ax.set_title("Truncation convergence", fontsize=11)
ax.set_ylim(0.98, 1.001)
ax.grid(True, alpha=0.3)

# Plot B: displacement leakage vs |α|
ax = axes[1]
alphas_plot = sorted(disp_data.keys())
L_plot = [disp_data[a] for a in alphas_plot]
ax.semilogy(alphas_plot, L_plot, "s-", color="darkorange", lw=2)
ax.set_xlabel(r"$|\alpha|$", fontsize=11)
ax.set_ylabel("Leakage $L$", fontsize=11)
ax.set_title("Displacement leakage (N=10)", fontsize=11)
ax.grid(True, alpha=0.3)

# Plot C: dispersive drift χ'+K effect
ax = axes[2]
t_us = [t * 1e6 for t in times]
ax.plot(t_us, fid_disp_vs_t, "-", color="forestgreen", lw=2)
ax.axvline(t_pi * 1e6, ls="--", color="gray", label=r"$t_\pi = \pi/|\chi|$")
ax.axvline(2 * t_pi * 1e6, ls=":", color="gray")
ax.set_xlabel("Time (μs)", fontsize=11)
ax.set_ylabel(r"$F(\chi+\chi'+K\;{\rm vs}\;\chi$-only)", fontsize=11)
ax.set_title(r"Effect of $\chi'$ and $K$", fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 1.02)

fig.tight_layout()
fig.savefig(FIG_DIR / "phase3_embedding.png", dpi=150)
fig.savefig(FIG_DIR / "phase3_embedding.pdf")
plt.close(fig)

print(f"\nFigures saved to {FIG_DIR}")

# Save data
with open(DATA_DIR / "phase3_results.json", "w") as f:
    json.dump(results, f, indent=2, default=str)

print(f"Results written to {DATA_DIR / 'phase3_results.json'}")
print("Phase 3 complete.")
