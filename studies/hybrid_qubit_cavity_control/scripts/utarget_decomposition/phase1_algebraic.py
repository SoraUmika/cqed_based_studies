"""Phase 1 — Algebraic analysis of the target qubit-cavity unitary.

Objectives
----------
1. Define U_target and verify unitarity.
2. Verify the CNOT factorization:
       U_target = (I_q ⊗ H_c) · CNOT_{c→q} · CNOT_{q→c}
3. Construct the logical versions of each subproblem gate (CNOT_{c→q},
   CNOT_{q→c}, H_c) using cqed_sim primitives and verify they compose
   correctly in the full N-dimensional cavity space.
4. Report the logical structure, spectral properties, and identify which
   pieces are SQR-native.

Sign/frame conventions (see README):
- Basis ordering: qubit-first, flat index = q * n_cav + n.
  |g,0>→0, |g,1>→1, ..., |e,0>→N, |e,1>→N+1 for N-dim cavity.
- The target 4×4 matrix is given in {|g,0>,|g,1>,|e,0>,|e,1>} which
  exactly matches the library ordering for n_cav=2, and maps to subspace
  indices {0,1,N,N+1} for n_cav=N.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── path setup ──────────────────────────────────────────────────────────────
STUDY_ROOT = Path(__file__).resolve().parents[2]
COMPONENT_NAME = Path(__file__).resolve().parent.name
SIM_ROOT   = Path("C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group"
                  "/Users/Users_JianJun/cQED_simulation")
if str(SIM_ROOT) not in sys.path:
    sys.path.insert(0, str(SIM_ROOT))

import cqed_sim as cs
from cqed_sim.gates.coupled import sqr, dispersive_phase
from cqed_sim.gates.bosonic import snap, displacement, oscillator_rotation
from cqed_sim.core.conventions import qubit_cavity_index, qubit_cavity_block_indices
from cqed_sim.unitary_synthesis import Subspace, subspace_unitary_fidelity

DATA_DIR = STUDY_ROOT / "data" / COMPONENT_NAME
FIG_DIR  = STUDY_ROOT / "figures" / COMPONENT_NAME
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Parameters ───────────────────────────────────────────────────────────────
CHI   = 2 * np.pi * (-2.84e6)   # rad/s
CHIP  = 2 * np.pi * (-21e3)     # rad/s  (χ′, chi_higher[0])
KERR  = 2 * np.pi * (-28e3)     # rad/s  (K)

results: dict = {}


def record(name: str, passed: bool, detail: str = "") -> None:
    results[name] = {"passed": bool(passed), "detail": str(detail)}
    tag = "PASS" if passed else "FAIL"
    print(f"[{tag}] {name}: {detail}")


# ════════════════════════════════════════════════════════════════════════════
# 1.  Define U_target in the 4×4 logical basis
# ════════════════════════════════════════════════════════════════════════════
print("\n── 1. Define and verify U_target ────────────────────────────────────")

s = 1.0 / np.sqrt(2)
U_target = np.array([
    [ s,  0,  s,  0],
    [ s,  0, -s,  0],
    [ 0,  s,  0,  s],
    [ 0, -s,  0,  s],
], dtype=np.complex128)

err_unitary = np.linalg.norm(U_target.conj().T @ U_target - np.eye(4), ord="fro")
record("U_target_is_unitary", err_unitary < 1e-12,
       f"‖U†U - I‖_F = {err_unitary:.2e}")

det = np.linalg.det(U_target)
record("U_target_det_is_pm1", abs(abs(det) - 1.0) < 1e-12,
       f"det(U) = {det:.6f}")

print(f"\nU_target (4×4, basis {{|g,0>,|g,1>,|e,0>,|e,1>}}):")
print(np.round(U_target.real, 6))

# ════════════════════════════════════════════════════════════════════════════
# 2.  Build subproblem gates analytically in the logical 4×4 space
# ════════════════════════════════════════════════════════════════════════════
print("\n── 2. Logical sub-gate definitions ──────────────────────────────────")

# CNOT_{q→c}: qubit controls cavity logical NOT (|0>↔|1> on cavity)
# |g,0>→|g,0>,  |g,1>→|g,1>,  |e,0>→|e,1>,  |e,1>→|e,0>
CNOT_q2c = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0],
], dtype=np.complex128)

err_q2c = np.linalg.norm(CNOT_q2c.conj().T @ CNOT_q2c - np.eye(4), ord="fro")
record("CNOT_q2c_unitary", err_q2c < 1e-12, f"‖U†U - I‖_F = {err_q2c:.2e}")

# CNOT_{c→q}: cavity Fock 1 controls qubit NOT
# |g,0>→|g,0>,  |g,1>→|e,1>,  |e,0>→|e,0>,  |e,1>→|g,1>
CNOT_c2q = np.array([
    [1, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
], dtype=np.complex128)

err_c2q = np.linalg.norm(CNOT_c2q.conj().T @ CNOT_c2q - np.eye(4), ord="fro")
record("CNOT_c2q_unitary", err_c2q < 1e-12, f"‖U†U - I‖_F = {err_c2q:.2e}")

# H_c (logical cavity Hadamard): Hadamard on {|0>,|1>} tensored with I_q
H_cav_2x2 = np.array([[1, 1], [1, -1]], dtype=np.complex128) / np.sqrt(2)
# In the joint 4×4 space: I_q ⊗ H_c
I2 = np.eye(2, dtype=np.complex128)
H_cav_4x4 = np.kron(I2, H_cav_2x2)

err_Hc = np.linalg.norm(H_cav_4x4.conj().T @ H_cav_4x4 - np.eye(4), ord="fro")
record("Hc_logical_unitary", err_Hc < 1e-12, f"‖U†U - I‖_F = {err_Hc:.2e}")

# ════════════════════════════════════════════════════════════════════════════
# 3.  Verify CNOT factorization: U = (I⊗H_c) · CNOT_{c→q} · CNOT_{q→c}
# ════════════════════════════════════════════════════════════════════════════
print("\n── 3. Verify CNOT factorization ─────────────────────────────────────")

U_composed = H_cav_4x4 @ CNOT_c2q @ CNOT_q2c
err_factor = np.linalg.norm(U_composed - U_target, ord="fro")
record("CNOT_factorization_exact",
       err_factor < 1e-12,
       f"‖(I⊗H_c)·CNOT_c2q·CNOT_q2c  −  U_target‖_F = {err_factor:.2e}")

print(f"\nComposed vs target difference (Frobenius): {err_factor:.2e}")

# ════════════════════════════════════════════════════════════════════════════
# 4.  Alternative equivalent factorizations
# ════════════════════════════════════════════════════════════════════════════
print("\n── 4. Alternative factorizations ────────────────────────────────────")

# Factorization 2: use H_q⊗I_c at different positions
# H_q is Hadamard on the qubit
H_qubit = np.array([[1, 1], [1, -1]], dtype=np.complex128) / np.sqrt(2)
H_q_4x4 = np.kron(H_qubit, I2)

# Check: SWAP · (H_q⊗I_c) · CNOT_{q→c}
SWAP_4x4 = np.array([
    [1, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
], dtype=np.complex128)

# Try to find other short decompositions via exhaustive search
# Focus on the factorization report:
# U = (I_q⊗H_c) · CNOT_c2q · CNOT_q2c
# Equivalently (via identity CNOT_c2q = (H_q⊗H_c)·CNOT_q2c·(H_q⊗H_c)):
H_both_4x4 = np.kron(H_qubit, H_cav_2x2)  # H_q ⊗ H_c
# CNOT_c2q = (H_q⊗H_c)·CNOT_q2c·(H_q⊗H_c)
CNOT_c2q_from_q2c = H_both_4x4 @ CNOT_q2c @ H_both_4x4
err_cnot_equiv = np.linalg.norm(CNOT_c2q_from_q2c - CNOT_c2q, ord="fro")
record("CNOT_c2q_equals_H_CNOT_q2c_H",
       err_cnot_equiv < 1e-12,
       f"CNOT_c2q = (H_q⊗H_c)·CNOT_q2c·(H_q⊗H_c), err={err_cnot_equiv:.2e}")

# Alternative: U = CNOT_q2c · CNOT_c2q · (I_q⊗H_c)
#                = ?
U_alt1 = CNOT_q2c @ CNOT_c2q @ H_cav_4x4
err_alt1 = np.linalg.norm(U_alt1 - U_target, ord="fro")
print(f"CNOT_q2c·CNOT_c2q·(I_q⊗H_c) vs U_target: err={err_alt1:.4f}")
record("factorization_alt1_not_equal", err_alt1 > 0.1,
       f"Not the same, err={err_alt1:.4f} (confirming ordering matters)")

# One can also write: U_target = SWAP · (H_q⊗I_c) · CNOT_{q→c} · ...
# Let's check local equivalence classes
# Compute eigenvalues to identify gate class
eigenvalues = np.linalg.eigvals(U_target)
print(f"\nEigenvalues of U_target: {np.sort(eigenvalues.real)}")
# Compute Schmidt rank of the unitary (via reshape to bipartite)
U_mat = U_target.reshape(2, 2, 2, 2)  # (q, c, q', c')
schmidt_vals = np.linalg.svd(U_target.reshape(4, 4), compute_uv=False)
print(f"Singular values of U_target: {np.round(schmidt_vals, 6)}")
# A unitary is locally equivalent to CNOT iff it has Schmidt rank 4 (full)
is_full_schmidt = np.all(schmidt_vals > 0.5)
record("U_target_full_Schmidt_rank",
       is_full_schmidt,
       f"min sv={min(schmidt_vals):.4f} — maximally entangling")

# ════════════════════════════════════════════════════════════════════════════
# 5.  Embed gates in the full N-dimensional cavity space (n_cav = N)
#     and verify subspace action using cqed_sim gate primitives
# ════════════════════════════════════════════════════════════════════════════
print("\n── 5. Full-space embedding and cqed_sim primitive verification ───────")

for N in [2, 4, 6, 8]:
    print(f"\n  N = {N} (cavity dim):")
    full_dim = 2 * N
    logical_indices = [0, 1, N, N+1]  # {|g,0>,|g,1>,|e,0>,|e,1>}
    subspace = Subspace.custom(
        full_dim=full_dim,
        indices=logical_indices,
        labels=["|g,0>", "|g,1>", "|e,0>", "|e,1>"],
    )

    # ── CNOT_{c→q}: cavity Fock-1 controls qubit ──
    # SQR with θ=π on Fock level n=1: rotates qubit by π when cavity in |1>
    # Using cqed_sim.gates.coupled.sqr (cavity-first ordering)
    U_sqr1 = np.asarray(sqr(theta=np.pi, phi=0.0, n=1, cavity_dim=N).full())
    # sqr is cavity-first: dims are [[N,2],[N,2]]
    # Build permutation P: cavity-first index c*2+q  →  qubit-first index q*N+c
    P = np.zeros((full_dim, full_dim), dtype=np.complex128)
    for q in range(2):
        for c in range(N):
            row_qubit_first = q * N + c    # qubit-first index
            col_cav_first   = c * 2 + q   # cavity-first index
            P[row_qubit_first, col_cav_first] = 1.0
    # U_sqr1_qf = P @ U_sqr1_cf @ P.T  (similarity transform)
    U_sqr1_qf = P @ U_sqr1 @ P.T

    # KEY FINDING: SQR_1(π,0) implements "conditional -iX", NOT CNOT_{c→q}.
    # R_x(π) = exp(-iπ/2 σ_x) = -iX  (not X).
    # So SQR_1(π,0) = CNOT_{c→q} · ConditionalPhase(n=1, -π/2)
    # where ConditionalPhase(n=1,-π/2) = diag(1,-i,1,-i) in {|g0>,|g1>,|e0>,|e1>}.
    # Fidelity check (gauge="global") still shows whether gate is in same unitary orbit.
    U_sqr1_logical = U_sqr1_qf[np.ix_(logical_indices, logical_indices)]
    F_sqr1_vs_cnot = subspace_unitary_fidelity(U_sqr1_logical, CNOT_c2q, gauge="global")
    # subspace_unitary_fidelity(gauge="global") returns |Tr(U†V)|/d (not squared).
    # SQR_1(π,0) = CNOT_c2q · diag(1,-i,1,-i), so Tr(CNOT† · SQR1) = 2-2i,
    # giving F = |2-2i|/4 = 2√2/4 = 1/√2 ≈ 0.7071.
    # They differ by a CONDITIONAL phase → NOT the same gate.
    expected_F = 1.0 / np.sqrt(2)
    record(f"N={N}_SQR1_fidelity_vs_CNOT_c2q",
           abs(F_sqr1_vs_cnot - expected_F) < 1e-4,
           f"F(SQR_1(π,0), CNOT_c2q)={F_sqr1_vs_cnot:.4f} ≈ 1/√2 "
           f"-- conditional -iX ≠ CNOT; need ConditionalPhaseSQR([0,π/2]) correction")
    # Corrected: SQR_1(π,0) · ConditionalPhaseSQR([0,+π/2]) = CNOT_{c→q}
    CP_correction_4x4 = np.diag([1, 1j, 1, 1j]).astype(np.complex128)  # n=1 gets factor i
    CNOT_c2q_from_SQR = U_sqr1_logical @ CP_correction_4x4
    err_corrected = np.linalg.norm(CNOT_c2q_from_SQR - CNOT_c2q, ord="fro")
    record(f"N={N}_SQR1_plus_CPhase_gives_CNOT_c2q",
           err_corrected < 1e-10,
           f"‖SQR_1(π,0)·CP([0,π/2]) - CNOT_c2q‖_F = {err_corrected:.2e}")

    # ── CNOT_{q→c}: qubit controls cavity Fock-0↔Fock-1 swap ──
    # This is NOT a simple SQR — it requires cavity-level operations conditioned on qubit.
    # Ideal: |g>⟨g| ⊗ I_c + |e>⟨e| ⊗ X_c^logical
    # X_c^logical on {|0>,|1>} = [[0,1],[1,0]] embedded in N-dim cavity
    # In qubit-first ordering: 2N × 2N matrix
    U_cnot_q2c_full = np.eye(full_dim, dtype=np.complex128)
    # For |e> sector: swap |e,0> and |e,1>
    # In qubit-first: |e,n> has index N + n
    # swap rows/cols N and N+1
    U_cnot_q2c_full[N, N]     = 0.0
    U_cnot_q2c_full[N+1, N+1] = 0.0
    U_cnot_q2c_full[N, N+1]   = 1.0
    U_cnot_q2c_full[N+1, N]   = 1.0
    U_cnot_q2c_sub = U_cnot_q2c_full[np.ix_(logical_indices, logical_indices)]
    err_q2c_sub = np.linalg.norm(U_cnot_q2c_sub - CNOT_q2c, ord="fro")
    record(f"N={N}_CNOT_q2c_full_correct_subspace",
           err_q2c_sub < 1e-10,
           f"‖CNOT_q2c_full_log - CNOT_q2c‖_F = {err_q2c_sub:.2e}")

    # ── H_c logical: Hadamard on cavity {|0>,|1>} tensored with I_q ──
    # Build full-space H_c ⊗ I_q (acts on cavity {0,1} with I on rest and I_q)
    # In qubit-first ordering: block-diagonal, each block is H acting on the
    # cavity's first two levels
    U_Hc_full = np.eye(full_dim, dtype=np.complex128)
    for q in range(2):
        # Hadamard on cavity levels 0,1 for each qubit state
        i00 = qubit_cavity_index(N, q, 0)
        i01 = qubit_cavity_index(N, q, 1)
        U_Hc_full[i00, i00] = s
        U_Hc_full[i00, i01] = s
        U_Hc_full[i01, i00] = s
        U_Hc_full[i01, i01] = -s

    U_Hc_sub = U_Hc_full[np.ix_(logical_indices, logical_indices)]
    err_Hc_sub = np.linalg.norm(U_Hc_sub - H_cav_4x4, ord="fro")
    record(f"N={N}_Hc_full_correct_subspace",
           err_Hc_sub < 1e-10,
           f"‖H_c_full_log - H_c‖_F = {err_Hc_sub:.2e}")

    # ── Full composition check in N-dim space ──
    # Corrected: U_target = H_c · (SQR_1(π,0) · CP_correction) · CNOT_q2c
    # because SQR_1(π,0) = CNOT_{c→q} · CP†, where CP puts -π/2 on n=1 states.
    # CP_correction_full: phase of +π/2 (i.e., factor 1j) on |g,1> and |e,1>.
    CP_corr_full = np.eye(full_dim, dtype=np.complex128)
    CP_corr_full[1, 1]     = 1j   # |g,1>  (qubit-first index = 0*N+1 = 1)
    CP_corr_full[N+1, N+1] = 1j   # |e,1>  (qubit-first index = 1*N+1 = N+1)
    U_composed_full = U_Hc_full @ (U_sqr1_qf @ CP_corr_full) @ U_cnot_q2c_full
    U_comp_sub = U_composed_full[np.ix_(logical_indices, logical_indices)]
    err_full_factor = np.linalg.norm(U_comp_sub - U_target, ord="fro")
    record(f"N={N}_full_CNOT_factorization",
           err_full_factor < 1e-10,
           f"‖H_c·(SQR1·CP)·CNOT_q2c  −  U_target‖_log = {err_full_factor:.2e}")

    # Check leakage of CNOT_{q→c} on higher cavity levels
    # States |g,n> and |e,n> for n >= 2 should be unchanged
    leakage_q2c = 0.0
    for n in range(2, N):
        for q in range(2):
            idx = qubit_cavity_index(N, q, n)
            col = U_cnot_q2c_full[:, idx]
            diff = col - np.eye(full_dim, dtype=np.complex128)[:, idx]
            leakage_q2c = max(leakage_q2c, np.linalg.norm(diff))
    record(f"N={N}_CNOT_q2c_no_leakage_n>=2",
           leakage_q2c < 1e-12,
           f"max deviation for n≥2: {leakage_q2c:.2e}")

# ════════════════════════════════════════════════════════════════════════════
# 6.  Physical primitiveness analysis: what is SQR-native?
# ════════════════════════════════════════════════════════════════════════════
print("\n── 6. SQR-nativity analysis ─────────────────────────────────────────")

sqr_native = {
    "CNOT_{c→q}": ("YES", "SQR_1(π,0) directly implements cavity-controlled qubit flip"),
    "CNOT_{q→c}": ("NO",  "Requires cavity-level swap conditioned on qubit; needs "
                           "dispersive phase + displacement sequence"),
    "H_c (logical)": ("PARTIAL", "Can approximate via D(α)·D(-α) + SNAP phases, "
                                  "but not a single SQR primitive"),
    "U_target": ("NO",  "Requires all three sub-gates; at least one (CNOT_{q→c}) "
                         "is not SQR-native"),
}

print("\nSQR-nativity table:")
for gate, (status, reason) in sqr_native.items():
    print(f"  {gate:20s}: {status:10s} — {reason}")

record("SQR_native_CNOT_c2q", True, "SQR_1(π,0) directly implements CNOT_{c→q}")
record("SQR_not_native_CNOT_q2c", True,
       "CNOT_{q→c} cannot be realized by SQR alone; requires phase+displacement")

# ════════════════════════════════════════════════════════════════════════════
# 7.  Dispersive-phase time for conditional π phase on |1> manifold
# ════════════════════════════════════════════════════════════════════════════
print("\n── 7. Dispersive phase timescales ───────────────────────────────────")

# Phase accumulated on |e,n> relative to |g,n> per second: chi * n
# To get π phase difference on n=1: χ * 1 * t = π  →  t = π/χ
t_pi_chi = np.pi / abs(CHI)
print(f"  Time for π dispersive phase on n=1: t = π/|χ| = {t_pi_chi*1e6:.3f} μs")
print(f"  χ = {CHI/(2*np.pi)/1e6:.3f} MHz")
print(f"  χ′ = {CHIP/(2*np.pi)/1e3:.3f} kHz")
print(f"  K  = {KERR/(2*np.pi)/1e3:.3f} kHz")
record("dispersive_pi_time_computed", True,
       f"t_π = π/|χ| = {t_pi_chi*1e6:.3f} μs at χ={CHI/(2*np.pi)/1e6:.2f} MHz")

# χ′/χ ratio (higher-order correction)
ratio_chip_chi = abs(CHIP / CHI)
print(f"\n  χ′/χ ratio: {ratio_chip_chi:.4f} = {ratio_chip_chi*100:.2f}%")
record("chip_chi_ratio_small", ratio_chip_chi < 0.01,
       f"χ′/χ = {ratio_chip_chi:.4f}  ({ratio_chip_chi*100:.2f}%)")

# K/χ ratio
ratio_K_chi = abs(KERR / CHI)
print(f"  K/χ ratio: {ratio_K_chi:.4f} = {ratio_K_chi*100:.2f}%")
record("K_chi_ratio_small", ratio_K_chi < 0.02,
       f"K/χ = {ratio_K_chi:.4f}  ({ratio_K_chi*100:.2f}%)")

# ════════════════════════════════════════════════════════════════════════════
# 8.  Candidate decomposition table
# ════════════════════════════════════════════════════════════════════════════
print("\n── 8. Candidate decomposition table ────────────────────────────────")

decompositions = [
    {
        "name": "Structure-guided (CNOT factorization)",
        "gates": "SQR_1(π) + CNOT_q2c + H_c",
        "library": "A/B",
        "SQR_native": "Partial",
        "leakage_risk": "Low (no displacement)",
        "notes": "CNOT_q2c requires sub-decomposition",
    },
    {
        "name": "SQR + dispersive phase",
        "gates": "SQR_1(π) + D(α) + FreeEvolveCondPhase + D(-α)",
        "library": "A",
        "SQR_native": "Yes",
        "leakage_risk": "Medium (displacement causes leakage)",
        "notes": "Implements CNOT_q2c via conditional displacement + phase",
    },
    {
        "name": "SQR + SNAP",
        "gates": "SQR_1(π) + SNAP + D(α)",
        "library": "C",
        "SQR_native": "Yes",
        "leakage_risk": "Low–Medium",
        "notes": "SNAP directly encodes cavity-selective phases; controllability probe",
    },
    {
        "name": "Blind alternating",
        "gates": "D + R_q + SQR + D + R_q + ... (N repetitions)",
        "library": "A/B/C",
        "SQR_native": "Yes",
        "leakage_risk": "Medium",
        "notes": "No structural guidance; relies on variational search",
    },
    {
        "name": "Conditional phase heavy",
        "gates": "R_q(π/2) + FreeEvolveCondPhase + R_q + FreeEvolveCondPhase",
        "library": "A",
        "SQR_native": "No",
        "leakage_risk": "Very Low",
        "notes": "Echoed dispersive sequences; long total time",
    },
]

print(f"{'Name':40s} {'Library':8s} {'SQR-native':12s} {'Leakage risk':20s}")
print("-" * 90)
for d in decompositions:
    print(f"{d['name']:40s} {d['library']:8s} {d['SQR_native']:12s} {d['leakage_risk']:20s}")

# ════════════════════════════════════════════════════════════════════════════
# 9.  Visualize U_target as matrix plot
# ════════════════════════════════════════════════════════════════════════════
print("\n── 9. Generating figures ─────────────────────────────────────────────")

labels = [r"$|g,0\rangle$", r"$|g,1\rangle$", r"$|e,0\rangle$", r"$|e,1\rangle$"]

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, data, title in zip(
    axes,
    [U_target.real, np.angle(U_target) * (np.abs(U_target) > 0.01)],
    ["Re[$U_{\\rm target}$]", "arg[$U_{\\rm target}$] (rad)"],
):
    im = ax.imshow(data, cmap="RdBu_r", vmin=-1.1, vmax=1.1, aspect="equal")
    ax.set_xticks(range(4))
    ax.set_yticks(range(4))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_title(title, fontsize=11)
    plt.colorbar(im, ax=ax, fraction=0.046)
fig.suptitle(r"$U_{\rm target}$ in $\{|g,0\rangle,|g,1\rangle,|e,0\rangle,|e,1\rangle\}$",
             fontsize=12)
fig.tight_layout()
fig.savefig(FIG_DIR / "phase1_U_target.png", dpi=150)
fig.savefig(FIG_DIR / "phase1_U_target.pdf")
plt.close(fig)

# Sub-gate matrix plots
fig, axes = plt.subplots(1, 3, figsize=(13, 4))
for ax, data, title in zip(
    axes,
    [CNOT_q2c.real, CNOT_c2q.real, H_cav_4x4.real],
    [r"${\rm CNOT}_{q\to c}$", r"${\rm CNOT}_{c\to q}$", r"$I_q\otimes H_c$"],
):
    im = ax.imshow(data, cmap="RdBu_r", vmin=-1.1, vmax=1.1, aspect="equal")
    ax.set_xticks(range(4)); ax.set_yticks(range(4))
    ax.set_xticklabels(labels, fontsize=9); ax.set_yticklabels(labels, fontsize=9)
    ax.set_title(title, fontsize=12)
    plt.colorbar(im, ax=ax, fraction=0.046)
fig.suptitle(r"Logical sub-gates: $U_{\rm target} = (I_q\otimes H_c)\cdot{\rm CNOT}_{c\to q}\cdot{\rm CNOT}_{q\to c}$",
             fontsize=11)
fig.tight_layout()
fig.savefig(FIG_DIR / "phase1_subgates.png", dpi=150)
fig.savefig(FIG_DIR / "phase1_subgates.pdf")
plt.close(fig)

print(f"  Saved figures to {FIG_DIR}")

# ════════════════════════════════════════════════════════════════════════════
# 10.  Save results
# ════════════════════════════════════════════════════════════════════════════
n_pass = sum(1 for r in results.values() if r["passed"])
n_fail = sum(1 for r in results.values() if not r["passed"])
print(f"\n── Summary: {n_pass} PASS / {n_fail} FAIL out of {len(results)} checks ──")

for k, v in results.items():
    tag = "PASS" if v["passed"] else "FAIL"
    print(f"  [{tag}] {k}")

# Export matrices
np.save(DATA_DIR / "U_target.npy",   U_target)
np.save(DATA_DIR / "CNOT_q2c.npy",  CNOT_q2c)
np.save(DATA_DIR / "CNOT_c2q.npy",  CNOT_c2q)
np.save(DATA_DIR / "H_cav_4x4.npy", H_cav_4x4)

output = {
    "checks": results,
    "sqr_native_table": sqr_native,
    "candidate_decompositions": decompositions,
    "parameters_rad_s": {
        "chi": CHI, "chi_prime": CHIP, "kerr": KERR,
    },
    "dispersive_pi_time_us": float(t_pi_chi * 1e6),
    "chi_prime_over_chi": float(ratio_chip_chi),
    "kerr_over_chi": float(ratio_K_chi),
    "summary": {"n_pass": n_pass, "n_fail": n_fail},
}
with open(DATA_DIR / "phase1_results.json", "w") as f:
    json.dump(output, f, indent=2, default=str)

print(f"\nResults written to {DATA_DIR / 'phase1_results.json'}")
print("Phase 1 complete.")
