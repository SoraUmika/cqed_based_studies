"""Cluster-state holographic simulation study — main execution script.

Phases:
  1. Verify the cqed_sim cluster target against the canonical definition.
  2. Compute ideal cluster-state observables (Pauli, stabilizers, ZXZ, string-order).
  3. Decompose the per-site unitary into QubitRotation + Displacement + SNAP.
  4. GRAPE-optimize each SNAP gate for minimum pulse duration.
  5. Timing analysis (SNAP-only + full sequence).
  6. Waveform-level simulation of the decomposed sequence.
"""

from __future__ import annotations

import json
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Any

os.environ["PYTHONIOENCODING"] = "utf-8"
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import qutip as qt

# ── Paths ────────────────────────────────────────────────────────────────────
STUDY_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = STUDY_ROOT.parents[1]
DATA_DIR = STUDY_ROOT / "data"
FIG_DIR = STUDY_ROOT / "figures"
SIM_ROOT = Path(
    "C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group"
    "/Users/Users_JianJun/cQED_simulation"
)
STYLE_PATH = WORKSPACE_ROOT / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"

DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

if str(SIM_ROOT) not in sys.path:
    sys.path.insert(0, str(SIM_ROOT))
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))

TWO_PI = 2.0 * np.pi
TOL_BRIGHT = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']

# ── cqed_sim imports ─────────────────────────────────────────────────────────
import cqed_sim  # noqa: E402
from cqed_sim import (  # noqa: E402
    DispersiveTransmonCavityModel,
    FrameSpec,
    ModelControlChannelSpec,
    PiecewiseConstantTimeGrid,
    SequenceCompiler,
    SimulationConfig,
    simulate_sequence,
)
from cqed_sim.core.ideal_gates import snap_op, embed_qubit_op, embed_cavity_op  # noqa: E402
from cqed_sim.gates.bosonic import snap, displacement  # noqa: E402
from cqed_sim.gates.qubit import rotation_xy  # noqa: E402
from cqed_sim.optimal_control import (  # noqa: E402
    GrapeConfig,
    GrapeSolver,
    LeakagePenalty as OCLeakagePenalty,
    UnitaryObjective as OCUnitaryObjective,
    build_control_problem_from_model,
)
from cqed_sim.unitary_synthesis import (  # noqa: E402
    Displacement as SynthDisplacement,
    GateSequence,
    QubitRotation as SynthQubitRotation,
    SNAP as SynthSNAP,
    Subspace,
    TargetUnitary,
    UnitarySynthesizer,
    subspace_unitary_fidelity,
)
from cqed_sim.unitary_synthesis.targets import make_target  # noqa: E402

# ── Device parameters (from AGENTS.md) ───────────────────────────────────────
OMEGA_Q = TWO_PI * 6.150e9   # rad/s
OMEGA_C = TWO_PI * 5.241e9   # rad/s
ALPHA   = TWO_PI * (-255e6)  # rad/s
CHI     = TWO_PI * (-2.84e6) # rad/s
CHI_P   = TWO_PI * (-21e3)   # rad/s
KERR    = TWO_PI * (-28e3)   # rad/s

N_CAV = 8
N_TR  = 2
N_SITES = 6    # number of cluster-state sites to simulate

# SNAP gate pulse optimization
SNAP_FIDELITY_THRESHOLD = 0.999
GRAPE_AMP_BOUND = TWO_PI * 50e6
GRAPE_MAXITER = 300
GRAPE_SEEDS = [17, 42, 73]

# Timing assumptions (from user specification)
ROTATION_DURATION_NS = 16.0
DISPLACEMENT_DURATION_NS = 48.0

# ── Helpers ──────────────────────────────────────────────────────────────────
def json_dump(path: Path, payload: Any) -> None:
    """Save JSON with numpy-safe serialization."""
    def default(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.float64, np.float32)):
            return float(o)
        if isinstance(o, (np.int64, np.int32)):
            return int(o)
        if isinstance(o, complex):
            return {"re": o.real, "im": o.imag}
        if isinstance(o, np.complex128):
            return {"re": float(o.real), "im": float(o.imag)}
        return str(o)
    path.write_text(json.dumps(payload, indent=2, sort_keys=False, default=default), encoding="utf-8")


def build_model(*, n_cav: int = N_CAV, n_tr: int = N_TR):
    return DispersiveTransmonCavityModel(
        omega_c=OMEGA_C, omega_q=OMEGA_Q, alpha=ALPHA,
        chi=CHI, chi_higher=(CHI_P,), kerr=KERR,
        n_cav=n_cav, n_tr=n_tr,
    )


def build_frame(model):
    return FrameSpec(omega_c_frame=model.omega_c, omega_q_frame=model.omega_q)


def logical_subspace(*, n_cav: int = N_CAV) -> Subspace:
    return Subspace.custom(
        full_dim=2 * n_cav,
        indices=(0, 1, n_cav, n_cav + 1),
        labels=("|g,0>", "|g,1>", "|e,0>", "|e,1>"),
    )


def state_fidelity(psi1: np.ndarray, psi2: np.ndarray) -> float:
    """Pure-state fidelity |<psi1|psi2>|^2."""
    return float(abs(np.vdot(psi1, psi2))**2)


def unitary_fidelity_4x4(U: np.ndarray, V: np.ndarray) -> float:
    """Process fidelity |Tr(U†V)|^2 / d^2 for 4×4 unitaries."""
    d = U.shape[0]
    return float(abs(np.trace(U.conj().T @ V))**2 / d**2)


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 1: Target Verification
# ═══════════════════════════════════════════════════════════════════════════════

def phase1_verify_target() -> dict[str, Any]:
    """Verify the cqed_sim cluster target against the canonical definition."""
    print("\n" + "=" * 70)
    print("  PHASE 1: Target Unitary Verification")
    print("=" * 70)
    results: dict[str, Any] = {}

    # 1a. Get the implemented target from cqed_sim
    U_impl = make_target("cluster", n_match=1)
    results["U_impl"] = U_impl
    print(f"  make_target('cluster', 1) shape: {U_impl.shape}")
    print(f"  Unitarity check: ||U†U - I|| = {np.linalg.norm(U_impl.conj().T @ U_impl - np.eye(4)):.2e}")

    # 1b. Construct the canonical per-site transfer matrix
    # Method: CZ_{q,c} · (H_q ⊗ I_c) — the "no-SWAP" version
    H = np.array([[1, 1], [1, -1]], dtype=np.complex128) / np.sqrt(2)
    CZ = np.diag([1.0, 1.0, 1.0, -1.0]).astype(np.complex128)
    SWAP = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], dtype=np.complex128)

    U_canonical_noswap = CZ @ np.kron(H, np.eye(2))  # CZ · (H_q ⊗ I_c)
    U_canonical_swap = SWAP @ CZ @ np.kron(H, np.eye(2))  # SWAP · CZ · (H_q ⊗ I_c)
    results["U_canonical_noswap"] = U_canonical_noswap
    results["U_canonical_swap"] = U_canonical_swap

    # Check if U_impl matches the SWAP version
    match_swap = np.allclose(U_impl, U_canonical_swap, atol=1e-12)
    match_noswap = np.allclose(U_impl, U_canonical_noswap, atol=1e-12)
    results["matches_swap_convention"] = bool(match_swap)
    results["matches_noswap_convention"] = bool(match_noswap)
    print(f"  Matches SWAP·CZ·(H⊗I): {match_swap}")
    print(f"  Matches CZ·(H⊗I) (no swap): {match_noswap}")

    # Process fidelity between the two conventions
    fid_impl_swap = unitary_fidelity_4x4(U_impl, U_canonical_swap)
    fid_impl_noswap = unitary_fidelity_4x4(U_impl, U_canonical_noswap)
    results["fidelity_impl_vs_swap"] = fid_impl_swap
    results["fidelity_impl_vs_noswap"] = fid_impl_noswap
    print(f"  Process fidelity (impl vs SWAP): {fid_impl_swap:.10f}")
    print(f"  Process fidelity (impl vs noSWAP): {fid_impl_noswap:.10f}")

    # 1c. Extract MPS tensors from both conventions and verify they produce
    #     the same cluster state
    # Convention 1 (no SWAP): qubit = physical, cavity = bond
    #   A^s_αβ = <s_q, β_c | U | 0_q, α_c>
    # Convention 2 (SWAP): after SWAP, cavity = physical output, qubit = bond
    #   B^s_αβ = <β_q, s_c | U_swap | 0_q, α_c>

    def extract_mps_tensors_noswap(U):
        """Extract MPS tensors with qubit = physical, cavity = bond."""
        A = {}
        for s in range(2):
            A[s] = np.zeros((2, 2), dtype=np.complex128)
            for alpha in range(2):
                for beta in range(2):
                    # |s, β> has flat index s*2 + β (qubit first)
                    # |0, α> has flat index 0*2 + α = α
                    A[s][alpha, beta] = U[s * 2 + beta, 0 * 2 + alpha]
        return A

    def extract_mps_tensors_swap(U):
        """Extract MPS tensors with cavity = physical (after SWAP), qubit = bond."""
        B = {}
        for s in range(2):
            B[s] = np.zeros((2, 2), dtype=np.complex128)
            for alpha in range(2):
                for beta in range(2):
                    # After SWAP: qubit=bond β, cavity=physical s
                    # |β_q, s_c> has flat index β*2 + s
                    # |0_q, α_c> has flat index 0*2 + α = α
                    B[s][alpha, beta] = U[beta * 2 + s, 0 * 2 + alpha]
        return B

    A_noswap = extract_mps_tensors_noswap(U_canonical_noswap)
    B_swap = extract_mps_tensors_swap(U_impl)

    # The canonical cluster-state MPS tensors are A^0 = I/√2, A^1 = Z/√2
    A0_expected = np.eye(2) / np.sqrt(2)
    A1_expected = np.diag([1, -1]).astype(np.complex128) / np.sqrt(2)

    mps_noswap_correct = (
        np.allclose(A_noswap[0], A0_expected, atol=1e-12) and
        np.allclose(A_noswap[1], A1_expected, atol=1e-12)
    )
    mps_swap_correct = (
        np.allclose(B_swap[0], A0_expected, atol=1e-12) and
        np.allclose(B_swap[1], A1_expected, atol=1e-12)
    )
    results["mps_tensors_noswap_correct"] = bool(mps_noswap_correct)
    results["mps_tensors_swap_correct"] = bool(mps_swap_correct)
    print(f"\n  MPS tensors (no-SWAP convention):")
    print(f"    A^0 = I/√2: {np.allclose(A_noswap[0], A0_expected, atol=1e-12)}")
    print(f"    A^1 = Z/√2: {np.allclose(A_noswap[1], A1_expected, atol=1e-12)}")
    print(f"  MPS tensors (SWAP convention, measuring cavity):")
    print(f"    B^0 = I/√2: {np.allclose(B_swap[0], A0_expected, atol=1e-12)}")
    print(f"    B^1 = Z/√2: {np.allclose(B_swap[1], A1_expected, atol=1e-12)}")

    # 1d. Generate cluster state via holographic protocol and compare
    #     Use the no-SWAP convention (simpler: qubit = physical, cavity = bond)
    #     and also the SWAP convention from make_target
    print(f"\n  Generating {N_SITES}-site cluster states via holographic protocol...")

    def generate_cluster_state_holographic(U, n_sites, convention="noswap"):
        """Generate N-site cluster state via repeated per-site unitary.

        Convention 'noswap': qubit = physical, cavity = bond.
          MPS tensor A^s = <s_q, ·| U |0_q, ·>
        Convention 'swap': cavity = physical, qubit = bond.
          MPS tensor B^s = <·_q, s_c| U |0_q, ·_c>
        """
        # Build the N-qubit state amplitude tensor
        # For each sequence of physical outcomes (s1, s2, ..., sN),
        # the amplitude is Tr(A^{s1} A^{s2} ... A^{sN}) for open boundary.
        # With open boundary: start from bond |0>, apply unitaries, trace out bond.
        if convention == "noswap":
            tensors = extract_mps_tensors_noswap(U)
        else:
            tensors = extract_mps_tensors_swap(U)

        # State in the full 2^N Hilbert space
        dim = 2**n_sites
        state = np.zeros(dim, dtype=np.complex128)
        for idx in range(dim):
            bits = [(idx >> (n_sites - 1 - k)) & 1 for k in range(n_sites)]
            # Open-boundary MPS: left boundary vector v_L = (1, 0) = bond |0>
            # Right boundary: trace over bond = sum over bond states
            v = np.array([1, 0], dtype=np.complex128)  # bond |0> initial
            for s in bits:
                v = tensors[s].T @ v  # A^s maps bond alpha -> bond beta
            # Open boundary: sum over final bond state (inner product with all bonds)
            state[idx] = np.sum(v)
        # Normalize
        norm = np.linalg.norm(state)
        if norm > 1e-15:
            state /= norm
        return state

    psi_noswap = generate_cluster_state_holographic(U_canonical_noswap, N_SITES, "noswap")
    psi_swap = generate_cluster_state_holographic(U_impl, N_SITES, "swap")

    # Generate canonical cluster state directly: (∏ CZ) H^N |0...0>
    def canonical_cluster_state(n_sites):
        """Build the canonical N-site 1D cluster state."""
        dim = 2**n_sites
        # Start with |0...0>
        state = np.zeros(dim, dtype=np.complex128)
        state[0] = 1.0

        # Apply H^⊗N
        H_single = np.array([[1, 1], [1, -1]], dtype=np.complex128) / np.sqrt(2)
        H_full = np.array([1.0], dtype=np.complex128)
        for _ in range(n_sites):
            H_full = np.kron(H_full, H_single)
        state = H_full @ state

        # Apply ∏ CZ_{i,i+1}
        for i in range(n_sites - 1):
            for idx in range(dim):
                bits = [(idx >> (n_sites - 1 - k)) & 1 for k in range(n_sites)]
                if bits[i] == 1 and bits[i + 1] == 1:
                    state[idx] *= -1
        return state

    psi_canonical = canonical_cluster_state(N_SITES)

    # Fidelities
    fid_noswap = state_fidelity(psi_noswap, psi_canonical)
    fid_swap = state_fidelity(psi_swap, psi_canonical)
    results["holostate_fidelity_noswap"] = fid_noswap
    results["holostate_fidelity_swap"] = fid_swap
    results["n_sites"] = N_SITES

    print(f"  F(holographic no-SWAP, canonical) = {fid_noswap:.12f}")
    print(f"  F(holographic SWAP,    canonical) = {fid_swap:.12f}")

    # Determine which convention to use going forward
    if fid_noswap > 0.999:
        print("  ✓ No-SWAP convention generates correct cluster state")
    if fid_swap > 0.999:
        print("  ✓ SWAP convention also generates correct cluster state")

    # Store the canonical cluster state for observable computation
    results["psi_canonical"] = psi_canonical
    results["psi_holographic"] = psi_swap if fid_swap > 0.999 else psi_noswap

    # The target unitary for decomposition is the SWAP version (from make_target)
    # since that's what cqed_sim implements. Both produce the correct cluster state.
    results["target_4x4"] = U_impl
    results["verification_passed"] = bool(fid_swap > 0.999 or fid_noswap > 0.999)

    print(f"\n  PHASE 1 RESULT: {'PASS' if results['verification_passed'] else 'FAIL'}")
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 2: Ideal Observables
# ═══════════════════════════════════════════════════════════════════════════════

def phase2_ideal_observables(psi: np.ndarray, n_sites: int) -> dict[str, Any]:
    """Compute all requested cluster-state observables."""
    print("\n" + "=" * 70)
    print("  PHASE 2: Ideal Cluster-State Observables")
    print("=" * 70)
    results: dict[str, Any] = {"n_sites": n_sites}
    dim = 2**n_sites

    # Pauli matrices
    I2 = np.eye(2, dtype=np.complex128)
    X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)

    def multi_qubit_op(op_list):
        """Build N-qubit operator from list of single-qubit ops."""
        result = np.array([1.0], dtype=np.complex128)
        for op in op_list:
            result = np.kron(result, op)
        return result

    def site_op(op, site, n):
        """Single-site operator embedded in N-qubit space."""
        ops = [I2] * n
        ops[site] = op
        return multi_qubit_op(ops)

    def expectation(op, state):
        return float(np.real(state.conj() @ op @ state))

    # 2a. Single-site Pauli expectations
    print("\n  2a. Single-site Pauli expectations:")
    pauli_table = []
    for i in range(n_sites):
        xi = expectation(site_op(X, i, n_sites), psi)
        yi = expectation(site_op(Y, i, n_sites), psi)
        zi = expectation(site_op(Z, i, n_sites), psi)
        pauli_table.append({"site": i, "X": xi, "Y": yi, "Z": zi})
        print(f"    site {i}: <X>={xi:+.8f}  <Y>={yi:+.8f}  <Z>={zi:+.8f}")
    results["pauli_expectations"] = pauli_table

    # Verify all are approximately zero
    max_pauli = max(max(abs(r["X"]), abs(r["Y"]), abs(r["Z"])) for r in pauli_table)
    results["max_single_site_pauli"] = max_pauli
    print(f"    Max |<σ>| = {max_pauli:.2e} (expected: 0)")

    # 2b. Stabilizer expectations ⟨Ki⟩
    print("\n  2b. Stabilizer expectations ⟨Ki⟩:")
    stabilizer_table = []
    for i in range(n_sites):
        ops = [I2] * n_sites
        ops[i] = X
        if i > 0:
            ops[i - 1] = Z
        if i < n_sites - 1:
            ops[i + 1] = Z
        K_i = multi_qubit_op(ops)
        ki = expectation(K_i, psi)
        label = ""
        if i == 0:
            label = f"X_{i} Z_{i+1}"
        elif i == n_sites - 1:
            label = f"Z_{i-1} X_{i}"
        else:
            label = f"Z_{i-1} X_{i} Z_{i+1}"
        stabilizer_table.append({"site": i, "K_i": ki, "label": label})
        print(f"    K_{i} = {label}: {ki:+.10f}")
    results["stabilizer_expectations"] = stabilizer_table

    # Verify all are +1
    max_stab_error = max(abs(r["K_i"] - 1.0) for r in stabilizer_table)
    results["max_stabilizer_deviation"] = max_stab_error
    print(f"    Max |⟨Ki⟩ - 1| = {max_stab_error:.2e} (expected: 0)")

    # 2c. ZXZ correlators (same as bulk stabilizers, just explicit)
    print("\n  2c. ZXZ correlators ⟨Z_{i-1} X_i Z_{i+1}⟩:")
    zxz_table = []
    for i in range(1, n_sites - 1):
        ops = [I2] * n_sites
        ops[i - 1] = Z
        ops[i] = X
        ops[i + 1] = Z
        zxz = expectation(multi_qubit_op(ops), psi)
        zxz_table.append({"site": i, "ZXZ": zxz})
        print(f"    Z_{i-1} X_{i} Z_{i+1} = {zxz:+.10f}")
    results["zxz_correlators"] = zxz_table

    # 2d. String-order correlators: <Z_i (∏_{k=i+1}^{j-1} X_k) Z_j>
    print("\n  2d. String-order correlators ⟨Z_i (∏X_k) Z_j⟩:")
    string_table = []
    for i in range(n_sites):
        for j in range(i + 2, n_sites):
            ops = [I2] * n_sites
            ops[i] = Z
            ops[j] = Z
            for k in range(i + 1, j):
                ops[k] = X
            val = expectation(multi_qubit_op(ops), psi)
            string_table.append({"i": i, "j": j, "value": val})
            label = f"Z_{i} " + " ".join(f"X_{k}" for k in range(i+1, j)) + f" Z_{j}"
            print(f"    {label} = {val:+.10f}")
    results["string_order_correlators"] = string_table

    # For ideal cluster state, string-order correlators should all be 0 for
    # |i-j| > 2 and +1 for the stabilizer (|i-j| = 2)
    # Actually let me check: for the cluster state, the string correlator
    # <Z_i X_{i+1} Z_{i+2}> = 1 (stabilizer). Longer strings:
    # <Z_i X_{i+1} X_{i+2} Z_{i+3}> = <Z_i X_{i+1} Z_{i+2}> * 0 = ?
    # This needs more careful analysis. Let me just report the numerical values.

    print(f"\n  PHASE 2 COMPLETE")
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 3: Decomposition
# ═══════════════════════════════════════════════════════════════════════════════

def phase3_decomposition(target_4x4: np.ndarray) -> dict[str, Any]:
    """Decompose the per-site unitary into QubitRotation + Displacement + SNAP."""
    print("\n" + "=" * 70)
    print("  PHASE 3: Gate-Set Decomposition")
    print("=" * 70)
    results: dict[str, Any] = {}

    n_cav = N_CAV
    subspace = logical_subspace(n_cav=n_cav)

    # Strategy: Use UnitarySynthesizer with SNAP-based gate set
    # Template: D + Rq + SNAP + D + Rq + SNAP + D + Rq
    # This allows displacement to re-arrange cavity populations,
    # SNAP to apply number-selective phases,
    # and Rq to handle qubit rotations.
    DISP_S = DISPLACEMENT_DURATION_NS * 1e-9
    ROT_S = ROTATION_DURATION_NS * 1e-9
    SNAP_S = 200e-9  # initial guess for SNAP duration

    # Try different decomposition depths
    best_result = None
    best_fidelity = 0.0
    best_n_snap = 0

    for n_snap in [1, 2, 3]:
        print(f"\n  Trying decomposition with {n_snap} SNAP gate(s)...")
        gates = []
        # Build: (D + Rq + SNAP) × n_snap + D + Rq
        for k in range(n_snap):
            gates.append(SynthDisplacement(
                name=f"D{k+1}", alpha=0.1+0.0j,
                duration=DISP_S, optimize_time=False,
            ))
            gates.append(SynthQubitRotation(
                name=f"Rq{k+1}", theta=0.1, phi=0.0,
                duration=ROT_S, optimize_time=False,
            ))
            gates.append(SynthSNAP(
                name=f"SNAP{k+1}",
                phases=[0.0] * n_cav,
                duration=SNAP_S,
                optimize_time=False,
            ))
        # Final displacement and rotation
        gates.append(SynthDisplacement(
            name=f"D{n_snap+1}", alpha=0.1+0.0j,
            duration=DISP_S, optimize_time=False,
        ))
        gates.append(SynthQubitRotation(
            name=f"Rq{n_snap+1}", theta=0.1, phi=0.0,
            duration=ROT_S, optimize_time=False,
        ))

        target = TargetUnitary(target_4x4, ignore_global_phase=True)
        try:
            synth = UnitarySynthesizer(
                primitives=gates,
                subspace=subspace,
                target=target,
                seed=42,
                optimize_times=False,
                optimizer="powell",
            )
            result = synth.fit(multistart=5, maxiter=500)
            fid = float(result.fidelity)
            print(f"    {n_snap}-SNAP fidelity: {fid:.8f}")

            if fid > best_fidelity:
                best_fidelity = fid
                best_result = result
                best_n_snap = n_snap

            if fid > 0.999:
                print(f"    Sufficient fidelity reached with {n_snap} SNAP gate(s)")
                break
        except Exception as exc:
            print(f"    Failed: {exc}")
            traceback.print_exc()

    if best_result is None:
        print("  ERROR: No decomposition succeeded")
        return {"error": "No decomposition succeeded"}

    results["n_snap_gates"] = best_n_snap
    results["ideal_decomposition_fidelity"] = best_fidelity
    seq = best_result.sequence

    # Extract gate parameters
    gate_params = []
    snap_params_list = []
    for gate in seq.gates:
        gate_type = type(gate).__name__
        entry = {"name": gate.name, "type": gate_type}
        if gate_type == "SNAP":
            phases = list(gate.phases) if hasattr(gate, 'phases') else []
            entry["phases"] = phases
            snap_params_list.append({"name": gate.name, "phases": phases})
        elif gate_type == "QubitRotation":
            entry["theta"] = float(gate.theta)
            entry["phi"] = float(gate.phi)
        elif gate_type == "Displacement":
            entry["alpha_re"] = float(np.real(gate.alpha))
            entry["alpha_im"] = float(np.imag(gate.alpha))
        gate_params.append(entry)

    results["gate_sequence"] = gate_params
    results["snap_gates"] = snap_params_list
    results["total_gates"] = len(seq.gates)
    results["n_rotations"] = sum(1 for g in seq.gates if type(g).__name__ == "QubitRotation")
    results["n_displacements"] = sum(1 for g in seq.gates if type(g).__name__ == "Displacement")

    print(f"\n  Best decomposition: {best_n_snap} SNAP, fidelity = {best_fidelity:.8f}")
    print(f"  Gate sequence ({len(seq.gates)} gates):")
    for entry in gate_params:
        if entry["type"] == "SNAP":
            phases_str = ", ".join(f"{p:.4f}" for p in entry["phases"][:4])
            print(f"    {entry['name']}: phases=[{phases_str}, ...]")
        elif entry["type"] == "QubitRotation":
            print(f"    {entry['name']}: θ={entry['theta']:.4f}, φ={entry['phi']:.4f}")
        elif entry["type"] == "Displacement":
            print(f"    {entry['name']}: α={entry['alpha_re']:.4f}+{entry['alpha_im']:.4f}j")

    # Compute decomposition-level observables
    # Use the decomposed unitary to prepare the cluster state and check observables
    full_U = np.asarray(seq.unitary(backend="ideal"), dtype=np.complex128)
    sub_U = subspace.restrict_operator(full_U)
    results["subspace_unitary"] = sub_U

    # Store the synthesis result for later use
    results["synthesis_result"] = best_result
    results["gate_sequence_obj"] = seq

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 4: SNAP Pulse Optimization
# ═══════════════════════════════════════════════════════════════════════════════

def phase4_snap_optimization(snap_gates: list[dict], decomp_results: dict) -> dict[str, Any]:
    """GRAPE-optimize each SNAP gate for minimum pulse duration."""
    print("\n" + "=" * 70)
    print("  PHASE 4: SNAP Gate Pulse Optimization")
    print("=" * 70)
    results: dict[str, Any] = {}

    model = build_model()
    frame = build_frame(model)
    subspace = logical_subspace()

    # For each SNAP gate, build its ideal unitary and sweep duration
    for snap_info in snap_gates:
        snap_name = snap_info["name"]
        phases = snap_info["phases"]
        print(f"\n  Optimizing {snap_name}: phases = {[f'{p:.3f}' for p in phases[:4]]}...")

        # Build the ideal SNAP operator in the full Hilbert space
        # SNAP acts on cavity only: diag(e^{iφ_0}, e^{iφ_1}, ...)
        # Embedded: I_qubit ⊗ SNAP_cavity
        snap_cavity = snap(phases[:N_CAV], dim=N_CAV)
        snap_full = qt.tensor(qt.qeye(N_TR), snap_cavity)
        snap_matrix = np.asarray(snap_full.full(), dtype=np.complex128)

        # Restrict to logical subspace for target
        snap_sub = subspace.restrict_operator(snap_matrix)

        # Duration sweep: try progressively shorter durations
        durations_ns = [500, 400, 300, 250, 200, 150, 120, 100, 80, 60]
        sweep_results = []

        for dur_ns in durations_ns:
            dur_s = dur_ns * 1e-9
            n_slices = max(10, int(dur_ns / 4))  # 4 ns per slice
            dt_s = dur_s / n_slices

            best_fid_this_dur = 0.0
            best_run = None

            for seed in GRAPE_SEEDS:
                try:
                    problem = build_control_problem_from_model(
                        model, frame=frame,
                        time_grid=PiecewiseConstantTimeGrid.uniform(steps=n_slices, dt_s=dt_s),
                        channel_specs=(
                            ModelControlChannelSpec(
                                name="qubit", target="qubit",
                                quadratures=("I", "Q"),
                                amplitude_bounds=(-GRAPE_AMP_BOUND, GRAPE_AMP_BOUND),
                                export_channel="qubit",
                            ),
                        ),
                        objectives=(
                            OCUnitaryObjective(
                                target_operator=snap_sub,
                                subspace=subspace,
                                ignore_global_phase=True,
                                name=f"{snap_name}_dur{dur_ns}",
                            ),
                        ),
                        penalties=(
                            OCLeakagePenalty(weight=0.05, subspace=subspace),
                        ),
                    )
                    result = GrapeSolver(GrapeConfig(
                        maxiter=GRAPE_MAXITER, seed=seed, random_scale=0.3
                    )).solve(problem)

                    fid = float(result.metrics.get("nominal_fidelity",
                                result.metrics.get("fidelity", 0.0)))
                    if fid > best_fid_this_dur:
                        best_fid_this_dur = fid
                        best_run = {
                            "seed": seed,
                            "fidelity": fid,
                            "converged": bool(result.metrics.get("converged", False)),
                            "iterations": int(result.metrics.get("iterations",
                                             result.metrics.get("n_iter", 0))),
                        }
                except Exception as exc:
                    print(f"      Seed {seed} @ {dur_ns}ns failed: {exc}")

            entry = {
                "duration_ns": dur_ns,
                "n_slices": n_slices,
                "best_fidelity": best_fid_this_dur,
                "above_threshold": best_fid_this_dur >= SNAP_FIDELITY_THRESHOLD,
            }
            if best_run:
                entry.update(best_run)
            sweep_results.append(entry)
            status = "✓" if best_fid_this_dur >= SNAP_FIDELITY_THRESHOLD else "✗"
            print(f"    {dur_ns:4d} ns: F = {best_fid_this_dur:.6f} {status}")

            # Stop if we've found the minimum duration that works
            if best_fid_this_dur < SNAP_FIDELITY_THRESHOLD and len(sweep_results) > 1:
                if sweep_results[-2]["above_threshold"]:
                    print(f"    → Minimum duration: {sweep_results[-2]['duration_ns']} ns")
                    break

        # Find the minimum duration above threshold
        passing = [r for r in sweep_results if r["above_threshold"]]
        if passing:
            min_dur = min(r["duration_ns"] for r in passing)
            best_entry = min((r for r in passing), key=lambda r: r["duration_ns"])
        else:
            min_dur = sweep_results[0]["duration_ns"]
            best_entry = max(sweep_results, key=lambda r: r["best_fidelity"])

        results[snap_name] = {
            "phases": phases,
            "duration_sweep": sweep_results,
            "minimum_duration_ns": min_dur,
            "best_fidelity": best_entry["best_fidelity"],
            "optimized_params": best_entry,
        }

    results["fidelity_threshold"] = SNAP_FIDELITY_THRESHOLD
    print(f"\n  PHASE 4 COMPLETE")
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 5: Timing Analysis
# ═══════════════════════════════════════════════════════════════════════════════

def phase5_timing(decomp_results: dict, snap_results: dict) -> dict[str, Any]:
    """Compute SNAP-dominated and full-sequence timing summaries."""
    print("\n" + "=" * 70)
    print("  PHASE 5: Timing Analysis")
    print("=" * 70)
    results: dict[str, Any] = {}

    n_rot = decomp_results["n_rotations"]
    n_disp = decomp_results["n_displacements"]
    n_snap = decomp_results["n_snap_gates"]

    # Collect optimized SNAP durations
    snap_durations = []
    for snap_info in decomp_results["snap_gates"]:
        snap_name = snap_info["name"]
        if snap_name in snap_results:
            dur = snap_results[snap_name]["minimum_duration_ns"]
        else:
            dur = 200.0  # fallback
        snap_durations.append({"name": snap_name, "duration_ns": dur})

    total_snap_ns = sum(d["duration_ns"] for d in snap_durations)
    total_rot_ns = n_rot * ROTATION_DURATION_NS
    total_disp_ns = n_disp * DISPLACEMENT_DURATION_NS
    total_full_ns = total_snap_ns + total_rot_ns + total_disp_ns

    results["snap_dominated"] = {
        "total_snap_time_ns": total_snap_ns,
        "per_snap_durations": snap_durations,
    }
    results["full_sequence"] = {
        "n_qubit_rotations": n_rot,
        "n_displacements": n_disp,
        "n_snap_gates": n_snap,
        "rotation_duration_each_ns": ROTATION_DURATION_NS,
        "displacement_duration_each_ns": DISPLACEMENT_DURATION_NS,
        "total_rotation_time_ns": total_rot_ns,
        "total_displacement_time_ns": total_disp_ns,
        "total_snap_time_ns": total_snap_ns,
        "total_sequence_time_ns": total_full_ns,
    }

    print(f"\n  5.1 SNAP-Dominated Timing:")
    print(f"    Number of SNAP gates: {n_snap}")
    for d in snap_durations:
        print(f"      {d['name']}: {d['duration_ns']:.1f} ns")
    print(f"    Total SNAP time: {total_snap_ns:.1f} ns")

    print(f"\n  5.2 Full Sequence Timing:")
    print(f"    {n_rot} × qubit rotation @ {ROTATION_DURATION_NS} ns = {total_rot_ns:.1f} ns")
    print(f"    {n_disp} × displacement @ {DISPLACEMENT_DURATION_NS} ns = {total_disp_ns:.1f} ns")
    print(f"    {n_snap} × SNAP gate (optimized) = {total_snap_ns:.1f} ns")
    print(f"    ──────────────────────────────────")
    print(f"    Total sequence time: {total_full_ns:.1f} ns")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 6: Waveform-Level Simulation
# ═══════════════════════════════════════════════════════════════════════════════

def phase6_waveform_simulation(target_4x4: np.ndarray, decomp_results: dict) -> dict[str, Any]:
    """Run pulse-level simulation of the decomposed sequence."""
    print("\n" + "=" * 70)
    print("  PHASE 6: Waveform-Level Simulation")
    print("=" * 70)
    results: dict[str, Any] = {}

    model = build_model()
    frame = build_frame(model)
    subspace = logical_subspace()
    seq = decomp_results.get("gate_sequence_obj")

    if seq is None:
        print("  No gate sequence available for waveform simulation")
        return {"error": "No gate sequence"}

    # Try waveform bridge
    try:
        from cqed_sim.unitary_synthesis.waveform_bridge import waveform_sequence_from_gates
        from cqed_sim.unitary_synthesis.systems import CQEDSystemAdapter

        system = CQEDSystemAdapter(model=model)
        wf_seq = waveform_sequence_from_gates(seq, frame=frame)

        # Simulate with pulse-level
        from cqed_sim.unitary_synthesis.backends import simulate_sequence as synth_simulate
        sim_result = synth_simulate(
            wf_seq, subspace,
            backend="cqed",
            system=system,
            need_operator=True,
            dt=4e-9,
            frame=frame,
        )

        if sim_result.subspace_operator is not None:
            pulse_sub_U = sim_result.subspace_operator
            pulse_fid = float(subspace_unitary_fidelity(
                pulse_sub_U, target_4x4, gauge="global"
            ))
            results["pulse_level_fidelity"] = pulse_fid
            print(f"  Pulse-level subspace fidelity: {pulse_fid:.8f}")
        else:
            results["pulse_level_fidelity"] = None
            print("  No subspace operator from pulse simulation")

        results["waveform_bridge_success"] = True
    except Exception as exc:
        print(f"  Waveform bridge failed: {exc}")
        traceback.print_exc()
        results["waveform_bridge_success"] = False
        results["waveform_bridge_error"] = str(exc)

    # Alternative: run GRAPE on the full per-site unitary directly
    print("\n  Running GRAPE for full per-site unitary (reference)...")
    try:
        n_slices = 50
        dt_s = 200e-9 / n_slices  # 200 ns total, 4 ns slices
        problem = build_control_problem_from_model(
            model, frame=frame,
            time_grid=PiecewiseConstantTimeGrid.uniform(steps=n_slices, dt_s=dt_s),
            channel_specs=(
                ModelControlChannelSpec(
                    name="storage", target="storage",
                    quadratures=("I", "Q"),
                    amplitude_bounds=(-GRAPE_AMP_BOUND, GRAPE_AMP_BOUND),
                    export_channel="storage",
                ),
                ModelControlChannelSpec(
                    name="qubit", target="qubit",
                    quadratures=("I", "Q"),
                    amplitude_bounds=(-GRAPE_AMP_BOUND, GRAPE_AMP_BOUND),
                    export_channel="qubit",
                ),
            ),
            objectives=(
                OCUnitaryObjective(
                    target_operator=target_4x4,
                    subspace=subspace,
                    ignore_global_phase=True,
                    name="full_per_site_grape",
                ),
            ),
            penalties=(
                OCLeakagePenalty(weight=0.02, subspace=subspace),
            ),
        )
        best_grape_fid = 0.0
        for seed in GRAPE_SEEDS:
            result = GrapeSolver(GrapeConfig(
                maxiter=GRAPE_MAXITER, seed=seed, random_scale=0.3
            )).solve(problem)
            fid = float(result.metrics.get("nominal_fidelity",
                        result.metrics.get("fidelity", 0.0)))
            if fid > best_grape_fid:
                best_grape_fid = fid
                best_grape_result = result

        results["grape_full_unitary_fidelity"] = best_grape_fid
        print(f"  GRAPE full-unitary fidelity: {best_grape_fid:.8f}")
    except Exception as exc:
        print(f"  GRAPE full-unitary failed: {exc}")
        results["grape_full_unitary_error"] = str(exc)

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# OBSERVABLE COMPARISON (Decomposed vs Ideal)
# ═══════════════════════════════════════════════════════════════════════════════

def compare_observables_decomposed(target_4x4: np.ndarray, sub_U: np.ndarray,
                                     n_sites: int) -> dict[str, Any]:
    """Compare observables from the decomposed unitary against ideal."""
    print("\n" + "=" * 70)
    print("  Observable Comparison: Decomposed vs Ideal")
    print("=" * 70)

    I2 = np.eye(2, dtype=np.complex128)
    X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)

    def multi_qubit_op(op_list):
        result = np.array([1.0], dtype=np.complex128)
        for op in op_list:
            result = np.kron(result, op)
        return result

    def site_op(op, site, n):
        ops = [I2] * n
        ops[site] = op
        return multi_qubit_op(ops)

    def expectation(op, state):
        return float(np.real(state.conj() @ op @ state))

    # Generate cluster state using the decomposed unitary
    def extract_mps_tensors_swap(U):
        B = {}
        for s in range(2):
            B[s] = np.zeros((2, 2), dtype=np.complex128)
            for alpha in range(2):
                for beta in range(2):
                    B[s][alpha, beta] = U[beta * 2 + s, 0 * 2 + alpha]
        return B

    tensors = extract_mps_tensors_swap(sub_U)

    dim = 2**n_sites
    state = np.zeros(dim, dtype=np.complex128)
    for idx in range(dim):
        bits = [(idx >> (n_sites - 1 - k)) & 1 for k in range(n_sites)]
        v = np.array([1, 0], dtype=np.complex128)
        for s in bits:
            v = tensors[s].T @ v
        state[idx] = np.sum(v)
    norm = np.linalg.norm(state)
    if norm > 1e-15:
        state /= norm

    results = {"n_sites": n_sites}

    # Stabilizers
    stab_table = []
    for i in range(n_sites):
        ops = [I2] * n_sites
        ops[i] = X
        if i > 0:
            ops[i - 1] = Z
        if i < n_sites - 1:
            ops[i + 1] = Z
        K_i = multi_qubit_op(ops)
        ki = expectation(K_i, state)
        stab_table.append({"site": i, "K_i": ki})
        print(f"    ⟨K_{i}⟩ = {ki:+.8f}")
    results["stabilizer_expectations"] = stab_table

    # Single-site Pauli
    pauli_table = []
    for i in range(n_sites):
        xi = expectation(site_op(X, i, n_sites), state)
        yi = expectation(site_op(Y, i, n_sites), state)
        zi = expectation(site_op(Z, i, n_sites), state)
        pauli_table.append({"site": i, "X": xi, "Y": yi, "Z": zi})
    results["pauli_expectations"] = pauli_table

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE GENERATION
# ═══════════════════════════════════════════════════════════════════════════════

def make_figures(all_results: dict) -> None:
    """Generate publication-quality figures."""
    print("\n  Generating figures...")

    # Fig 1: Stabilizer expectations (ideal vs decomposed)
    fig, ax = plt.subplots(figsize=(8, 4))
    ideal_stab = all_results.get("phase2", {}).get("stabilizer_expectations", [])
    decomp_stab = all_results.get("decomp_observables", {}).get("stabilizer_expectations", [])
    if ideal_stab:
        sites = [s["site"] for s in ideal_stab]
        ideal_vals = [s["K_i"] for s in ideal_stab]
        ax.plot(sites, ideal_vals, "o-", color=TOL_BRIGHT[0], label="Ideal", markersize=8)
    if decomp_stab:
        sites_d = [s["site"] for s in decomp_stab]
        decomp_vals = [s["K_i"] for s in decomp_stab]
        ax.plot(sites_d, decomp_vals, "s--", color=TOL_BRIGHT[1], label="Decomposed", markersize=6)
    ax.axhline(1.0, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("Site $i$")
    ax.set_ylabel(r"$\langle K_i \rangle$")
    ax.set_title("Cluster-State Stabilizer Expectations")
    ax.legend()
    ax.set_ylim(0.9, 1.05)
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"stabilizer_expectations.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Fig 2: SNAP duration sweep
    snap_data = all_results.get("phase4", {})
    snap_keys = [k for k in snap_data if k.startswith("SNAP")]
    if snap_keys:
        fig, ax = plt.subplots(figsize=(8, 4))
        for i, key in enumerate(snap_keys):
            sweep = snap_data[key]["duration_sweep"]
            durs = [s["duration_ns"] for s in sweep]
            fids = [s["best_fidelity"] for s in sweep]
            ax.plot(durs, fids, "o-", color=TOL_BRIGHT[i % len(TOL_BRIGHT)], label=key)
        ax.axhline(SNAP_FIDELITY_THRESHOLD, color="red", linestyle="--",
                   label=f"Threshold ({SNAP_FIDELITY_THRESHOLD})")
        ax.set_xlabel("Pulse Duration (ns)")
        ax.set_ylabel("GRAPE Fidelity")
        ax.set_title("SNAP Gate Pulse Optimization")
        ax.legend()
        for ext in ("png", "pdf"):
            fig.savefig(FIG_DIR / f"snap_duration_sweep.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # Fig 3: Single-site Pauli expectations
    pauli_data = all_results.get("phase2", {}).get("pauli_expectations", [])
    if pauli_data:
        fig, ax = plt.subplots(figsize=(8, 4))
        sites = [d["site"] for d in pauli_data]
        ax.plot(sites, [d["X"] for d in pauli_data], "o-", color=TOL_BRIGHT[0], label=r"$\langle X_i \rangle$")
        ax.plot(sites, [d["Y"] for d in pauli_data], "s-", color=TOL_BRIGHT[1], label=r"$\langle Y_i \rangle$")
        ax.plot(sites, [d["Z"] for d in pauli_data], "^-", color=TOL_BRIGHT[2], label=r"$\langle Z_i \rangle$")
        ax.axhline(0.0, color="gray", linestyle=":", alpha=0.5)
        ax.set_xlabel("Site $i$")
        ax.set_ylabel("Expectation Value")
        ax.set_title("Single-Site Pauli Expectations (Ideal Cluster State)")
        ax.legend()
        for ext in ("png", "pdf"):
            fig.savefig(FIG_DIR / f"pauli_expectations.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # Fig 4: String-order correlators heatmap
    string_data = all_results.get("phase2", {}).get("string_order_correlators", [])
    n_sites = all_results.get("phase2", {}).get("n_sites", 6)
    if string_data:
        matrix = np.full((n_sites, n_sites), np.nan)
        for entry in string_data:
            matrix[entry["i"], entry["j"]] = entry["value"]
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(matrix, cmap="RdBu_r", vmin=-1, vmax=1, origin="upper")
        ax.set_xlabel("Site $j$")
        ax.set_ylabel("Site $i$")
        ax.set_title(r"String-Order Correlators $\langle Z_i (\prod X_k) Z_j \rangle$")
        plt.colorbar(im, ax=ax, label="Correlator Value")
        for ext in ("png", "pdf"):
            fig.savefig(FIG_DIR / f"string_order_correlators.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    print("  Figures saved.")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    t0 = time.perf_counter()
    all_results = {}

    # Phase 1: Target verification
    p1 = phase1_verify_target()
    all_results["phase1"] = {k: v for k, v in p1.items()
                              if not isinstance(v, np.ndarray) or v.size < 100}
    all_results["phase1"]["verification_passed"] = p1["verification_passed"]

    if not p1["verification_passed"]:
        print("\n  CRITICAL: Target verification failed. Cannot proceed.")
        json_dump(DATA_DIR / "results.json", all_results)
        return

    # Phase 2: Ideal observables
    psi = p1["psi_canonical"]
    p2 = phase2_ideal_observables(psi, N_SITES)
    all_results["phase2"] = p2

    # Phase 3: Decomposition
    p3 = phase3_decomposition(p1["target_4x4"])
    all_results["phase3"] = {k: v for k, v in p3.items()
                              if k not in ("synthesis_result", "gate_sequence_obj", "subspace_unitary")}

    if "error" in p3:
        print("\n  CRITICAL: Decomposition failed. Cannot proceed to SNAP optimization.")
        json_dump(DATA_DIR / "results.json", all_results)
        return

    # Decomposition-level observables
    if "subspace_unitary" in p3:
        decomp_obs = compare_observables_decomposed(p1["target_4x4"], p3["subspace_unitary"], N_SITES)
        all_results["decomp_observables"] = decomp_obs

    # Phase 4: SNAP optimization
    p4 = phase4_snap_optimization(p3["snap_gates"], p3)
    all_results["phase4"] = p4

    # Phase 5: Timing
    p5 = phase5_timing(p3, p4)
    all_results["phase5"] = p5

    # Phase 6: Waveform simulation
    p6 = phase6_waveform_simulation(p1["target_4x4"], p3)
    all_results["phase6"] = {k: v for k, v in p6.items()
                              if not isinstance(v, np.ndarray) or v.size < 100}

    # Generate figures
    make_figures(all_results)

    # Save all results
    json_dump(DATA_DIR / "results.json", all_results)

    elapsed = time.perf_counter() - t0
    print(f"\n{'=' * 70}")
    print(f"  ALL PHASES COMPLETE in {elapsed:.1f}s")
    print(f"{'=' * 70}")

    # Print summary
    print(f"\n  === SUMMARY ===")
    print(f"  Target verification: {'PASS' if p1['verification_passed'] else 'FAIL'}")
    if "ideal_decomposition_fidelity" in p3:
        print(f"  Decomposition fidelity: {p3['ideal_decomposition_fidelity']:.8f}")
        print(f"  Number of SNAP gates: {p3['n_snap_gates']}")
    if "full_sequence" in p5:
        fs = p5["full_sequence"]
        print(f"  Total SNAP time: {fs['total_snap_time_ns']:.1f} ns")
        print(f"  Total sequence time: {fs['total_sequence_time_ns']:.1f} ns")
    if "grape_full_unitary_fidelity" in p6:
        print(f"  GRAPE full-unitary fidelity: {p6['grape_full_unitary_fidelity']:.8f}")


if __name__ == "__main__":
    main()
