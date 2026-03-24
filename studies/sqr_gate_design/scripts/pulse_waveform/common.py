"""
Common definitions for the SQR pulse-waveform design study.

Physical setup:
  - Dispersive transmon-cavity system with g/e/f transmon levels.
  - Truncated Fock space n = 0, ..., N_FOCK-1.
  - χ = 2π × (−2.84 MHz), α = 2π × (−255 MHz).
  - Dimensionless scan variable: χT/(2π) = |χ/(2π)| · T.

All frequencies in rad/s, times in seconds.
"""

import numpy as np

import runtime_compat  # noqa: F401

# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------
OMEGA_Q = 2 * np.pi * 6.150e9      # qubit frequency (rad/s)
OMEGA_C = 2 * np.pi * 5.241e9      # cavity frequency (rad/s)
ALPHA   = 2 * np.pi * (-255e6)     # anharmonicity (rad/s)
CHI     = 2 * np.pi * (-2.84e6)    # dispersive shift (rad/s)
CHI_PRIME = 2 * np.pi * (-21e3)   # second-order dispersive (rad/s), Phase 4
KERR    = 2 * np.pi * (-28e3)     # cavity self-Kerr (rad/s), Phase 4

# ---------------------------------------------------------------------------
# Decoherence parameters
# ---------------------------------------------------------------------------
T1 = 20e-6        # transmon relaxation time (s)
T2 = 20e-6        # transmon dephasing time (s)
T_PHI = 1.0 / (1.0 / T2 - 1.0 / (2 * T1))  # pure dephasing time (s) = 40 μs

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------
N_FOCK  = 4       # number of Fock levels in truncated subspace
N_CAV   = N_FOCK + 2   # cavity Hilbert space dimension (buffer)
N_TR    = 3       # transmon levels (g, e, f) for leakage modeling
DT      = 2e-9    # simulation time step (s)

# ---------------------------------------------------------------------------
# Scan parameters
# ---------------------------------------------------------------------------
# Scan variable: χT/(2π) = |χ/(2π)| · T  (number of chi-periods)
CHI_T_VALUES = np.array([0.5, 1, 1.5, 2, 3, 5, 7, 10], dtype=float)

# Target branch for SQR
TARGET_N0 = 1     # target Fock level

# Target rotation angles
THETA_TARGET = np.pi   # π rotation (X_π gate)
PHI_TARGET   = 0.0     # rotation axis azimuth (X axis)

# ---------------------------------------------------------------------------
# Colorblind-friendly palette (Tol's Bright)
# ---------------------------------------------------------------------------
TOL_BRIGHT = [
    '#4477AA', '#EE6677', '#228833', '#CCBB44',
    '#66CCEE', '#AA3377', '#BBBBBB',
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def duration_from_chi_t(chi_t_2pi, chi=CHI):
    """Gate duration T such that |χ/(2π)| · T = chi_t_2pi."""
    f_chi = abs(chi) / (2 * np.pi)  # Hz
    return chi_t_2pi / f_chi


def build_model(chi=CHI, chi_prime=0.0, kerr=0.0, n_cav=N_CAV, n_tr=N_TR):
    """Build a DispersiveTransmonCavityModel with the given parameters."""
    from cqed_sim.core import DispersiveTransmonCavityModel
    chi_higher = (chi_prime,) if chi_prime != 0.0 else ()
    return DispersiveTransmonCavityModel(
        omega_c=OMEGA_C,
        omega_q=OMEGA_Q,
        alpha=ALPHA,
        chi=chi,
        chi_higher=chi_higher,
        kerr=kerr,
        n_cav=n_cav,
        n_tr=n_tr,
    )


def build_noise_spec(t1=T1, tphi=T_PHI):
    """Build a NoiseSpec for transmon T1 relaxation and pure dephasing."""
    from cqed_sim.sim.noise import NoiseSpec
    return NoiseSpec(t1=t1, tphi=tphi)


def state_fidelity_dm(target_ket, rho):
    """State fidelity F = <ψ|ρ|ψ> for a pure target |ψ> and density matrix ρ.

    Also works if rho is a ket (returns |<ψ|ρ>|²).
    """
    if rho.isoper:
        t = target_ket.full().flatten()
        r = rho.full()
        return float(np.real(t.conj() @ r @ t))
    else:
        return float(abs(target_ket.dag() * rho) ** 2)


def build_frame(model):
    """Rotating frame at qubit and cavity bare frequencies."""
    from cqed_sim.core import FrameSpec
    return FrameSpec(omega_c_frame=model.omega_c, omega_q_frame=model.omega_q)


def branch_frequencies(model, frame, n_max):
    """Return array of rotating-frame qubit transition frequencies for n=0..n_max-1."""
    from cqed_sim.core.frequencies import manifold_transition_frequency
    return np.array([
        manifold_transition_frequency(model, n, frame)
        for n in range(n_max)
    ])


def target_qubit_unitary(theta, phi):
    """2×2 target qubit rotation R(θ, φ) = exp(-i θ/2 [cos φ σ_x + sin φ σ_y])."""
    c = np.cos(theta / 2)
    s = np.sin(theta / 2)
    return np.array([
        [c, -1j * s * np.exp(-1j * phi)],
        [-1j * s * np.exp(1j * phi), c],
    ], dtype=np.complex128)


def conditional_process_fidelity(U_actual, U_target):
    """Process fidelity |Tr(V† U)|² / 4 for 2×2 unitary blocks."""
    return abs(np.trace(U_target.conj().T @ U_actual))**2 / 4.0


def z_corrected_target_fidelity(U_actual, U_target):
    """Process fidelity maximized over a qubit Z-rotation.

    Computes max_α |Tr(U_target† Z(α) U_actual)|² / d²
    where Z(α) = diag(1, e^{iα}).

    The optimal α removes the dynamical Z-phase accumulated during the pulse,
    which in experiment is handled by a virtual-Z frame update.

    Returns (fidelity, alpha_opt).
    """
    Uc = U_target.conj()
    # Tr(U†ZU_a) = A + B*e^{iα}, max |Tr|² = (|A|+|B|)² at α = arg(A)-arg(B)
    A = Uc[0, 0] * U_actual[0, 0] + Uc[0, 1] * U_actual[0, 1]
    B = Uc[1, 0] * U_actual[1, 0] + Uc[1, 1] * U_actual[1, 1]
    alpha_opt = np.angle(A) - np.angle(B) if abs(B) > 1e-15 else 0.0
    fid = (abs(A) + abs(B))**2 / 4.0
    return fid, alpha_opt


def identity_fidelity_with_z(U_n, alpha):
    """Fidelity of U_n vs identity after applying global Z(alpha) = diag(1, e^{iα}).

    F = |Tr(Z(α) U_n)|² / 4 = |U_n[0,0] + e^{iα} U_n[1,1]|² / 4.
    """
    return abs(U_n[0, 0] + np.exp(1j * alpha) * U_n[1, 1])**2 / 4.0


def spectator_z_fidelity(U_n):
    """Best-fit Z-rotation fidelity for a spectator branch.

    Returns (fidelity, best_phase) where fidelity = max_φ |Tr(Z(φ)† U_n)|²/4.
    """
    # For Z(φ) = diag(e^{-iφ/2}, e^{iφ/2}), the overlap is
    # Tr(Z†U) = U_00 e^{iφ/2} + U_11 e^{-iφ/2}
    # Maximize |Tr|² over φ → max = (|U_00| + |U_11|)².
    u00 = U_n[0, 0]
    u11 = U_n[1, 1]
    phi_opt = np.angle(u11) - np.angle(u00)
    fid = (abs(u00) + abs(u11))**2 / 4.0
    return fid, phi_opt


def spectator_transverse_error(U_n):
    """Transverse error: 1 - max_φ F(U_n, Z(φ)).  Measures non-Z content."""
    fid, _ = spectator_z_fidelity(U_n)
    return 1.0 - fid


def extract_branch_unitaries(final_states, model, n_fock):
    """Extract 2×2 qubit unitaries for each Fock branch from simulated states.

    Parameters
    ----------
    final_states : list of qt.Qobj
        Final states from simulations with initial states |g,n⟩ and |e,n⟩ for each n.
        Expected order: [|g,0⟩_final, |e,0⟩_final, |g,1⟩_final, |e,1⟩_final, ...]
    model : DispersiveTransmonCavityModel
    n_fock : int
        Number of Fock levels.

    Returns
    -------
    list of 2×2 ndarray
        Effective qubit unitary for each branch n=0..n_fock-1.
    """
    unitaries = []
    for n in range(n_fock):
        # Initial states: |g,n⟩ and |e,n⟩
        # We extract the qubit part conditioned on Fock level n
        idx_g = 2 * n
        idx_e = 2 * n + 1
        psi_g_final = final_states[idx_g]
        psi_e_final = final_states[idx_e]

        # Extract qubit amplitudes conditioned on staying in Fock n
        n_cav = model.n_cav
        # Full state is |qubit⟩ ⊗ |cavity⟩, indices: q*n_cav + cav
        psi_g_arr = psi_g_final.full().flatten()
        psi_e_arr = psi_e_final.full().flatten()

        # Qubit amplitudes in Fock-n sector: indices n (q=0) and n_cav+n (q=1)
        U_n = np.array([
            [psi_g_arr[n], psi_e_arr[n]],             # |g,n⟩ column
            [psi_g_arr[n_cav + n], psi_e_arr[n_cav + n]],  # |e,n⟩ column
        ], dtype=np.complex128)

        unitaries.append(U_n)

    return unitaries


def extract_leakage(final_states, model, n_fock):
    """Leakage probability to the |f⟩ level for each branch.

    Returns array of shape (n_fock,) where entry n is the average leakage
    probability from |g,n⟩ and |e,n⟩ initial states into any |f,m⟩.
    Requires n_tr >= 3.
    """
    n_cav = model.n_cav
    n_tr = model.n_tr
    if n_tr < 3:
        return np.zeros(n_fock)

    leakages = np.zeros(n_fock)
    for n in range(n_fock):
        leak = 0.0
        for k in range(2):  # g, e initial states
            psi = final_states[2 * n + k].full().flatten()
            # |f,m⟩ indices: 2*n_cav + m for m = 0..n_cav-1
            f_pop = np.sum(np.abs(psi[2 * n_cav : 3 * n_cav]) ** 2)
            leak += f_pop
        leakages[n] = leak / 2.0  # average over g, e
    return leakages
