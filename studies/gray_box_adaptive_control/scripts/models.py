"""
models.py — Physical model definitions for the gray-box adaptive control study.

Defines the truth model (with chi_higher and kerr) and the learner model (without),
as well as helper factories for frame specs, target subspaces, and the target matrix.

Physical system: DispersiveTransmonCavityModel
    n_cav = 4 cavity Fock states, n_tr = 3 transmon levels
    Full Hilbert space dimension: 4 * 3 = 12
    Control subspace: {|g,0>, |g,1>, |g,2>, |g,3>, |e,0>, |e,1>, |e,2>, |e,3>} (8-dim)

Qubit-major basis ordering (repository convention):
    Index = qubit_level * n_cav + cavity_level
    |g,0>=0, |g,1>=1, |g,2>=2, |g,3>=3
    |e,0>=4, |e,1>=5, |e,2>=6, |e,3>=7
    |f,0>=8, |f,1>=9, |f,2>=10, |f,3>=11

Target gate: Qubit X gate simultaneously for all Fock sectors:
    |g,n> -> |e,n> and |e,n> -> |g,n> for n=0,1,2,3
    In the 8-dim subspace this is the swap block [[0,I],[I,0]].
"""

import numpy as np

from cqed_sim import DispersiveTransmonCavityModel, FrameSpec, NoiseSpec
from cqed_sim.unitary_synthesis import Subspace

# ---------------------------------------------------------------------------
# Physical constants (rad/s)
# ---------------------------------------------------------------------------

OMEGA_C = 2 * np.pi * 5.241e9       # Storage cavity angular frequency (rad/s)
OMEGA_Q = 2 * np.pi * 6.150e9       # Qubit angular frequency (rad/s)
ALPHA = 2 * np.pi * (-255e6)        # Transmon anharmonicity (rad/s)

# True model dispersive coupling parameters
CHI_TRUE = 2 * np.pi * (-2.84e6)     # Qubit-cavity dispersive coupling chi (rad/s)
CHI_HIGHER_TRUE = 2 * np.pi * (-21e3)  # Second-order chi: (chi/2)*n*(n-1) term (rad/s)
KERR_TRUE = 2 * np.pi * (-28e3)       # Self-Kerr of cavity (rad/s)

# Learner's prior (30% off chi_true, no knowledge of chi_higher or kerr)
CHI_PRIOR = 2 * np.pi * (-2.0e6)    # Learner's prior chi estimate (rad/s)

# ---------------------------------------------------------------------------
# Hilbert space dimensions
# ---------------------------------------------------------------------------

N_CAV = 4   # Cavity Fock truncation (states 0..3)
N_TR = 3    # Transmon level truncation (levels g, e, f)

# Full Hilbert space dimension
FULL_DIM = N_CAV * N_TR  # = 12

# ---------------------------------------------------------------------------
# Noise and readout parameters
# ---------------------------------------------------------------------------

# Noisy evaluation noise spec (Lindblad terms)
TRUTH_NOISE = NoiseSpec(
    t1=50e-6,         # T1 = 50 us (qubit energy relaxation)
    tphi=40e-6,       # T_phi = 40 us (qubit pure dephasing)
    kappa=1.0 / 50e-6,  # Cavity loss rate kappa = 1/T1_cav (rad/s)
)

# Readout confusion matrix: M[i,j] = P(observe i | true j)
# Row = observed state, Col = true state
# M @ [P_g, P_e] = [P_obs_g, P_obs_e]
CONFUSION_MATRIX = np.array([
    [0.97, 0.05],   # P(observe g | true g), P(observe g | true e)
    [0.03, 0.95],   # P(observe e | true g), P(observe e | true e)
])

# Default probe shot count for chi Ramsey experiment
N_SHOTS_PROBE = 1000

# ---------------------------------------------------------------------------
# Model factories
# ---------------------------------------------------------------------------


def make_truth_model() -> DispersiveTransmonCavityModel:
    """
    Build the truth model including chi_higher and Kerr nonlinearity.

    The truth model is used to:
    - Generate probe data (chi Ramsey measurement)
    - Evaluate all GRAPE pulses

    Returns
    -------
    DispersiveTransmonCavityModel
        Full truth model with chi, chi_higher, and kerr set to true values.
    """
    return DispersiveTransmonCavityModel(
        omega_c=OMEGA_C,
        omega_q=OMEGA_Q,
        alpha=ALPHA,
        chi=CHI_TRUE,
        chi_higher=(CHI_HIGHER_TRUE,),  # sequence of higher-order chi coefficients
        kerr=KERR_TRUE,
        n_cav=N_CAV,
        n_tr=N_TR,
    )


def make_learner_model(
    chi: float = CHI_PRIOR,
    chi_higher_val: float = 0.0,
    kerr_val: float = 0.0,
) -> DispersiveTransmonCavityModel:
    """
    Build the learner's model with a specified chi value.

    The learner model is used for GRAPE optimization. By default it omits
    chi_higher and kerr (the learner does not know these). After gray-box
    inference, chi is updated to chi_hat.

    Parameters
    ----------
    chi : float
        Dispersive coupling to use (rad/s). Default is CHI_PRIOR.
    chi_higher_val : float
        Second-order chi coefficient (rad/s). Default is 0.0 (unknown to learner).
    kerr_val : float
        Cavity self-Kerr (rad/s). Default is 0.0 (unknown to learner).

    Returns
    -------
    DispersiveTransmonCavityModel
        Learner's model with specified parameters.
    """
    chi_higher_seq = (chi_higher_val,) if chi_higher_val != 0.0 else ()
    return DispersiveTransmonCavityModel(
        omega_c=OMEGA_C,
        omega_q=OMEGA_Q,
        alpha=ALPHA,
        chi=float(chi),
        chi_higher=chi_higher_seq,
        kerr=float(kerr_val),
        n_cav=N_CAV,
        n_tr=N_TR,
    )


def make_frame(model: DispersiveTransmonCavityModel) -> FrameSpec:
    """
    Build a FrameSpec rotating at the model's own cavity and qubit frequencies.

    Both truth model and learner model frames use the model's own omega_c and
    omega_q. Since we only vary chi (not omega_c or omega_q), both frames are
    identical in this study.

    Parameters
    ----------
    model : DispersiveTransmonCavityModel
        Model whose frequencies define the rotating frame.

    Returns
    -------
    FrameSpec
    """
    return FrameSpec(
        omega_c_frame=float(model.omega_c),
        omega_q_frame=float(model.omega_q),
    )


def make_grape_subspace() -> Subspace:
    """
    Build the control subspace for GRAPE: the g+e qubit manifold over Fock 0-3.

    In the qubit-major basis ordering (q * n_cav + n):
        |g,0>=0, |g,1>=1, |g,2>=2, |g,3>=3
        |e,0>=4, |e,1>=5, |e,2>=6, |e,3>=7

    The f-level states (indices 8-11) are leakage states.

    Returns
    -------
    Subspace
        Custom subspace with full_dim=12, indices for g+e manifold.
    """
    # g manifold: indices 0,1,2,3  (g * N_CAV + n = 0*4 + n)
    # e manifold: indices 4,5,6,7  (e * N_CAV + n = 1*4 + n)
    indices = list(range(N_CAV)) + list(range(N_CAV, 2 * N_CAV))
    labels = (
        [f"|g,{n}>" for n in range(N_CAV)]
        + [f"|e,{n}>" for n in range(N_CAV)]
    )
    return Subspace.custom(
        full_dim=FULL_DIM,
        indices=indices,
        labels=labels,
    )


def make_target_matrix() -> np.ndarray:
    """
    Build the 8x8 target unitary matrix for the qubit X gate on all Fock sectors.

    The target is the simultaneous qubit flip:
        |g,n> -> |e,n>,  |e,n> -> |g,n>   for n = 0, 1, 2, 3

    In the subspace ordering [|g,0>, |g,1>, |g,2>, |g,3>, |e,0>, |e,1>, |e,2>, |e,3>]:

        U_target = [[0, I4],
                    [I4, 0]]

    where I4 is the 4x4 identity. This is block-off-diagonal (pure swap between
    g and e sectors).

    Returns
    -------
    np.ndarray, shape (8, 8), dtype complex128
    """
    I4 = np.eye(N_CAV, dtype=np.complex128)
    Z4 = np.zeros((N_CAV, N_CAV), dtype=np.complex128)
    return np.block([[Z4, I4], [I4, Z4]])


# ---------------------------------------------------------------------------
# Index lookup utilities
# ---------------------------------------------------------------------------


def basis_index(qubit_level: int, cavity_level: int) -> int:
    """
    Return the flat basis index for |q, n> in the 12-dim qubit-major ordering.

    Parameters
    ----------
    qubit_level : int
        Qubit level (0=g, 1=e, 2=f).
    cavity_level : int
        Cavity Fock level (0..N_CAV-1).

    Returns
    -------
    int
    """
    return int(qubit_level) * N_CAV + int(cavity_level)


def g_state_index(n: int) -> int:
    """Return flat index for |g, n>."""
    return basis_index(0, n)


def e_state_index(n: int) -> int:
    """Return flat index for |e, n>."""
    return basis_index(1, n)


if __name__ == "__main__":
    # Quick sanity check
    truth = make_truth_model()
    learner = make_learner_model()
    frame_t = make_frame(truth)
    frame_l = make_frame(learner)
    sub = make_grape_subspace()
    target = make_target_matrix()

    print("Truth model:")
    print(f"  omega_c = {truth.omega_c / (2*np.pi) / 1e9:.4f} GHz")
    print(f"  omega_q = {truth.omega_q / (2*np.pi) / 1e9:.4f} GHz")
    print(f"  chi     = {truth.chi / (2*np.pi) / 1e6:.4f} MHz")
    print(f"  chi_higher = {list(truth.chi_higher)}")
    print(f"  kerr    = {truth.kerr / (2*np.pi) / 1e3:.2f} kHz")
    print(f"  subsystem_dims = {truth.subsystem_dims}")
    print()
    print("Learner model:")
    print(f"  chi     = {learner.chi / (2*np.pi) / 1e6:.4f} MHz  (prior)")
    print()
    print("Subspace:")
    print(f"  indices = {sub.indices}")
    print(f"  labels  = {sub.labels}")
    print()
    print("Target matrix (8x8):")
    print(np.round(np.abs(target), 2))
    print()
    print(f"Basis indices: |g,2>={g_state_index(2)}, |e,2>={e_state_index(2)}")
