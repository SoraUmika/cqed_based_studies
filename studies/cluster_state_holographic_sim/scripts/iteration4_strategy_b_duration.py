"""Iteration 4 follow-up: duration-penalized bounded Strategy B optimisation.

This targets the user request to minimise active gate time for the
SQR-based decomposition while keeping the displacement bound fixed.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

STUDY_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STUDY_ROOT / "data"
ART_DIR = STUDY_ROOT / "artifacts"
for directory in (DATA_DIR, ART_DIR):
    directory.mkdir(parents=True, exist_ok=True)

SIM_ROOT = Path("C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group/Users/Users_JianJun/cQED_simulation")
if str(SIM_ROOT) not in sys.path:
    sys.path.insert(0, str(SIM_ROOT))

from cqed_sim.unitary_synthesis import (  # noqa: E402
    ConditionalPhaseSQR,
    Displacement,
    DriftPhaseModel,
    ExecutionOptions,
    LeakagePenalty,
    MultiObjective,
    SQR,
    Subspace,
    SynthesisConstraints,
    TargetUnitary,
    UnitarySynthesizer,
    GateSequence,
    simulate_sequence,
    subspace_unitary_fidelity,
)
from cqed_sim.unitary_synthesis.targets import make_target  # noqa: E402

NO_DRIFT = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)
U_TARGET = make_target("cluster", n_match=1)
TARGET = TargetUnitary(U_TARGET, ignore_global_phase=True)


def make_subspace(n_cav: int) -> Subspace:
    return Subspace.custom(2 * n_cav, [0, 1, n_cav, n_cav + 1], ["|g,0>", "|g,1>", "|e,0>", "|e,1>"])


def make_bounded_strategy_b(n_blocks: int, n_cav: int) -> GateSequence:
    gates = []
    theta = [0.0] * n_cav
    theta[0] = np.pi / 2
    if n_cav > 1:
        theta[1] = np.pi / 4
    phi = [0.0] * n_cav
    for block in range(n_blocks):
        gates.append(Displacement(f"D{block}", alpha=0.3 + 0.0j, duration=200e-9))
        gates.append(SQR(f"S{2 * block}", theta_n=theta[:], phi_n=phi[:], drift_model=NO_DRIFT, duration=400e-9))
        gates.append(ConditionalPhaseSQR(f"CP{block}", phases_n=[0.0] * n_cav, drift_model=NO_DRIFT, duration=200e-9))
        gates.append(SQR(f"S{2 * block + 1}", theta_n=theta[:], phi_n=phi[:], drift_model=NO_DRIFT, duration=400e-9))
    gates.append(Displacement(f"D{n_blocks}", alpha=0.3 + 0.0j, duration=200e-9))
    return GateSequence(gates=gates, n_cav=n_cav)


def rescale_sequence(sequence: GateSequence, target_n_cav: int) -> GateSequence:
    new_gates = []
    for gate in sequence.gates:
        if isinstance(gate, SQR):
            theta = list(gate.theta_n) + [0.0] * (target_n_cav - len(gate.theta_n))
            phi = list(gate.phi_n) + [0.0] * (target_n_cav - len(gate.phi_n))
            new_gates.append(
                SQR(
                    gate.name,
                    theta_n=theta[:target_n_cav],
                    phi_n=phi[:target_n_cav],
                    drift_model=gate.drift_model,
                    duration=gate.duration,
                )
            )
        elif isinstance(gate, ConditionalPhaseSQR):
            phases = list(gate.phases_n) + [0.0] * (target_n_cav - len(gate.phases_n))
            new_gates.append(
                ConditionalPhaseSQR(
                    gate.name,
                    phases_n=phases[:target_n_cav],
                    drift_model=gate.drift_model,
                    duration=gate.duration,
                )
            )
        elif isinstance(gate, Displacement):
            new_gates.append(Displacement(gate.name, alpha=gate.alpha, duration=gate.duration))
        else:
            new_gates.append(gate)
    return GateSequence(gates=new_gates, n_cav=target_n_cav)


def main() -> None:
    n_cav_opt = 2
    sequence = make_bounded_strategy_b(n_blocks=2, n_cav=n_cav_opt)
    subspace = make_subspace(n_cav_opt)
    constraints = SynthesisConstraints(max_amplitude=0.3)
    objectives = MultiObjective(fidelity_weight=1.0, leakage_weight=0.05, duration_weight=0.01)

    start = time.perf_counter()
    synthesizer = UnitarySynthesizer(
        primitives=sequence.gates,
        subspace=subspace,
        objectives=objectives,
        leakage_penalty=LeakagePenalty(weight=0.05),
        synthesis_constraints=constraints,
        execution=ExecutionOptions(engine="auto", use_fast_path=True),
    )
    result = synthesizer.fit(target=TARGET, init_guess="heuristic", multistart=3, maxiter=200)
    elapsed = time.perf_counter() - start

    subspace_fidelity = float(subspace_unitary_fidelity(result.simulation.subspace_operator, U_TARGET, gauge="global"))
    total_duration_ns = float(sum(gate.duration for gate in result.sequence.gates) * 1e9)

    seq_nc12 = rescale_sequence(result.sequence, 12)
    sim_nc12 = simulate_sequence(seq_nc12, make_subspace(12))
    nc12_fidelity = float(subspace_unitary_fidelity(sim_nc12.subspace_operator, U_TARGET, gauge="global"))
    full_operator = sim_nc12.full_operator
    full_operator_np = full_operator.full() if hasattr(full_operator, "full") else np.asarray(full_operator)
    logical_indices = [0, 1, 12, 13]
    leakage = 0.0
    for basis_index in logical_indices:
        basis_state = np.zeros(24, dtype=complex)
        basis_state[basis_index] = 1.0
        propagated = full_operator_np @ basis_state
        subspace_population = sum(abs(propagated[index]) ** 2 for index in logical_indices)
        leakage += 1.0 - subspace_population
    avg_leakage = float(leakage / len(logical_indices))

    output = {
        "label": "B2_amp0.3_duration_penalized",
        "success": bool(result.success),
        "objective": float(result.objective),
        "fidelity_nc2": subspace_fidelity,
        "fidelity_nc12": nc12_fidelity,
        "avg_leakage_nc12": avg_leakage,
        "elapsed_s": elapsed,
        "max_amp": 0.3,
        "duration_weight": 0.01,
        "total_duration_ns": total_duration_ns,
        "gate_durations_ns": {gate.name: float(gate.duration * 1e9) for gate in result.sequence.gates},
        "disp_amps": [float(abs(gate.alpha)) for gate in result.sequence.gates if hasattr(gate, "alpha")],
    }

    (DATA_DIR / "iteration4_strategy_b_duration.json").write_text(json.dumps(output, indent=2))
    result.save(str(ART_DIR / "best_strategy_B_duration_penalized.json"), include_history=True)
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()