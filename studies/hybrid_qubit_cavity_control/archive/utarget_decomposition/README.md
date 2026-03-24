# cQED Synthesis of a Hybrid Qubit-Cavity Unitary

## Problem Class
`OPT` | `DES` | `ANA`

## Motivation

We study whether the target qubit-cavity unitary

```
U_target = (1/√2) [[1,0,1,0],[1,0,-1,0],[0,1,0,1],[0,-1,0,1]]
```

in the logical subspace {|g,0⟩, |g,1⟩, |e,0⟩, |e,1⟩} can be synthesized
efficiently in a dispersive cQED system using physical primitives (SQR,
conditional phase, qubit rotations, displacements, SNAP).

The target equals `(I_q⊗H_c) · CNOT_{c→q} · CNOT_{q→c}` — a three-CNOT-equivalent
circuit mixing qubit and cavity logical degrees of freedom.  This study
determines the **best physically meaningful synthesis strategy** and assesses
feasibility at device-realistic parameters.

## Goals

1. Verify algebraic structure of U_target and the CNOT factorization.
2. Determine whether U_target is exactly or only approximately synthesizable
   on the logical subspace with available primitive families.
3. Find the shortest primitive sequence for each of Libraries A, B, C.
4. Quantify leakage outside {|0⟩,|1⟩} cavity, sensitivity to χ, K, χ′, and
   pulse amplitude errors.
5. Deliver projected process fidelity F_proj and leakage L for the best
   sequence found.

## Methods

### Device parameters (rad/s unless stated)
- χ  = 2π × (−2.84 MHz)
- χ′ = 2π × (−21 kHz)   [chi_higher[0] in DispersiveTransmonCavityModel]
- K  = 2π × (−28 kHz)   [kerr in DispersiveTransmonCavityModel]
- Transmon dim n_tr = 2 (ge only; f-level leakage excluded in this study)

### cqed_sim modules used
- `cqed_sim.core.DispersiveTransmonCavityModel`, `FrameSpec`
- `cqed_sim.gates.coupled` (SQR, dispersive_phase, conditional_displacement)
- `cqed_sim.gates.bosonic` (displacement, snap, oscillator_rotation)
- `cqed_sim.unitary_synthesis` (UnitarySynthesizer, GateSequence, Subspace,
  QubitRotation, SQR, Displacement, FreeEvolveCondPhase, ConditionalPhaseSQR,
  SNAP, DriftPhaseModel, TargetUnitary, LeakagePenalty, MultiObjective)
- `cqed_sim.unitary_synthesis.metrics` (subspace_unitary_fidelity,
  leakage_metrics)

### Primitive libraries tested
- **Library A**: {D(α), R_q(θ,φ), SQR_n, FreeEvolveCondPhase}
- **Library B**: {D(α), R_q(θ,φ), SQR_n, ConditionalPhaseSQR}
- **Library C**: {D(α), R_q(θ,φ), SQR_n, SNAP}  ← controllability probe

### Convergence / sanity checks
- Unitarity of target verified analytically (Phase 1)
- Logical subspace projector fidelity checked vs ideal gate products (Phase 1)
- Truncation convergence checked for N ∈ {4, 6, 8, 10} (Phase 3)
- Robustness sweep: ±5% χ, ±5% amplitude, ±10% pulse duration (Phase 6)

### Sign/frame conventions
- **Frame**: rotating at (ω_c, ω_q), so static Hamiltonian diagonal terms vanish
  and only dispersive/Kerr residuals remain.
- **H/ℏ** = −(χ/2) a†a σ_z − (χ′/2) a†²a² σ_z + (K/2) a†²a²
  In library coordinates: DriftPhaseModel(chi=χ, chi2=χ′, kerr=K) with the
  excitation-projector convention chi * n * |e⟩⟨e|.
- **σ_z** convention: σ_z|g⟩ = +|g⟩, σ_z|e⟩ = −|e⟩ (QuTiP default).
- **Basis ordering**: qubit first, cavity second: flat index = q * n_cav + n.
  Target matrix rows/cols ordered as {|g,0⟩,|g,1⟩,|e,0⟩,|e,1⟩}.
- **SQR** (cqed_sim.gates.coupled): cavity-first ordering |n⟩⟨n| ⊗ R(θ,φ).
  Note: `sqr_op` in `cqed_sim.core.ideal_gates` uses qubit-first ordering.
  We use `cqed_sim.gates.coupled.sqr` throughout.

## Expected Outcomes

- Phase 1: unitarity verified, CNOT factorization confirmed algebraically.
- Phase 2: short 3–5 gate logical decomposition found in 4D ideal space.
- Phase 3: cavity truncation N ≥ 6 sufficient for ≤ 0.1% convergence error.
- Phase 4: F_proj > 0.99 achievable with Library B or C at N = 8.
- Phase 5: finite pulse duration (SQR ~2 μs, displacement ~200 ns) causes
  ≲1% infidelity from χ′/K accumulation during pulses.
- Phase 6: ±5% χ calibration error causes ≲2% fidelity loss; robustness
  ordering Library C > Library B > Library A.

## Scripts

| Script | Phase | Description |
|---|---|---|
| `phase1_algebraic.py` | 1 | Verify unitarity, CNOT factorization, subproblem gates |
| `phase2_idealized.py` | 2 | 4D ideal synthesis — seed sequences for Phase 4 |
| `phase3_embedding.py` | 3 | Convergence vs cavity truncation N |
| `phase4_variational.py` | 4 | Full variational synthesis (main result) |
| `phase5_pulses.py` | 5 | Finite-duration pulse effects |
| `phase6_robustness.py` | 6 | Robustness analysis |

## Status

COMPLETE

## Results Summary

| Phase | Key Result |
|---|---|
| 1 (Algebraic) | U_target maximally entangling; CNOT factorization exact (38/38 checks pass) |
| 2 (Ideal 4D) | D+SQR+CP achieves F=1.0 at depth-2; SQR+CP alone limited to F=0.5 |
| 3 (Embedding N=8) | Truncation converges at N=6; analytic decomposition F=1.0 for all N |
| 4 (Physical) | **L1c (D+Rq+SQR, depth-11): F=0.9193, L_avg=0.101** (best library) |
| 5 (Pulses) | Drift fully compensated by optimizer; residual infidelity = leakage only |
| 6 (Robustness) | Zero sensitivity to χ and duration; 0.58% RMS from amplitude; 0.12% from phase |

See `report/report.pdf` for the full 9-page report.

## Suggested Upstreaming

- The `make_target` convenience function in `cqed_sim.unitary_synthesis.targets`
  could be extended with a `"cnot_factored"` entry for the specific U_target
  studied here.
- The logical CNOT_{q→c} gate (qubit-controlled cavity logical NOT) is not
  currently in the built-in target library; adding it would help future hybrid
  gate synthesis studies.
