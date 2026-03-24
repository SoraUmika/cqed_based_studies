# Pulse-Waveform Design for Selective Qubit Rotation (SQR) in Dispersive cQED

## Problem Class

OPT | ANA | DES

## Motivation

Selective qubit rotation (SQR) is a fundamental control primitive in dispersive cQED
that acts on one cavity Fock branch while leaving others unchanged. However, the naive
single-tone selective pulse achieves at best a conditional-phase SQR rather than a true
SQR, because off-resonant spectator branches accumulate unavoidable AC Stark phase.

This study systematically compares waveform families (single-tone, multitone,
composite/echoed, and numerically optimized) against two target types (true SQR and
conditional-phase SQR) to determine:

- the minimum χT required for useful branch-selective control,
- whether ideal SQR is the right primitive or conditional-phase SQR suffices,
- which waveform family offers the best fidelity-duration tradeoff on a truncated subspace.

## Goals

1. Establish baseline performance of a single-tone Gaussian selective pulse for true SQR
   and conditional-phase SQR targets on a truncated subspace n=0…3 (N=4 Fock levels).
2. Quantify branch-resolved process fidelity and spectator phase as a function of χT
   for χT ∈ [1, 50].
3. Compare four waveform families at equal χT: (a) single-tone Gaussian, (b) single-tone
   phase-modulated Gaussian, (c) one-segment multitone, (d) two-segment multitone echo.
4. Determine whether conditional-phase SQR fidelity is substantially higher than true SQR
   fidelity across all families.
5. Use piecewise-optimized control (GRAPE) as a theoretical upper bound on achievable
   fidelity for the truncated subspace.
6. Identify the minimum χT needed to reach F ≥ 0.999 for each waveform family and target type.
7. Test robustness of the best pulses to χ′, cavity Kerr K, amplitude error (±5%), and
   detuning error (±χ/10).
8. Determine whether composite/echoed pulses are necessary to approximate true SQR.
9. Extend the study from the fixed $R_X(\pi)$ target to generalized rotations with
   $\theta \in \{\pi/2, \pi\}$ and $\phi \in \{0, \pi/4, \pi/2\}$.
10. Quantify logical-operator performance for target branches $n_0 \in \{0, 1, 2\}$ and
   logical truncations $N \in \{3, 4\}$.
11. Separate strict logical fidelity from cavity-block-phase-relaxed fidelity to identify
   whether errors are removable block phases or genuine logical mismatch.
12. Compare echoed and non-echoed selective pulses on the same generalized logical metrics.

## Methods

### cqed_sim modules used

- `DispersiveTransmonCavityModel` — system model with χ, χ′, Kerr
- `FrameSpec` — rotating frame at qubit and cavity frequencies
- `Pulse`, `GaussianEnvelope`, `multitone_gaussian_envelope`, `MultitoneTone` — pulse construction
- `build_sqr_multitone_pulse`, `build_sqr_tone_specs` — standard SQR multitone builder
- `SequenceCompiler`, `CompiledSequence` — waveform compilation
- `simulate_sequence`, `SimulationConfig`, `SimulationSession`, `prepare_simulation` — time-domain simulation
- `conditioned_bloch_xyz`, `conditioned_qubit_state` — branch-resolved state extraction
- `manifold_transition_frequency` — branch-resolved qubit frequencies
- `sqr` (gates.coupled) — ideal SQR target operator
- `multi_sqr` — multi-branch SQR target
- `subspace_unitary_fidelity`, `leakage_metrics` — fidelity metrics
- `calibrate_sqr_gate`, `extract_effective_qubit_unitary` — SQR calibration
- `ConditionedQubitTargets`, `run_conditioned_multitone_validation` — conditioned multitone layer
- `run_targeted_subspace_multitone_validation` — full logical subspace validation
- `GrapeSolver`, `build_control_problem_from_model`, `UnitaryObjective` — GRAPE optimal control
- `logical_block_phase_diagnostics`, `logical_block_phase_op` — cavity-only logical block-phase extraction and compiled correction layers

### Standalone extensions

- Custom composite two-segment pulse builder (not in cqed_sim)
- Custom phase-modulated envelope function
- Custom conditional-phase SQR fidelity metric (closest Z-rotation on spectators)

### Generalized-target extension

- `scripts/run_extended_targets_and_echo.py` performs the closed-system extension scan over
   generalized $(\theta, \phi)$, target-branch, and logical-truncation cases.
- Extended waveform families: single-tone Gaussian, one-segment multitone baseline,
   echoed single-tone Gaussian, echoed one-segment multitone.
- Joint logical metrics: strict logical fidelity, cavity-block-phase-relaxed fidelity,
   same-block population, state-transfer fidelity, and leakage.
- Extended figures: `extended_representative_family_comparison`,
   `extended_axis_angle_heatmaps`, `extended_branch_truncation_sensitivity`, and
   `extended_cavity_phase_effects`.

## Assumptions

- Transmon: ω_q = 2π × 6.150 GHz, α = 2π × (−255 MHz), n_tr = 3 (g/e/f manifold)
- Storage: ω_c = 2π × 5.241 GHz
- Dispersive shift: χ = 2π × (−2.84 MHz)
- Higher-order: χ′ = 2π × (−21 kHz), K = 2π × (−28 kHz) [Phase 4]
- Decoherence: T1 = T2 = 20 μs [Phase 5]
- Truncated Fock space: N = 4 (n = 0, 1, 2, 3)
- Hilbert space: n_cav = 6, n_tr = 3 → dim = 18
- Convergence: fidelity stable to 1×10⁻⁵ when doubling n_cav
- Gate duration scan: χT/(2π) ∈ {0.5, 1, 1.5, 2, 3, 5, 7, 10} (true dimensionless parameter)
- Simulation dt = 2 ns
- Leakage to |f⟩ tracked at each χT/(2π) value
- Phase-compilation follow-up: enlarged logical window N = 8 implemented with n_cav = 10
- Compiled comparisons allow one global qubit virtual-Z in addition to any fitted cavity-only phase layer

## Expected Outcomes

- Single-tone Gaussian achieves high conditional-phase SQR fidelity (F > 0.999)
  at χT/(2π) ≥ 3 but poor true SQR fidelity due to spectator phase.
- Cosine-squared (Hann) envelope offers the best selectivity at intermediate χT/(2π).
- GRAPE with per-branch phase freedom achieves F ≈ 1 (proper upper bound).
- Conditional-phase SQR is significantly easier than true SQR for all waveform families.
- The physically natural primitive is conditional-phase SQR.
- Optimal χT/(2π) ≈ 2–3 balancing selectivity and decoherence (F_net ≈ 0.98).
- Leakage to |f⟩ negligible (< 10⁻⁴) at all drive amplitudes used.

## Status

COMPLETE — Original five-phase study, generalized-target extension, optimized
single-branch multitone study, all-branch simultaneous multitone study, and
the phase-compilation follow-up are complete and merged into the report.

---

## Study Extensions

### Optimized Single-Branch Multitone SQR Gates

Three optimized multitone families were studied (independent-tone 8p, detuned 12p,
smooth-basis 16p) targeting R_X(π) on branch n₀ = 1. Key finding: optimized multitone
improves F_block by up to 5.6× at χT/2π = 0.5 but offers ≤ 1% improvement for
χT/2π ≥ 3. Results merged into the main report.

**Files:** `scripts/run_followup_multitone.py`, `data/followup_multitone_results.npz`,
`figures/followup_*.{png,pdf}`

### All-Branch Simultaneous Multitone SQR Study

Targets R_X(π) on ALL Fock branches n = 0…3 simultaneously, using a 4-tone multitone
pulse with independent amplitudes, phases, and detunings (12 params). Optimization via
differential evolution (DE) global search + Nelder-Mead local refinement. GRAPE provides
the optimal control upper bound. χT/2π scanned over {0.5, 1.0, 1.5, 2.0, 3.0, 5.0}.

**Key result:** DE+NM achieves F_block ≥ 0.987 (per-branch rotations accurate) but
F_strict ≤ 0.27 (inter-branch phases incoherent). GRAPE achieves F ≥ 0.9999 at all
χT/2π, confirming the limitation is the parametric Gaussian ansatz, not the system.

| χT/2π | Baseline F_b | DE+NM F_b | DE+NM F_s | GRAPE F |
|-------|-------------|-----------|-----------|---------|
| 0.5   | 0.529       | 0.987     | 0.023     | 0.9999  |
| 1.0   | 0.148       | 0.994     | 0.034     | 1.0000  |
| 1.5   | 0.820       | 0.994     | 0.036     | 1.0000  |
| 2.0   | 0.480       | 1.000     | 0.266     | 1.0000  |
| 3.0   | 0.986       | 1.000     | 0.022     | 1.0000  |
| 5.0   | 0.983       | 1.000     | 0.228     | 1.0000  |

**Files:** `scripts/run_allbranch_multitone.py`, `scripts/run_allbranch_fast.py`,
`scripts/plot_allbranch.py`, `data/allbranch_multitone_results.npz`,
`figures/allbranch_*.{png,pdf}`

### Phase-Compilation Follow-Up

This extension explicitly tests the hypothesis that the remaining strict-logical
error is mostly a correctable cavity-only Fock-phase profile. Two complementary
cases were analyzed:

1. Enlarged-window single-target SQR on logical levels n = 0…7 for Gaussian and
   cosine-squared baselines.
2. Representative short-gate structured all-branch multitone pulses obtained
   from the existing Phase-B Nelder-Mead optimization.

**Key findings:**

- For single-target SQR, the extracted cavity-phase profile is extremely smooth
  and nearly linear in n (worst linear-fit RMS < 4.6 × 10⁻⁵ rad across the scan).
- On the enlarged N = 8 window, linear cavity-phase compilation improves the
  global-Z-gauged strict logical fidelity by 0.017–0.023 for both Gaussian and
  cosine-squared envelopes, and improves coherent-superposition benchmarks.
- The compiled cavity phase does not close the full gap to the branch-local-Z
  relaxed fidelity; spectator branch-local qubit-Z structure remains the larger
  error source.
- In the all-branch short-gate structured multitone case, cavity-only
  compilation fails dramatically: at χT/2π = 0.5, the raw strict fidelity is
  0.138, cavity-compiled fidelity is 0.028, branch-local-Z-relaxed fidelity is
  0.972, and the saved GRAPE reference is 0.9999. That regime is therefore not
  a bosonic-phase-compilation problem.

**Files:** `scripts/phase_compilation_common.py`,
`scripts/run_phase_compilation_followup.py`,
`scripts/plot_phase_compilation_followup.py`,
`scripts/validate_phase_compilation_followup.py`,
`data/phase_compilation_results.npz`,
`data/phase_compilation_summary.json`,
`figures/phase_compilation_*.{png,pdf}`

## Suggested Upstreaming

- Promote the reusable logical phase-compilation helpers from
  `scripts/phase_compilation_common.py` into a `cqed_sim`-level analysis utility
  alongside the targeted-subspace workflow.
- Expose a library helper for “global qubit virtual-Z + cavity-only block-phase”
  compiled-fidelity evaluation, since this gauge is the natural experiment-facing
  comparison used in this follow-up.
