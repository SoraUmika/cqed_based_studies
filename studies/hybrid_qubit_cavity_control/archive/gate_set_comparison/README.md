# Hybrid Universal Control Gate-Set Comparison in a 2x2 Qubit+Cavity Logical Subspace

## Problem Class
`OPT` | `DES` | `ANA`

## Motivation
This study asks a practical control-design question for the local cQED platform:
which gate family is the best realistic route to hybrid qubit+cavity universal
control when we start from the smallest nontrivial logical space
`{|g,0>, |g,1>, |e,0>, |e,1>}`.

The point is not to restate abstract universality. The point is to compare
candidate control libraries under the device Hamiltonian and timing hierarchy
already used in this repo: fast qubit rotations, fast cavity displacements,
slower selective dispersive operations, native dispersive drift, native
exchange/sideband-style interactions, and direct waveform optimization through
the `cqed_sim` GRAPE layer.

The expanded benchmark is still deliberately narrow:

1. Fock logical encoding `|0>_L = |0>`, `|1>_L = |1>`.
2. A 4D logical subspace embedded in a larger truncated Hilbert space.
3. At least one maximally entangling hybrid unitary as the main target.
4. Explicit comparison between selective libraries, native-gate libraries,
   and waveform-level GRAPE references.
5. Gate-library ranking by fidelity, leakage, duration, robustness, and
   implementation burden.

## Goals
1. Audit the local `cqed_sim` framework for hybrid gates, logical-subspace
   tools, unitary synthesis, native exchange primitives, and GRAPE support.
2. Write a compact literature-backed design memo covering SNAP-based control,
   SQR-like selective hybrid control, ECD-like conditional-displacement
   control, SWAP-/sideband-/exchange-style native control, CNOD-like fast
   conditional displacement, direct optimal control, and native activated
   nonlinear routes.
3. Extend the local `cqed_sim` copy where needed so native hybrid primitives
   can be benchmarked through the synthesis stack instead of only discussed
   conceptually.
4. Build a reusable benchmark for logical 2x2 qubit+cavity unitary synthesis
   on the Fock-encoded subspace `{|0>, |1>}`.
5. Compare six candidate routes:
   `A = {R_q, D, SNAP}` with native chi-wait allowed,
   `B = {R_q, D, SQR-like selective hybrid control}`,
   `C = {R_q, D, ConditionalDisplacement / ECD-like control}`,
   `D = {R_q, D, native chi-wait entangler}`,
   `E = GRAPE waveform optimization`,
   `F = {R_q, JaynesCummingsExchange, BlueSidebandExchange}`.
6. Use at least one maximally entangling hybrid target and report logical
   fidelity, leakage, duration, robustness, and resource counts.
7. Produce a final recommendation for the best practical route on this device,
   not just the most expressive route in principle.

## Methods
### Local device and model assumptions
- `omega_q / 2pi = 6.150 GHz`
- `omega_c / 2pi = 5.241 GHz`
- `alpha / 2pi = -255 MHz`
- `chi / 2pi = -2.84 MHz`
- `chi' / 2pi = -21 kHz`
- `K / 2pi = -28 kHz`
- Initial logical encoding: cavity Fock `|0>`, `|1>`
- Initial qubit model: two-level `g/e` logical manifold unless a script
  explicitly upgrades to `n_tr = 3` for validation

### `cqed_sim` modules/functions used
- `DispersiveTransmonCavityModel`, `FrameSpec`
- `cqed_sim.gates`: `displacement`, `snap`, `controlled_snap`,
  `conditional_displacement`, `dispersive_phase`, `sqr`, `jaynes_cummings`,
  `blue_sideband`
- `cqed_sim.unitary_synthesis`: `Subspace`, `TargetUnitary`,
  `UnitarySynthesizer`, `GateSequence`, `QubitRotation`, `SNAP`,
  `Displacement`, `SQR`, `ConditionalPhaseSQR`, `FreeEvolveCondPhase`,
  `ConditionalDisplacement`, `JaynesCummingsExchange`,
  `BlueSidebandExchange`, `LeakagePenalty`, `MultiObjective`,
  `ExecutionOptions`, `subspace_unitary_fidelity`, `leakage_metrics`
- `cqed_sim.optimal_control`: `GrapeSolver`, `GrapeConfig`,
  `ModelControlChannelSpec`, `PiecewiseConstantTimeGrid`,
  `UnitaryObjective`, `build_control_problem_from_model`,
  `evaluate_control_with_simulator`
- `SequenceCompiler`, `simulate_sequence`, `NoiseSpec` for replay/validation

### Framework gap audit
1. The local `cqed_sim` copy has been extended during this study to expose
   first-class native synthesis primitives:
   `ConditionalDisplacement`, `JaynesCummingsExchange`, and
   `BlueSidebandExchange`. They are wired into the synthesis API,
   fast-evaluation path, API reference, and dedicated unit tests.
2. The remaining native-control gap is not operator availability but calibrated
   pulse export. The new native primitives are ideal effective-interaction
   models inside `unitary_synthesis`; they are not yet backed by a full
   waveform bridge or experiment-ready calibration stack.
3. No study-ready cat-encoding workflow is currently needed for the first pass,
   so the initial benchmark stays in the Fock `{|0>, |1>}` encoding.
4. The study uses `cqed_sim` GRAPE as the waveform-level reference, but the
   decomposition libraries are compared first through the `unitary_synthesis`
   layer with explicit duration priors.

### Candidate gate libraries
- **Gate Set A**: baseline dispersive bosonic control
  `{QubitRotation, Displacement, SNAP}` with native chi-wait allowed where
  physically natural for the entangler
- **Gate Set B**: hybrid selective control
  `{QubitRotation, Displacement, SQR / ConditionalPhaseSQR}`
- **Gate Set C**: ECD-like fast hybrid control
  `{QubitRotation, Displacement, ConditionalDisplacement}`
- **Gate Set D**: minimal native entangler library
  `{QubitRotation, Displacement, FreeEvolveCondPhase}`
- **Gate Set E**: direct waveform control
  piecewise-constant qubit and storage drives under the full modeled Hamiltonian
- **Gate Set F**: native SWAP-/sideband-style control
  `{QubitRotation, JaynesCummingsExchange, BlueSidebandExchange}`

### Main targets
1. Hybrid logical `CZ(q,c)` on `{|g,0>, |g,1>, |e,0>, |e,1>}`
2. Hybrid logical `CX(c->q)` obtained either directly or by local-basis
   conversion from `CZ`
3. Local cavity basis change `I_q \otimes H_c`

### Evaluation metrics
- Logical unitary fidelity on the selected 4D subspace
- Average and worst-case leakage outside the logical subspace
- Total duration from assigned primitive durations or GRAPE schedule duration
- Robustness under small perturbations in `chi`, `chi'`, `K`, and control scale
- Resource counts: primitive count, optimized parameters, waveform channels
- Composite score `J = F_logical - lambda_leak * L - lambda_T * T/T_ref`

### Validation plan
- Sanity: exact algebraic checks for target operators, native entangler
  equivalence, and native-vs-selective ordering checks
- Convergence: stability versus cavity truncation and GRAPE time-grid resolution
- Literature comparison: qualitative and quantitative comparison to published
  speed/selectivity claims for SNAP, ECD-like control, sideband/SWAP-style
  native control, and direct optimal control

## Expected Outcomes
- A native dispersive entangler plus fast local controls should remain a strong
  practical baseline for 2x2 Fock-encoded hybrid entanglers.
- SNAP-based synthesis should remain expressive but likely slower due to
  selective operations and displacement-induced leakage.
- SQR-like selective hybrid control may reduce circuit depth but may be
  duration-limited relative to fast native entanglers.
- ECD-like conditional displacements should become a fairer benchmark once
  treated as first-class synthesis primitives, but the strict Fock benchmark may
  still penalize them.
- SWAP-/sideband-style native control should be competitive on speed and
  hybridity, especially for transfer-like tasks, while risking more leakage than
  exact dispersive logical control.
- GRAPE should provide an upper benchmark and may reveal whether the best
  decomposition library is already near-optimal.

## Known Limitations
- The benchmark is still restricted to the Fock `{|0>, |1>}` logical encoding.
- `ConditionalDisplacement`, `JaynesCummingsExchange`, and
  `BlueSidebandExchange` are ideal effective-interaction primitives, not yet
  pulse-exported calibrated schedules.
- The first decomposition comparison still emphasizes closed-system logical
  synthesis with robustness sweeps; richer open-system replay is a validation
  step, not the optimization loop itself.
- The first pass keeps the transmon as an effective two-level ancilla for the
  decomposition study, which omits direct `|f>` leakage in the synthesis loop.
- The GRAPE benchmark required best-of-seed fine-grid validation to stabilize
  its time-grid comparison, so the GRAPE numbers should still be treated as a
  strong benchmark rather than a fully exhaustive optimum.

## Results Summary
- Best practical route for the present 2x2 Fock benchmark:
  `{R_q, D, SNAP, native chi-wait}`.
- Best local cavity action: Gate Set A (`D-SNAP-D`) with
  `F_strict = 0.9887`, leakage `0.0185`, duration `1260 ns`.
- Best native local alternative: Gate Set F (blue-sideband + JC exchange) with
  `F_strict = 0.8784`, leakage `0.0788`, duration `440 ns`; it beats the
  selective SQR route on both fidelity and duration but not the SNAP baseline.
- Best exact entangler: native dispersive wait plus fast qubit rotations
  (Gate Set D, also recovered inside Gate Set A) with
  `F_strict = 0.9999999`, negligible leakage, duration `256 ns`.
- Best GRAPE local benchmark: `F_strict = 0.9618` at `320 ns`.
- Best GRAPE entangler benchmark: `F_strict = 0.9458` at `400 ns`.
- Selective SQR entangler reached block-gauge fidelity `1.0` but strict
  fidelity `0.7071`, indicating a missing local cavity-phase correction rather
  than missing entangling structure.
- The first-class ECD-like conditional-displacement primitive was still not
  competitive in the strict `{|0>, |1>}` encoding and should be revisited after
  changing the cavity encoding or adding calibrated echoed-pulse compilation.
- Native SWAP-/sideband-style control is now benchmarked explicitly and emerges
  as a serious fast candidate, but in this smallest Fock block it remains
  leakier and less exact than the native chi-wait entangler and the SNAP local
  cavity gate.

## Validation Summary
- Sanity checks passed, including the new native-vs-selective SWAP benchmarks.
- Truncation convergence passed for the winning decomposition libraries at
  `n_cav = 6, 8, 10`.
- Higher-budget optimizer reruns reproduced the winning decomposition results.
- GRAPE time-grid validation passed in a best-of-seed fine-grid regime.
- See `data/validation_summary.json`, `data/study_summary.json`, and
  `report/report.pdf` for the full validation record.

## Suggested Upstreaming
- Completed locally in `cqed_sim`: add first-class
  `ConditionalDisplacement`, `JaynesCummingsExchange`, and
  `BlueSidebandExchange` primitives to `unitary_synthesis`, together with
  fast-eval support, API documentation, and unit tests.
- Remaining upstream target: add waveform-bridge / pulse-export coverage for
  native hybrid primitives so decomposition winners can be replayed through the
  pulse compiler without study-local glue.
- Add a helper for hybrid logical benchmark targets on qubit+cavity 2x2 blocks.
- Add a documented end-to-end recipe for synthesis-library comparison against
  GRAPE on a retained subspace, including native-gate families.

## Status
COMPLETE
