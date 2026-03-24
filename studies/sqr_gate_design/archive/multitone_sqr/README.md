# Simultaneous Multitone SQR Design

## Problem Class

OPT | ANA | DES

## Motivation

The earlier SQR study established that a single selective multitone waveform can
look reasonable under relaxed branchwise metrics while still failing as a strict
logical gate. The next natural question is more ambitious: can one common
multitone waveform implement multiple selective qubit rotations at the same
time, and if not, is the limiting defect a correctable phase structure or a
deeper waveform-control limitation?

This study focuses on simultaneous multitone SQR primitives acting on several
Fock branches in one shot. The main experiment-facing question is whether
simultaneous multi-branch SQR is viable as:

1. a direct small-angle logical primitive,
2. a compiled primitive with a cavity block-phase correction,
3. a pulse family that only becomes useful after waveform-side correction, or
4. a regime that requires richer segmented or optimal-control waveforms.

## Goals

1. Benchmark simultaneous multitone SQR on logical levels `n = 0..3` for
   multi-target sets `{0,1}`, `{0,2}`, `{0,1,2}`, and `{0,1,2,3}`.
2. Separate reduced conditioned-qubit performance from strict logical
   targeted-subspace performance.
3. Determine how simultaneous-gate performance depends on the requested
   rotation angle `theta in {pi/8, pi/4, pi/2, pi}` and on gate duration
   `chi T / 2 pi in {1, 2, 3, 5}`.
4. Test whether the dominant strict-logical defect is removable by an ideal
   cavity-only logical block-phase correction extracted from the simulated
   restricted operator.
5. Check whether simple waveform corrections can rescue hard large-angle cases:
   target-tone amplitude correction and multistart amplitude/phase/detuning
   optimization.
6. Validate that representative conclusions are stable to time-step and cavity
   truncation changes.
7. Replay representative optimized or baseline waveforms on a qutrit transmon
   model to estimate `|f>` leakage that is not captured by the two-level
   targeted-subspace workflow.
8. Produce a report that identifies when simultaneous multitone SQR is
   experiment-facing, when it is only useful as a small-angle primitive, and
   when it should be abandoned in favor of richer control families.

## Methods

### cqed_sim modules used

- `DispersiveTransmonCavityModel` for the dispersive transmon-cavity model
- `FrameSpec` and `manifold_transition_frequency` for frame-aware branch
  frequencies
- `ConditionedQubitTargets`, `ConditionedMultitoneRunConfig`,
  `run_conditioned_multitone_validation`, and `optimize_conditioned_multitone`
  for branch-local conditioned multitone analysis
- `build_block_rotation_target_operator`,
  `build_spanning_state_transfer_set`, and
  `run_targeted_subspace_multitone_validation` for strict logical-subspace
  validation
- `ConditionedMultitoneCorrections` and `ConditionedOptimizationConfig` for
  waveform-side pulse corrections
- `logical_block_phase_op` for ideal cavity-only compiled phase layers
- `SimulationConfig` and `prepare_simulation` for qutrit replay spot-checks

### Workflow

1. Run a baseline grid over target set, `theta`, and `chi T / 2 pi`.
2. Record reduced conditioned fidelity, strict logical fidelity, best-fit
   cavity block-phase compiled fidelity, state-transfer metrics, block
   preservation, and logical block-phase diagnostics.
3. Study one representative hard case with explicit target-tone amplitude scans
   and multistart waveform correction.
4. Replay representative compiled sequences on a qutrit transmon model to
   estimate transmon leakage.
5. Validate convergence in `dt` and `n_cav`.
6. Write and compile a manuscript-style report.

### API Check

The required simultaneous multitone workflow exists in `cqed_sim`:

- `cqed_sim.calibration.conditioned_multitone`
- `cqed_sim.calibration.targeted_subspace_multitone`

These provide the needed common-waveform builder, conditioned-qubit metrics,
strict targeted-subspace validation, and ideal cavity block-phase correction
analysis.

### Documented Gap

`targeted_subspace_multitone` currently evaluates conditioned qubit metrics via
two-level reduced qubit states, so the strict targeted-subspace workflow is not
directly compatible with `n_tr = 3`. The main logical-fidelity study therefore
uses `n_tr = 2`, and representative waveforms are replayed on `n_tr = 3` using
the general simulation stack to estimate `|f>` leakage.

## Assumptions & Convergence

- Storage cavity frequency: `omega_c = 2 pi x 5.241 GHz`
- Qubit frequency: `omega_q = 2 pi x 6.150 GHz`
- Anharmonicity: `alpha = 2 pi x (-255 MHz)`
- Dispersive shift: `chi = 2 pi x (-2.84 MHz)`
- Main logical window: `n = 0..3`
- Baseline targeted-subspace model: `n_tr = 2`, `n_cav = 6`
- Qutrit leakage replay model: `n_tr = 3`, `n_cav = 6`
- Time-step for production grid: `dt = 4 ns`
- Convergence checks: rerun representative cases at `dt = 2 ns` and
  `n_cav = 7`
- Objective emphasis for strict logical validation:
  `qubit = 0.3`, `subspace = 1.0`, `preservation = 0.25`, `leakage = 0.25`

## Expected Outcomes

- Small-angle simultaneous multitone SQR may be feasible for pairs and triples
  of branches, but large-angle one-shot simultaneous `X_pi` gates are likely to
  fail.
- Cavity-only block-phase compilation may provide only modest improvement unless
  the baseline waveform is already close to the intended logical operator.
- If large-angle failure persists even after simple pulse correction, then the
  natural next step is a richer segmented ansatz or direct optimal control.
- Representative qutrit leakage is expected to be small compared with the
  coherent logical mismatch.

## Status

COMPLETE

## Suggested Upstreaming

- Add an officially supported qutrit-compatible targeted-subspace multitone
  diagnostic path so simultaneous logical benchmarking does not need a separate
  replay pass for `|f>` leakage.
- Expose a public helper that converts a targeted-subspace validation result
  into a compiled logical-step operator suitable for repeated-step compilation
  studies.
