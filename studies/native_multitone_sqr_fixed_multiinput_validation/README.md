# Native Multitone SQR Fixed Multi-Input Validation

## Problem Class
REP | ANA

## Motivation
The recent `cqed_sim` fixes changed the effective SQR conventions and, in at least one prior study, exposed a second configuration issue involving `fock_fqs_hz`. The older multitone SQR studies mostly reported blockwise or reduced metrics, plus a small number of state checks, but did not uniformly establish that the same native multitone pulse reaches the correct target state for multiple distinct qubit inputs on each addressed Fock manifold.

## Goals
1. Re-evaluate the representative native multitone pulse from each of the 4 most recent SQR studies under the fixed `cqed_sim` package.
2. Define a probe-state validation ladder that distinguishes success on a single favorable input from success on a spanning set of qubit inputs.
3. Measure whether each study's native multitone pulse implements the intended target state map for each addressed Fock level when started from multiple initial states.
4. Compare the fixed-package probe-state results to the original saved study-level metrics to identify which earlier conclusions remain convincing.
5. Provide a cross-study conclusion about whether native multitone SQR pulses now have robust evidence for input-consistent target-state preparation.

## Methods
- Use the current `cqed_sim` package to re-simulate representative native multitone pulses saved by:
  - `parameterized_waveform_residual_z_cancellation`
  - `multitone_sqr_arbitrary_fock_conditional_rotations`
  - `ideal_sqr_direct_vs_echoed_multitone`
  - `corrected_sqr_conditioned_rotation_metric`
- Attempt an analytic answer first:
  - For a block-diagonal intended qubit action on each Fock manifold, matching the image of a single input state `|g,n>` is not sufficient because infinitely many distinct SU(2) blocks share that same image up to an orthogonal action and phase.
  - Matching `{|g,n>, |+x,n>}` is stronger but still does not fully determine the action on the orthogonal direction.
  - Matching the spanning quartet `{|g,n>, |e,n>, |+x,n>, |+y,n>}` on every addressed manifold is a practical reduced-state analogue of qubit process tomography. If the quartet and the restricted process fidelity agree, that is strong evidence of the intended conditioned qubit action.
- Rebuild each run configuration with the fixed-package interpretation of `fock_fqs_hz`, i.e. do not pass the older frame-shifted override.
- For each study, simulate the saved native multitone pulse on the probe-state ladder and compare the full final state to the ideal target-state output generated from the study's stated target operator.

## Expected Outcomes
- Per-study tables of full-state fidelities for each probe input.
- A clear separation between:
  - favorable single-state success,
  - limited two-state success,
  - and strong multi-input evidence.
- A cross-study statement about whether native multitone SQR is convincingly demonstrated under the fixed simulator.

## Known Limitations
- This study re-validates representative native multitone pulse artifacts under the fixed simulator rather than re-optimizing every historical case from scratch.
- The four source studies use different target families, so cross-study comparisons are qualitative unless normalized by the common probe-state ladder.
- The validation intentionally ignores detailed cavity-subspace behavior beyond the addressed Fock labeling, except insofar as leakage lowers full-state fidelity.

Representative cases used in the completed run:
- `parameterized_waveform_residual_z_cancellation`: native baseline multitone, random target-D, `N_active = 4`, `|chi|T/2pi = 3`, `chi + chi'`
- `multitone_sqr_arbitrary_fock_conditional_rotations`: direct native multitone, structured family-C, `N_active = 4`, `|chi|T/2pi = 3`, `chi + chi'`
- `ideal_sqr_direct_vs_echoed_multitone`: direct native multitone, smooth x-axis target, `N_active = 3`, `|chi|T/2pi = 3`, `chi + chi'`
- `corrected_sqr_conditioned_rotation_metric`: corrected unitary-optimized direct multitone representative case, `N_active = 4`, `|chi|T/2pi = 3`

## Status
COMPLETE
