# Cluster-State Unitary Decomposition with Native Free Evolution and SNAP

## Problem Class
OPT | ANA

## Motivation
The cluster-state transfer unitary on the logical subspace
`{|g,0>, |g,1>, |e,0>, |e,1>}` is the key per-site operation for the
holographic cluster-state workflow already studied elsewhere in this repo.
Earlier work established that both selective phase gates and native
free-evolution gates can generate number-dependent phases, but it did not cleanly
separate which portion of the cluster-state unitary is carried by native
dispersive waiting and which portion must be supplied by explicit phase
corrections. This study isolates that question by first restricting the search
to `Displacement + QubitRotation + FreeEvolveCondPhase`, then extending the gate
set with `SNAP` to test whether SNAP materially reduces gate depth, reduces total
required entangling wait time, or mostly introduces redundant phase freedom.

## Goals
1. Quantify the best achievable fidelity for cluster-unitary decompositions using
   only `Displacement`, `QubitRotation`, and `FreeEvolveCondPhase` with
   wait-time optimization.
2. Compare that native-wait-only family against matched gate families that also
   include `SNAP`.
3. Attribute phase generation gate-by-gate so the report can identify which
   logical phases come from free evolution and which come from explicit phase
   gates.
4. Compare fidelity, gate depth, total sequence duration, total entangling wait
   time, explicit phase-gate count, and an implementation-complexity proxy
   across the gate-set choices.
5. Validate the main ranking against cavity truncation changes and block-phase
   diagnostics.

## Methods
- Reuse the cluster target, device parameters, and shared helpers from
  `studies/cluster_state_holographic_unified/scripts/common.py`.
- Use `cqed_sim.unitary_synthesis` primitives:
  `Displacement`, `QubitRotation`, `FreeEvolveCondPhase`, `SNAP`,
  `GateSequence`, `TargetUnitary`, `Subspace`, and `UnitarySynthesizer`.
- Use `GateSequence.phase_decomposition(...)`, `drift_phase_table(...)`, and
  `logical_block_phase_diagnostics(...)` to separate native drift phases from
  explicit programmed phase layers.
- Run structured ansatz sweeps over block count, initial wait-time priors, and
  multistart seeds for both the no-SNAP and SNAP-extended families.
- Run `scripts/run_improvement_pass.py` to re-optimize the winning native and
  interleaved-SNAP depth-6 families at `n_cav=6`, replay them with the
  `cqed_sim` pulse backend, and score a documented Lindblad noise surrogate
  built from the gate pulse unitaries and nominal `T1`, `Tphi`, and cavity-loss
  rates.
- Save machine-readable artifacts for the best sequences and the full comparison
  table.

## Expected Outcomes
- A defensible statement about whether `SNAP` is essential, helpful, or largely
  redundant once native free evolution is available.
- Best-sequence comparisons that report, at minimum: fidelity, gate depth, total
  duration, total free-evolution wait time, explicit phase-gate count, and a
  concise implementation-complexity assessment.
- Gate-by-gate attribution tables showing how the target unitary’s number-
  dependent phases are distributed between native waiting and explicit SNAP
  layers.

## Known Limitations
- The original `n_cav=4` winners have now been extended by a targeted
  `n_cav=6` re-optimization pass on the best native and best interleaved-SNAP
  depth-6 families. That pass improves the native case from replay fidelity
  `0.9177` to `0.9234` at `n_cav=6`, and the interleaved-SNAP case from
  `0.9473` to `0.9908` at `n_cav=6`. The family-wide search, however, has not
  yet been re-run at larger truncation, so threshold-optimal depth claims still
  rest on the original `n_cav=4` sweep.
- Coherent replay now exists, but it uses the `cqed_sim` pulse-unitary backend
  as a control surrogate rather than a compiled waveform bridge. On this build,
  the waveform bridge does not support `FreeEvolveCondPhase` or `SNAP`, so the
  current pulse-level check is a model-backed coherent replay rather than a
  hardware-distorted waveform simulation.
- Noisy validation is now available through a documented local Lindblad
  surrogate built from the gate pulse unitaries, nominal `T1=30 us`,
  `Tphi=20 us`, and cavity `T1=200 us`. This is informative for ranking, but it
  is not yet a direct `cqed_sim` compiled noisy replay of the full FE/SNAP
  sequence.
- The baseline study uses the dispersive model and logical `n=0,1` cavity block,
  so conclusions are specific to this encoding and simulator convention.
- Global optimality is not guaranteed; the search uses structured ansatz sweeps,
  one screening seed, and two refinement seeds rather than an exhaustive global
  search.

## Status
COMPLETE
