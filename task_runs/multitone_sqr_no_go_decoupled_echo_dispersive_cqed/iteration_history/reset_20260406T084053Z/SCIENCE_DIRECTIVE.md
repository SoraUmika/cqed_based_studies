# Science Directive - Iteration 1
Date: 2026-04-06
Study: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed`
Run: `task_runs/multitone_sqr_no_go_decoupled_echo_dispersive_cqed`

## Classification
- Primary class: `ANA`
- Secondary classes: `REP`, `DES`

## Core Question
Is an exact ideal simultaneous shared-line multitone SQR realizable in dispersive cQED when tones are fixed at the block transition frequencies, artificial per-tone detuning is forbidden, and the only free corrections are per-tone amplitude and azimuth?

## Prior-Work Audit Findings
1. Earlier nearby “ideal SQR” studies in this repository do not directly answer the strict question because they explicitly optimized `d_omega` or moved to richer composite families.
2. Some earlier echoed analyses were physically useful but narrower than the present prompt because they focused on `x`-axis targets or on waveform families outside the strict simultaneous shared-line model.
3. This study should cite prior work only after tracing each claim to its actual control knobs and metric.

## Hypotheses
1. In the strict simultaneous shared-line no-detuning model, off-resonant spectator tones induce blockwise effective `Z` generators that cannot be canceled generically with amplitude and azimuth knobs alone.
2. The stronger decoupled-block approximation removes that obstruction because each block retains only one resonant transverse drive.
3. The echoed sequence can suppress some `Z` accumulation only when the toggling assumptions are strong enough; it is not an exact universal rescue of ideal SQR.

## Ordered Tasks
1. **Formal no-go derivation**
   - Derive the blockwise interaction-frame Hamiltonian with all spectator tones retained.
   - Produce a two-block proof first, then generalize to multiple addressed blocks.
   - State clearly what “generic impossibility” means and list every assumption.
2. **Strict numerical falsification attempt**
   - Use `cqed_sim` full shared-line propagation.
   - Disable `d_omega`; optimize only amplitude and azimuth corrections.
   - Compare the exact operator against the ideal SQR using multiple metrics, not only a scalar fidelity.
3. **Decoupled-block counterexample**
   - Build a documented reduced helper that drops spectator tones block-by-block.
   - Show analytically and numerically that ideal blockwise SQR becomes achievable there.
   - Explain why this is no longer the same physical simultaneous shared-line problem.
4. **Echo analysis**
   - Define the exact `pi`-pulse convention.
   - Derive the ideal toggling-frame cancellation conditions.
   - Replay the finite-duration echoed sequence numerically and compare it against plain simultaneous multitone SQR.
5. **Validation and reporting**
   - Run sanity and convergence checks.
   - Write the report, compile the PDF, and create the reproducibility notebook.

## Numerical Design
- Representative model variants:
  - `chi_only`
  - `chi_plus_chiprime`
- Active windows:
  - `N_active = 2, 3, 4`
- Duration grid:
  - `|chi| T / 2pi = 1, 3, 5`
- Target families:
  - one structured `XY` family with nontrivial azimuth variation
  - one aligned-`x` family to give the echo sequence its best chance
  - a small random `XY` ensemble for falsification pressure

## Acceptance Criteria
1. The analytical section must identify the obstruction explicitly as a spectator-induced blockwise `Z` contribution, not merely as “optimization failed”.
2. The numerical section must report blockwise behavior and residual `Z` diagnostics, not only overall gate fidelity.
3. The decoupled-block section must prove and numerically confirm realizability.
4. The echo section must say whether the sequence works exactly, approximately, or only in special aligned regimes.
5. The final verdict must separate strict shared-line claims from stronger approximations and richer control families.
