# Execution Summary
Date: 2026-04-06
Study: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed`
Run: `task_runs/multitone_sqr_no_go_decoupled_echo_dispersive_cqed`

## Scope
This run tested the strict simultaneous shared-line multitone SQR question in dispersive cQED under the following control restrictions:
- one resonant tone per addressed Fock block,
- no artificial per-tone detuning,
- amplitude and azimuth corrections only,
- ideal target defined as a direct sum of pure blockwise `XY` rotations with no residual block-dependent `Z` phase.

The run also analyzed two alternatives:
- a stronger decoupled-block approximation that drops spectator tones exactly,
- the echoed sequence `half-SQR -> pi -> half-SQR -> pi`.

## Analytical Outcome
- Derived a controlled two-block Magnus no-go statement for the square shared-line pulse.
- The second-order blockwise `Z` coefficients are explicit:
  - `zeta_0 = -lambda_1^2 K + lambda_0 lambda_1 L`
  - `zeta_1 = +lambda_0^2 K - lambda_0 lambda_1 L`
- Canceling both requires `lambda_0 = lambda_1` together with the additional tuned relation `L = K`.
- Once the transverse target has already fixed the available amplitude and azimuth knobs, those extra conditions are nongeneric.
- Generalized the conclusion to many blocks as a controlled generic-impossibility statement: the exact-cancellation set has no open interior in target-and-duration space.

## Numerical Outcome
- Strict shared-line production sweep:
  - `44` base cases
  - `24` echo cases
  - runtime about `529.4 s`
- Strict shared-line results:
  - mean restricted average gate fidelity: `0.6094`
  - best case: `0.8058`
  - worst case: `0.3011`
  - mean best-fit block-gauge improvement in process fidelity: about `0.0023`
- Exact reduced blockwise replay:
  - matched the full strict shared-line result to machine precision
  - confirms the failure is already present in the block-resolved shared-line dynamics rather than being caused by leakage
- Decoupled-block model:
  - reproduced the ideal target with fidelity `1.0` in every tested case
- Echo:
  - plain matched-set mean fidelity: `0.7133`
  - ideal instantaneous echo mean fidelity: `0.2018`
  - finite Gaussian echo mean fidelity: `0.3366`
  - ideal echo reduced matched-set mean max residual-`Z` from `0.0786 rad` to `0.0135 rad`
  - finite `40 ns` Gaussian echo increased matched-set mean max residual-`Z` to `0.4146 rad`
  - neither echoed construction outperformed the plain strict pulse in any matched case

## Validation
- Sanity:
  - exact reduced replay versus full strict model: fidelity `1.0`
  - decoupled-block target match: fidelity `1.0`
- Convergence:
  - representative strict case stayed at `0.696774` restricted average gate fidelity under a larger optimization budget
  - changing `dt` from `2 ns` to `1 ns` shifted that case only to `0.697201`
  - increasing the transmon truncation from `2` to `3` levels shifted it only to `0.696958`
- Local prior-work audit:
  - earlier nearby "positive" ideal-SQR studies relied on `d_omega` or richer waveform families and therefore do not directly contradict the strict no-detuning conclusion

## Main Conclusion
The study supports a skeptical but clear verdict:
- exact ideal simultaneous shared-line multitone SQR is not generically available under the strict no-detuning amplitude-plus-azimuth ansatz studied here,
- the stronger decoupled-block approximation does permit ideal SQR, but it is a different physical problem,
- the echoed multitone sequence suppresses some residual `Z` accumulation only approximately and does not rescue the strict gate.

## Output Package
- Report: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/report/report.pdf`
- Reproducibility notebook: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/scripts/reproducibility_notebook.ipynb`
- Headline summary: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/study_summary.json`
- Validation summary: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/validation_summary.json`
