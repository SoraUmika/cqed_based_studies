# Execution Summary
Date: 2026-04-06
Study: `studies/multitone_sqr_echo_rigorous_followup`
Run: `task_runs/multitone_sqr_echo_rigorous_followup`

## Scope
This run repeated the strict no-detuning simultaneous shared-line multitone-SQR study with a more rigorous echo focus:
- same strict physical model and no-detuning restriction,
- same controlled no-go and decoupled-block baseline,
- but now with a genuinely optimized echoed ansatz,
- explicit phase-sensitive metrics,
- a total-duration-matched direct comparator,
- and a representative manifold-aware refocusing benchmark.

## Main Quantitative Results
- Strict direct shared-line baseline:
  - mean restricted average gate fidelity: `0.7133`
  - best case: `0.8058`
  - mean probe-state fidelity: `0.7251`
- Exact reduced blockwise replay:
  - matched the full strict result to machine precision
  - minimum reduced-versus-full restricted process fidelity: `1.0`
- Stronger decoupled-block approximation:
  - fidelity `1.0` in every tested case
- Replayed ideal instantaneous echo:
  - mean restricted average gate fidelity: `0.2006`
  - mean max residual-`Z`: `0.0098 rad`
  - mean probe-state fidelity: `0.0886`
  - key lesson: residual-`Z` suppression alone is not evidence of gate success
- Jointly optimized ideal instantaneous echo:
  - mean restricted average gate fidelity: `0.5912`
  - best case: `0.8098`
  - beat the plain direct pulse in `8/16` cases, but only at `|chi|T/2pi = 5`
  - never realized the ideal gate
- Jointly optimized finite Gaussian echo:
  - mean restricted average gate fidelity: `0.3176`
  - mean max residual-`Z`: `0.4913 rad`
  - mean probe-state fidelity: `0.4177`
  - no case beat the active-duration direct pulse
  - only `3/16` cases beat the total-duration-matched direct comparator, and all at poor absolute fidelity
- Manifold-aware shared-line multitone refocusing subset:
  - refocusing pulse alone: mean fidelity `0.1805`
  - echoed construction using that pulse: mean fidelity `0.5073` on the four-case hard subset
  - better than the total-duration-matched direct pulse on that subset (`0.4627`) but still far below the plain direct baseline (`0.7080`)

## Interpretation
The stricter echo verdict is now more nuanced and more defensible:
- naive replayed ideal echo can make the residual-`Z` metric look excellent while the gate itself is terrible,
- genuine sequence-level ideal-echo optimization can help in some long-duration special cases,
- but that idealized improvement does not translate into a practical physical echoed gate once finite refocusing pulses are included,
- and the manifold-aware shared-line refocusing benchmark does not overturn that physical conclusion.

## Validation
- Sanity:
  - exact reduced replay equals the full strict shared-line result to machine precision
  - decoupled-block target match remains exact
- Convergence:
  - representative higher-budget direct and finite-echo reruns remained consistent with the production verdict
- Compute:
  - full production run wall clock about `1810 s`

## Main Conclusion
The earlier strict no-go conclusion survives, but the echo statement is now sharper:
- exact ideal simultaneous shared-line no-detuning SQR is still unsupported,
- ideal instantaneous echo is a useful upper bound and can modestly help some long-duration cases,
- physical finite echoed constructions still do not rescue the gate.
