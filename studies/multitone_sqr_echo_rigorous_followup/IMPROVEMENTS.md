# Improvement Log: Rigorous Echo-Ansatz Follow-Up for Multitone SQR in Dispersive cQED

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- [P1][MEDIUM] **The ideal instantaneous echo is only an upper bound**: This follow-up showed that a jointly optimized ideal instantaneous echo can help in some long-duration special cases. Any future claim about a physical echoed gate must therefore keep the idealized and finite-pulse stories strictly separated.

## Recommended Improvements (P2)
- [P2][LOW] **Upstream the explicit logical probe metrics**: The replayed ideal echo looked successful on residual-`Z` alone but failed badly on the explicit probe states. Those diagnostics should be easy to reuse elsewhere.
- [P2][MEDIUM] **Map the long-duration ideal-echo sweet spots more carefully**: Optimized ideal instantaneous echo beat the plain direct pulse only at `|chi|T/2pi = 5`. The structure of that narrow regime deserves a dedicated study if idealized upper bounds remain interesting.

## Nice-to-Haves (P3)
- [P3][LOW] **Try richer but still physically plausible refocusing families**: The shared-line multitone `X_pi` benchmark did not work well, but a different finite refocusing family might clarify how much of the failure is intrinsic versus ansatz-limited.

## Open Questions
- Why does the ideal instantaneous echo help only at the longer duration and not at `|chi|T/2pi = 3`?
- Is there any finite refocusing family that can inherit the idealized long-duration improvement without collapsing on the physical metrics?
- Which of the surviving idealized improvements are robust and which are optimizer-tuned coincidences?

## What Was Tried and Did Not Work
- **Replayed ideal instantaneous echo**: The toggling-consistent replay crushed the mean maximum residual-`Z` error to `0.0098 rad`, but the mean restricted average gate fidelity was only `0.2006` and the mean explicit probe fidelity was only `0.0886`. It is a strong warning against using residual-`Z` alone as a success metric.
- **Finite Gaussian echo as a practical rescue**: Even after joint sequence-level optimization, the mean restricted average gate fidelity was only `0.3176`; no case beat the active-duration direct pulse.
- **Shared-line multitone refocusing pulse**: The separately optimized manifold-aware refocusing benchmark had mean fidelity only `0.1805` on its own, so it was already too poor to support an exact or near-exact echoed gate.

## Compute & Resource Notes
- Single favorable pilot case with light budgets: about `62 s`.
- Full production run with the final study budget: about `1810 s` on CPU.
- No extra packages were installed for this follow-up.

## Resolved
- **True optimized echo verdict obtained**: The follow-up replaced inherited-half replay with genuine sequence-level echo optimization and showed that the physical finite-pulse verdict remains negative.
- **Phase-sensitive metric gap closed**: The study now reports explicit probe-state fidelities alongside the framework gate and leakage metrics, which made the replayed-echo failure mode unambiguous.
- **Manifold-aware refocusing benchmark applied**: The strongest requested future-work item was tested on a representative hard subset and did not overturn the main conclusion.
