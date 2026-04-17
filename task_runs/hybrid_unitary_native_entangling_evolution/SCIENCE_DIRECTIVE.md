# Science Directive: Hybrid Unitary Native Entangling Evolution

## Decision
STOP

## Objective

The symbolic native-entangling study milestone is complete. Further work is now optional extension work rather than a required next phase for the established
`U_target = (I otimes H_c) CNOT_{c->q} CNOT_{q->c}`
on the logical subspace `{|g,0>, |g,1>, |e,0>, |e,1>}`.

## Planning Note

A dedicated high-reasoning planning model from research_config.json would produce better planning quality for this phase. The execution work below is still well-defined enough to proceed in the current session.

This directive is retained only as a handoff note for future extension work.

## Scientific Rationale

- Phase 2 established that the two-native-wait architecture is the correct native-heavy regime. One-wait candidates plateau at fidelity 0.25, while the best two-wait candidates remain near-unit fidelity.
- Phase 3.1 showed that direct waveform-bridge replay is blocked for the shortlisted native-heavy candidates.
- Phase 3.2 / 3.3 now completed the missing symbolic depth diagnostics. These confirm that:
   - `N2_exact_hc_to_exact_hc` is the symbolic upper bound and remains structurally faithful across the logical probe set.
   - `N2_A_local_to_A_local` is the best experimentally grounded symbolic candidate, with average final probe fidelity ~0.969 but visibly larger transient Bloch and Wigner distortions than the exact upper bound.
- The next scientific decision is no longer about missing diagnostics. It is whether the report should be finalized around symbolic gate-model validation with an explicit pulse-backed limitation, or whether one more engineering pass should be spent on building a pulse-backed validation path.

## Ordered Tasks

1. Promote a single report focal candidate.
   - Unless a better pulse-backed route appears immediately, use `N2_A_local_to_A_local` as the experimental focal point and `N2_exact_hc_to_exact_hc` as the symbolic upper bound comparator.
2. Decide the validation scope explicitly.
   - If a clean model-backed or pulse-backed replay path can be added quickly, do it.
   - Otherwise, finalize the report as a symbolic cqed_sim gate-model validation study with the pulse-backed gap called out in both `IMPROVEMENTS.md` and the report limitations section.
3. Draft the report section set around the now-established conclusions:
   - one native wait is insufficient,
   - two native waits are sufficient,
   - the archive `A_local` route remains the best experimentally grounded symbolic candidate,
   - depth diagnostics reveal larger transient distortion for the `A_local` route even when final fidelity stays high.
4. Record the exact pulse-backed blocker as a future-work item rather than leaving it implicit.

## Success Criteria for This Iteration

- The report can name a single experimental focal candidate and a single symbolic upper-bound comparator.
- The symbolic depth diagnostics are integrated into the study narrative rather than remaining only as standalone figure files.
- The remaining pulse-backed validation gap is described explicitly enough that a future agent can pick it up without re-deriving the current state.