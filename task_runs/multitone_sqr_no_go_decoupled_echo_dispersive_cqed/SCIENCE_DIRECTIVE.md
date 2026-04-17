# Science Directive - Iteration 2 (Extension Pass)
Date: 2026-04-06
Study: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed`
Run: `task_runs/multitone_sqr_no_go_decoupled_echo_dispersive_cqed`

## Study Objective
Use the archived completed pass as the accepted baseline and extend it only where the governing study objective remains under-resolved. The extension pass should strengthen two specific parts of the existing conclusion without redoing the validated baseline sweep: (1) make the nongeneric cancellation set behind the strict no-detuning no-go statement explicit, and (2) test whether the negative echoed-SQR verdict is robust once the comparison is restricted to symmetry-aligned targets and better-refocusing pulse conventions.

## Accepted Prior Baseline
The following results from the archived completed pass are accepted unless a targeted regression check contradicts them:
- The controlled block-resolved no-go derivation for strict simultaneous shared-line multitone SQR with amplitude and azimuth corrections only.
- The strict shared-line falsification sweep and its headline result that the best restricted average gate fidelity stayed well below unity.
- The exact reduced blockwise replay sanity check showing machine-precision agreement with the full strict shared-line evolution.
- The decoupled-block result showing unit-fidelity ideal SQR under the stronger spectator-dropped approximation.
- The baseline echo result: ideal instantaneous echo reduces some residual `Z` while degrading fidelity, and the tested finite vacuum-calibrated Gaussian echo performs worse than the plain strict pulse.
- The archived report, figures, data products, and reproducibility notebook as the validated starting point for this extension.

## Exact Extension Scope
1. Turn the existing "generic impossibility" claim into an explicit tuned-set map for the strict two-block problem so that the exceptional parameter relations are charted rather than only described abstractly.
2. Revisit the echoed construction only where the archived pass left a material limitation: whether symmetry-aligned targets and more manifold-aware finite refocusing choices materially change the practical verdict.
3. Preserve the archived 44-case strict sweep, validation package, and report as baseline evidence. Do not rerun the full baseline study unless a focused code change requires a regression recheck.

## Problem Classification
- Primary class: `ANA`
- Secondary class: `DES`
- Supporting class: `REP` for regression against the archived baseline results

## Physics Context
The physical question remains the same as in the archived pass. In the dispersive regime, a shared qubit-control line carries one tone per addressed Fock block at the corresponding manifold transition frequency. With no artificial per-tone detuning and no explicit `Z` compensation term, each block sees one resonant transverse contribution plus off-resonant spectator tones from the other addressed blocks. The extension work does not reopen that model choice. Instead, it tightens the interpretation of two already-identified consequences of that model: the structure of the accidental cancellation set and the exact limits of echoed refocusing.

## Analytic Preliminary
Start from the same block-resolved dispersive Hamiltonian used in the archived report,
```text
H_0 = sum_n (omega_n / 2) sigma_z tensor |n><n|
```
with a shared-line multitone drive at the exact manifold frequencies,
```text
H_d(t) = sum_m (Omega_m f(t) / 2)
    [exp(-i(omega_m t - phi_m)) sigma_+ + exp(+i(omega_m t - phi_m)) sigma_-].
```

For the strict two-block square-pulse problem, the archived pass already established
```text
zeta_0 = -lambda_1^2 K(Delta,T) + lambda_0 lambda_1 L(Delta,T,delta)
zeta_1 = +lambda_0^2 K(Delta,T) - lambda_0 lambda_1 L(Delta,T,delta).
```
The extension should use these equations more constructively:
- Solve or visualize the exact-cancellation relations as functions of duration, phase difference, and amplitude ratio.
- Distinguish the lower-dimensional tuned set from the open target-and-duration regions where the no-go is generic.
- Identify which special cases survive once the first-order transverse target constraints are enforced.

For the echo follow-up, use the toggling identities
```text
X_pi X X_pi = X,
X_pi Y X_pi = -Y,
X_pi Z X_pi = -Z,
```
so the half-SQR -> pi -> half-SQR -> pi sequence can preserve only aligned-`x` transverse structure exactly at first order. Generic `XY` targets are not symmetry matched to the echo. Even for aligned-`x` targets, noncommutation between the two half segments means that `Z` cancellation is only first-order unless stronger assumptions are added. A better finite refocusing design can test robustness of the archived practical verdict, but it does not change the exact toggling-frame algebra.

Controlled approximations for the extension:
- Analytic tuned-set mapping may continue to use the dispersive block-resolved model and Magnus/average-Hamiltonian reasoning.
- Exact shared-line follow-up checks must use the full compiled-waveform propagation and logical-subspace diagnostics already validated in the archived pass.
- Any manifold-aware finite echo design must state clearly whether it remains a shared-line physical construction or introduces a stronger effective-control assumption.

## Hypotheses
1. The exact-cancellation subset of the strict two-block problem is lower-dimensional and intersects nontrivial target families only on finely tuned relations rather than on open parameter regions.
2. Exact shared-line simulations near the tuned set will either confirm isolated accidental near-successes or show that full logical-subspace metrics still reject apparently good reduced diagnostics.
3. Symmetry-aligned or manifold-aware echoed constructions may reduce residual `Z` for aligned-`x` targets, but they will not restore exact generic `XY` SQR under the same strict shared-line no-detuning control model.
4. Any apparent echo improvement that survives must be judged on both restricted fidelity and blockwise residual-generator diagnostics; a single state-transfer score is not enough.

## Experiment Design
### Experiment 1: Tuned-Set Map For The Strict Two-Block Problem
- Purpose: Make the nongeneric cancellation set explicit and convert the qualitative no-go statement into a charted exceptional-set result.
- Method: Reuse the archived analytic formulas for `K` and `L`, construct amplitude-ratio / phase-difference / duration maps for the two-block case, and identify intersections with first-order target-matching constraints.
- Parameters: Start with `N_active = 2`, `chi_only`, then confirm representative slices under `chi_plus_chiprime`. Prioritize equal-angle and aligned-`x` slices because they are the most plausible accidental special cases.
- Expected outcome: Exact cancellation survives only on isolated or lower-dimensional tuned sets after the transverse target constraints are imposed.
- Success criterion: Produce machine-readable maps plus at least one figure that separates generic-failure regions from any analytically admissible tuned loci.

### Experiment 2: Exact Shared-Line Checks Near Tuned And Off-Tuned Points
- Purpose: Test whether the analytically identified tuned loci remain genuinely special under the exact compiled-waveform logical-subspace simulation.
- Method: Select a small checkpoint set from Experiment 1, run exact shared-line propagation with the existing targeted-subspace workflow, and compare full vs reduced replay, restricted fidelity, blockwise residual `Z`, and best-fit block-phase diagnostics.
- Parameters: Use a focused checkpoint set rather than rerunning the full archived sweep. Include both tuned and nearby off-manifold points.
- Expected outcome: The exact shared-line simulation either confirms isolated accidental behavior or shows that the apparent tuned points are not robust once the full operator metric is applied.
- Success criterion: Every claimed special case is documented with exact operator-level diagnostics and compared against a nearby off-manifold control point.

### Experiment 3: Echo Robustness For Symmetry-Aligned And Better-Refocused Cases
- Purpose: Resolve the strongest open limitation from the archived pass: whether the negative finite-echo result was mainly a vacuum-calibrated pulse artifact or a robust consequence of the strict shared-line model.
- Method: Start from the archived aligned-`x` cases, compare plain strict pulses, ideal instantaneous echo, the archived Gaussian finite echo, and at least one more symmetry-aware or manifold-aware finite refocusing construction if it can be expressed cleanly with the existing `cqed_sim` pulse/sequence runtime.
- Parameters: Restrict to a small aligned-`x` subset first. Only expand to broader target families if the aligned test reveals a credible improvement regime.
- Expected outcome: Better finite refocusing may improve residual `Z` for aligned targets, but any benefit should remain limited, approximate, and clearly narrower than an exact generic echoed-SQR rescue claim.
- Success criterion: Determine whether any finite echoed construction beats the corresponding plain strict pulse on both fidelity and residual `Z`; if not, strengthen the negative verdict with a sharper symmetry-based explanation.

## Execution Plan
1. **[ANALYZE]** Add focused extension analysis entry points without overwriting the archived baseline outputs.
   - Files to create: `scripts/run_extension_study.py`, optionally `scripts/validate_extension.py` if a dedicated validator keeps the baseline validator untouched.
   - Files to extend only if necessary: `scripts/common.py`.
   - Expected output: extension-specific JSON summaries, case artifacts, and figures saved under new extension-prefixed filenames.
2. **[RUN]** Compute the analytic tuned-set map for the two-block strict model.
   - Expected output: `data/extension_tuned_set_map.json` and at least one tuned-set figure in both `.png` and `.pdf`.
3. **[RUN]** Execute exact shared-line checkpoints on tuned and nearby off-tuned cases.
   - Expected output: extension checkpoint artifacts and a compact comparison table/summary.
4. **[RUN]** Execute the aligned-`x` echo robustness follow-up.
   - Expected output: extension echo comparison artifacts and at least one figure comparing plain vs echoed constructions on matched cases.
5. **[VALIDATE]** Run focused convergence and regression checks only on the new extension cases.
   - Required checks: timestep sensitivity for the strongest extension claim, full-vs-reduced consistency where applicable, and regression against archived baseline numbers for any reused case.
6. **[DOCUMENT]** Preserve the existing report and append an extension section only after the new results are validated.
   - Expected output in later phases: appended report section, updated notebook cells only if needed for the extension artifacts, and a new review handoff.

## Assumptions And Approximations
- Preserve the archived strict shared-line no-detuning model as the baseline physical question.
- Preserve the archived dispersive block-resolved analytic picture as the first-principles starting point.
- Treat the archived baseline report, figures, notebook, and machine-readable summaries as validated reference outputs rather than as work to redo.
- Any finite echo follow-up must state whether it remains a physically shared-line construction or whether it effectively adds extra control resources.
- Do not reinterpret the decoupled-block success as evidence against the strict no-go; it remains a stronger approximation used only as a contrast case.

## Known Risks
- The tuned-set map may reveal isolated special cases that require careful wording so they do not weaken the generic no-go statement by overextension.
- A manifold-aware finite echo may require additional calibration logic not already exposed cleanly in `cqed_sim`; if so, document the gap before adding local helpers.
- It is easy to over-read a small echo improvement if fidelity and blockwise `Z` are not judged together.
- Reusing archived code paths without isolated extension outputs could accidentally overwrite accepted baseline artifacts; avoid that.

## Stopping Criteria For This Iteration
- Stop after the tuned-set map is explicit, the echo robustness limitation is resolved as far as the existing runtime allows, and all new claims are backed by focused exact-simulation checks.
- Do not rerun the archived full baseline sweep unless a code change invalidates a baseline comparison.
- If the `cqed_sim` runtime lacks a clean public path for the finite refocusing experiment, document that gap, record the blocker, and continue with the tuned-set mapping and any echo analysis that remains valid.

## Compute Budget Estimate
- Analytic tuned-set mapping: low cost, expected minutes.
- Focused exact shared-line checkpoints: moderate cost, expected tens of minutes on CPU if kept to a small checkpoint set.
- Echo follow-up: moderate cost, similar to the checkpoint study if restricted to aligned-`x` cases.
- Total expected extension-pass runtime before reporting: comfortably below the archived full-sweep cost if scope remains focused.