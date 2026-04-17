# Progress Log

## 2026-04-06 - Initialization and prior-work audit

- Created the study directory `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed` and the matching run directory.
- Wrote the mandatory `README.md`, `IMPROVEMENTS.md`, `study_state.json`, `TASK_CHECKLIST.md`, and this progress log.
- Audited the closest prior studies and confirmed a major scope mismatch with the present prompt:
  - `ideal_sqr_direct_vs_echoed_multitone` optimized `d_omega`.
  - `multitone_sqr_arbitrary_fock_conditional_rotations` also optimized `d_omega` in both the direct and echoed workflows.
  - `parameterized_waveform_residual_z_cancellation` explored richer waveform families and also allowed `d_omega`.
  - `the_definitive_ideal_sqr_gate_study` summarized broader native-rich and echoed constructions rather than the strict simultaneous no-detuning ansatz.
- Audited the live `cqed_sim` source in the sibling editable checkout and confirmed that the strict shared-line targeted-subspace evaluator uses the full compiled waveform and exact propagator machinery, so it is suitable for the present study once `d_omega` is disabled.

## 2026-04-06 - Implementation and production sweep

- Implemented the strict no-detuning study code on top of `cqed_sim`.
- Added local helpers for:
  - exact reduced blockwise replay of the compiled shared-line waveform,
  - a stronger decoupled-block approximation with spectator tones removed by construction,
  - blockwise axis-angle and residual-generator diagnostics,
  - the ideal and finite echoed multitone sequence replay.
- Ran the full production sweep and saved machine-readable summaries, per-case artifacts, and waveform dumps.
- Generated the report figures:
  - `duration_fidelity_tradeoff`
  - `blockwise_residual_z_vs_duration`
  - `addressed_subspace_scaling`
  - `plain_vs_echo_comparison`

## 2026-04-06 - Main findings

- Derived the controlled two-block no-go result:
  - second-order blockwise `Z` coefficients are explicit,
  - canceling both requires `lambda_0 = lambda_1` together with an additional tuned phase-duration relation,
  - therefore exact cancellation is nongeneric once the transverse target has already fixed the available amplitude and azimuth knobs.
- Strict shared-line falsification attempt:
  - mean restricted average gate fidelity: `0.6094`
  - best case: `0.8058`
  - worst case: `0.3011`
- Exact reduced blockwise replay matched the full strict model to machine precision, so the failure is already present in the block-resolved shared-line dynamics rather than being caused by leakage.
- The stronger decoupled-block model reproduced the ideal target with unit fidelity in every tested case.
- Echo findings:
  - ideal instantaneous echo reduced the matched-set mean maximum residual-`Z` error from `0.0786 rad` to `0.0135 rad`,
  - but mean fidelity fell from `0.7133` for the plain strict pulse to `0.2018`,
  - finite `40 ns` Gaussian echoes performed worse still in every matched case.

## 2026-04-06 - Validation, notebook, and report

- Ran the validation workflow and documented the representative convergence spot checks.
- Built the required reproducibility notebook in `scripts/reproducibility_notebook.ipynb`.
- Wrote the LaTeX report and compiled `report/report.pdf`.
- Updated the README and improvement log with final findings and limitations.
- Prepared the execution summary and review request for science-director handoff.

## 2026-04-06T08:46:30Z - Archived stale extension signal files
- SCIENCE_DIRECTIVE.md
- EXECUTION_SUMMARY.md
- REVIEW_REQUEST.md

## 2026-04-06T19:37:14Z - Extension planning pass completed
- Re-read the study README, improvement log, study state, archived directive, archived execution summary, archived review request, and the current report to preserve continuity with the completed pass.
- Treated the archived report, figures, data, notebook, and validation package as accepted baseline work rather than work to redo.
- Wrote a fresh extension-scoped `SCIENCE_DIRECTIVE.md` for the current run.
- Narrowed the new implementation scope to two justified follow-ups tied to the existing objective: explicit tuned-set mapping for the strict two-block no-go and a sharper echo robustness study for symmetry-aligned and better-refocused cases.
- Marked the study state as `PLANNED` so the next staged phase is implementation.

## 2026-04-06T20:00:38.8996994Z - Extension implementation and focused validation completed

- Implemented the extension-only runner `scripts/run_extension_study.py` and the focused validator `scripts/validate_extension.py` so the archived baseline scripts, data products, figures, report, and notebook remain untouched.
- Generated extension-scoped outputs only:
  - `data/extension_tuned_set_map.json`
  - `data/extension_checkpoint_summary.json`
  - `data/extension_echo_summary.json`
  - `data/extension_results.json`
  - `data/extension_summary.json`
  - `data/extension_validation_summary.json`
  - `figures/extension_tuned_set_map.{png,pdf}`
  - `figures/extension_exact_checkpoint_comparison.{png,pdf}`
  - `figures/extension_echo_followup_comparison.{png,pdf}`
- The equal-angle aligned-`x` tuned-set map identified the first nontrivial accidental root at `|chi| T / (2 pi) = 0.7151483265621014`, but the exact shared-line tuned checkpoint still achieved only restricted average gate fidelity `0.4269864615` with maximum residual-`Z` error `1.3824023561 rad`.
- Nearby off-tuned checkpoints stayed non-ideal: the `off-minus` case dropped to fidelity `0.3466337350`, while the `off-plus` case rose only to `0.5175766363` with residual `Z` still `0.9133177244 rad`. Including `chi'` in the tuned checkpoint left the result essentially unchanged.
- Echo follow-up sharpened the earlier conclusion rather than reversing it: a finite manifold-aware multitone refocusing pulse beat the plain strict pulse on both fidelity and residual `Z` for the tuned and `off-plus` aligned-`x` cases, but the best tuned result still stalled at fidelity `0.4799250330` with residual `Z` `0.7172409433 rad`.
- The focused validator initially failed because it read archived scalar metrics from a raw case artifact instead of the archived `data/study_results.json` row. After patching that bookkeeping bug, the validation completed successfully.
- Focused validation results:
  - archived baseline regression matched exactly (`0.0` difference in the tracked fidelity and residual-`Z` metrics),
  - reduced-vs-full tuned-case fidelity checks remained `1.0` across `dt = 1, 2, 4 ns`,
  - finite echo can partially help on aligned-`x` checkpoints, so the defended conclusion is now `partial rescue only`, not `no finite echo can ever help`.
- Updated the run-state files so the next staged phase is `report`, where the preserved LaTeX report can be extended with the new evidence and the review handoff can be refreshed.

## 2026-04-06T20:13:26.5190388Z - Report extension and review handoff completed

- Extended the preserved `report/report.tex` rather than replacing it.
- Updated the abstract, main echo conclusion, final verdict, and appendix so the report now states the tuned-set map explicitly and narrows the echoed claim to partial, symmetry-aligned rescue only.
- Added extension-focused report coverage for:
  - the explicit tuned-set map,
  - exact checkpoint comparisons near the tuned locus,
  - the symmetry-aligned echo follow-up,
  - focused extension validation tables,
  - extension reproducibility artifacts and scripts.
- Ran a report preflight pass to check the new figures, label references, and main-text prose constraints relevant to this update.
- Attempted the standard in-place `report.pdf` rebuild, but Windows prevented overwrite because `report.pdf` was locked by another process.
- Verified the updated LaTeX source by compiling successfully to `report/report_build.pdf` instead.
- Wrote `REVIEW_REQUEST.md`, refreshed `EXECUTION_SUMMARY.md`, updated `TASK_CHECKLIST.md`, advanced `study_state.json` to `REVIEW_REQUESTED`, and refreshed `RESUME_PROMPT.md` so the next staged phase is science-director review.

## 2026-04-06T20:38:35.1509754Z - Polish pass completed

- Re-read the approved review directive, checklist, progress log, blocker file, and preserved report before making only presentation-level changes.
- Backed up the approved report source to `report/report.tex.prepolish`.
- Polished `report/report.tex` by adding the missing introduction roadmap sentence, a compact Methods parameter table, a clearer local-probe sentence in the extension section, and a tighter reproducibility appendix layout.
- Tried the canonical `report.pdf` build path first, but MiKTeX still could not write `report.pdf` because another process holds a Windows file lock.
- Verified the polished report through a successful alternate build to `report/report_build.pdf`.
- Re-attempted `report_build.pdf -> report.pdf` after compilation; the lock still persisted, so `report_build.pdf` remains the verified handoff artifact for this pass.
- Cleaned up the path-heavy appendix paragraphs enough to remove the earlier severe reproducibility-section underfull-box warnings; the remaining log output is limited to a small set of non-blocking underfull boxes in dense narrative lines and the long `cqed_sim` reference entry.
- Updated `TASK_CHECKLIST.md`, `study_state.json`, and `POLISH_COMPLETE.md` so the study and run now terminate in the `COMPLETE` state.
