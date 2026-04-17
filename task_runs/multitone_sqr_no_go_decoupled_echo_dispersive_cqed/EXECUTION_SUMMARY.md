# Execution Summary — Iteration 2

## Tasks Completed
- [x] `ext1`: mapped the strict two-block tuned cancellation set and selected exact checkpoints.
  - Output: `data/extension_tuned_set_map.json`
  - Key result: the first nontrivial equal-angle aligned-`x` tuned root occurs at `|chi| T / (2 pi) = 0.7151483265621014`.
- [x] `ext2`: ran focused exact shared-line checkpoint comparisons for tuned and nearby off-tuned cases.
  - Output: `data/extension_checkpoint_summary.json`, `data/extension_results.json`, `data/extension_summary.json`, and extension-prefixed case artifacts under `artifacts/cases/`
  - Key result: the tuned shared-line checkpoint still reached only restricted average gate fidelity `0.4269864615` with maximum residual-`Z` error `1.3824023561 rad`; the best nearby off-tuned case rose only to `0.5175766363`.
- [x] `ext3`: tested symmetry-aligned echoed refocusing, including a manifold-aware finite multitone refocusing pulse.
  - Output: `data/extension_echo_summary.json`, `figures/extension_echo_followup_comparison.png`, `figures/extension_echo_followup_comparison.pdf`
  - Key result: `echo_finite_manifold_aware_multitone` beat the plain pulse on both fidelity and residual `Z` for the tuned and `off-plus` aligned-`x` checkpoints, but still stalled far from ideal (`0.4799250330` tuned fidelity, `0.5374049217` off-plus fidelity).
- [x] `ext4`: completed focused extension validation and archived baseline regression.
  - Output: `data/extension_validation_summary.json`
  - Key result: archived baseline regression matched exactly, reduced-vs-full tuned-case fidelity remained `1.0` across `dt = 1, 2, 4 ns`, and the extension pass is ready to be written into the preserved report.
- [x] `ext5`: extended the preserved report, ran a report preflight pass, and refreshed the science-director handoff files.
  - Output: `report/report.tex`, `report/report_build.pdf`, `task_runs/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/REVIEW_REQUEST.md`
  - Key result: the paper now states the tuned-set map explicitly, softens the echo verdict to `partial rescue only` for symmetry-aligned checkpoints, and is ready for science-director review.

## Tasks Failed
- `ext4-initial`: first validator run.
  - Error: `KeyError: 'restricted_average_gate_fidelity'`
  - Attempted fixes: traced the archived scalar metrics to `data/study_results.json` rather than the raw per-case artifact.
  - Resolution: patched `scripts/validate_extension.py` to read the archived result row and reran the validator successfully.

## Key Results
- The accidental second-order tuned set is real in the analytic map but does not translate into a high-fidelity exact shared-line gate.
- Tuned checkpoint: restricted average gate fidelity `0.4269864615`, maximum residual-`Z` error `1.3824023561 rad`.
- `off-minus` checkpoint: fidelity `0.3466337350`, maximum residual-`Z` error `1.7996181140 rad`.
- `off-plus` checkpoint: fidelity `0.5175766363`, maximum residual-`Z` error `0.9133177244 rad`.
- Adding `chi'` to the tuned aligned-block spacing changed the tuned checkpoint negligibly (`0.4269918200` fidelity, `1.3822314089 rad` residual `Z`).
- Best tuned finite echo construction: `echo_finite_manifold_aware_multitone` with fidelity `0.4799250330` and maximum residual-`Z` error `0.7172409433 rad`.
- Best off-plus finite echo construction: `echo_finite_manifold_aware_multitone` with fidelity `0.5374049217` and maximum residual-`Z` error `0.6908428112 rad`.
- Archived baseline regression differences were exactly `0.0` for the tracked fidelity and residual-`Z` metrics.
- The preserved report now includes the extension figures, focused validation summary, and expanded reproducibility appendix for the extension artifacts.

## Result Digest
- The tuned equal-angle aligned-`x` root occurs at `|chi| T / (2 pi) = 0.7151483265621014`.
- The exact shared-line tuned checkpoint remains poor despite that root, so the accidental second-order cancellation does not rescue the full gate.
- The nearby `off-plus` point outperforms the tuned point in the full exact model, but still remains far from an ideal SQR gate.
- The extension tightens the earlier claim about finite echo: some finite echoes can help on both fidelity and residual `Z`, but the help is only partial.
- The defended conclusion is now `partial rescue only`, not `universal finite-echo failure`.
- The archived baseline products were preserved throughout; all new outputs are extension-scoped.
- The updated LaTeX source compiled successfully as `report_build.pdf`; the existing `report.pdf` could not be overwritten because it was locked by another process.

## Reviewer Pre-Check
| Required Action | Addressed? | Evidence |
|----------------|-----------|---------|
| No prior `REVIEW_DIRECTIVE.md` existed for this extension pass; execute the new `SCIENCE_DIRECTIVE.md` scope end-to-end without overwriting the accepted baseline. | Yes | `scripts/run_extension_study.py`, `scripts/validate_extension.py`, `data/extension_summary.json`, `data/extension_validation_summary.json`, `figures/extension_tuned_set_map.pdf`, `figures/extension_exact_checkpoint_comparison.pdf`, `figures/extension_echo_followup_comparison.pdf` |

## Anomalies / Concerns
- The analytically tuned checkpoint is not the numerically strongest exact case; the `off-plus` checkpoint performs better while remaining clearly non-ideal.
- A finite manifold-aware echo can improve both fidelity and residual `Z` at aligned-`x` checkpoints, so any future report text must avoid the stronger claim that `no finite echo can ever help`.
- The extension validation is focused rather than exhaustive; it anchors the archived baseline, the new tuned checkpoint, and the new echo comparison family.
- The updated PDF path for this handoff is `report/report_build.pdf`; `report/report.pdf` remains the locked baseline file until that lock is released.

## Updated File Manifest
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/scripts/run_extension_study.py`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/scripts/validate_extension.py`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/extension_tuned_set_map.json`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/extension_checkpoint_summary.json`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/extension_echo_summary.json`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/extension_results.json`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/extension_summary.json`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/extension_validation_summary.json`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/figures/extension_tuned_set_map.png`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/figures/extension_tuned_set_map.pdf`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/figures/extension_exact_checkpoint_comparison.png`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/figures/extension_exact_checkpoint_comparison.pdf`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/figures/extension_echo_followup_comparison.png`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/figures/extension_echo_followup_comparison.pdf`
- Modified: `task_runs/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/TASK_CHECKLIST.md`
- Modified: `task_runs/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/PROGRESS_LOG.md`
- Created: `task_runs/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/REVIEW_REQUEST.md`
- Modified: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/report/report.tex`
- Created: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/report/report_build.pdf`
- Modified: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/IMPROVEMENTS.md`
- Modified: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/study_state.json`

## Compute Notes
- `run_extension_study.py`: about `47.6 s` wall clock on CPU.
- No extra packages or GPU backends were required.
- The focused validator included one bookkeeping fix and then completed cleanly.
- The updated report source compiled cleanly under the alternate job name `report_build`; the only in-place rebuild failure was the external lock on `report.pdf`.