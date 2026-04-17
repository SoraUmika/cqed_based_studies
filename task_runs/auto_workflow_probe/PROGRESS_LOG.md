# Progress Log

## 2026-04-06T06:45:21Z - Study initialized via quickstart
- Objective: # Auto Workflow Probe Study

Study Name: auto_workflow_probe

Run a minimal cQED ANA study that starts from first principles and answers this question:

Can free dispersive evolution alone gener... (truncated, see prompt file)
- Study path: studies/auto_workflow_probe
- Run path:   task_runs/auto_workflow_probe

## 2026-04-06T01:51:11.6939693-05:00 - Recovery invocation resumed at plan phase
- Confirmed that the study was still in the bootstrap state with only B0.1 completed.
- Read research_config.json, AGENTS.md guidance, the probe README, study_state.json, checklist, progress log, blockers, and the recovery prompt before editing any files.
- Located the accessible local cqed_sim API reference under the current user's Box path after the stale configured path failed on this machine.

## 2026-04-06T01:51:11.6939693-05:00 - Plan phase completed for the lightweight probe
- Wrote task_runs/auto_workflow_probe/SCIENCE_DIRECTIVE.md with an analytic-first ANA plan, concrete hypotheses, and a minimal execution plan.
- Expanded studies/auto_workflow_probe/README.md to the required study format, including the first-principles dispersive derivation and a zero-gap cqed_sim plan.
- Normalized studies/auto_workflow_probe/IMPROVEMENTS.md to the required section structure and added lightweight future-agent notes.
- Updated studies/auto_workflow_probe/study_state.json to status PLANNED with explicit pending implementation tasks and compute budget.
- Added the missing studies/auto_workflow_probe/artifacts/ directory and initialized the report skeleton files needed by the loop.

## 2026-04-06T07:05:00Z - Implement phase started for the lightweight probe
- Read the implementation-phase context again, including the study state, science directive, checklist, progress log, blockers, and the applicable script and task-run instruction files.
- Re-checked the local cqed_sim API reference for `DispersiveTransmonCavityModel` and `dispersive_phase` before writing any simulation code.
- Chose the smallest implementation that still exercises the loop correctly: one self-contained study script, one CSV data file, one JSON artifact, and one dual-format figure.

## 2026-04-06T07:05:00Z - Implement phase edits prepared
- Added `studies/auto_workflow_probe/scripts/free_dispersive_pi_probe.py` to compute the analytic pi time, evaluate the cqed_sim helper unitary on a minimal qubit-cavity Hilbert space, and save the probe outputs.
- Updated `studies/auto_workflow_probe/IMPROVEMENTS.md` with the concrete implementation footprint and expected runtime envelope.
- Marked the two implementation tasks complete in `task_runs/auto_workflow_probe/TASK_CHECKLIST.md`.
- Reserved validation, report-writing, notebook, and review-request work for later phases to keep this invocation aligned with `phase=implement`.

## 2026-04-06T02:04:54.5164145-05:00 - Implement phase script executed successfully
- Ran `python scripts/free_dispersive_pi_probe.py` from the study root with no runtime errors.
- Generated `data/phase_difference_samples.csv`, `artifacts/free_dispersive_pi_probe_summary.json`, and `figures/phase_difference_vs_idle_time.{png,pdf}`.
- Measured analytic-versus-helper wrapped phase error of 8.882e-16 rad for `n_cav=2` and zero wrapped difference between `n_cav=2` and `n_cav=3` over the sampled scan.
- Recorded the representative pi-crossing time as 176.056338 ns and total wall time as 1.120 s.

## 2026-04-06T02:09:20.1367352-05:00 - Validate phase started for the lightweight probe
- Re-read the loop configuration, study README, study state, science directive, execution summary, checklist, progress log, blockers, and the saved probe outputs before changing any validation markers.
- Confirmed directly from the saved JSON artifact and CSV data that the probe still reports t_pi = 176.056338028169 ns, 51 sampled points, a maximum wrapped analytic-versus-helper mismatch of 8.882e-16 rad, and a zero wrapped difference between n_cav=2 and n_cav=3.
- Identified that study_state.json had drifted to the unsupported status `IMPLEMENTED`; planned a loop-compatible correction during the validation checkpoint so the next suggested phase becomes report rather than unknown.

## 2026-04-06T02:09:20.1367352-05:00 - Validate phase completed for the lightweight probe
- Marked V3.1-V3.3 complete in TASK_CHECKLIST.md using the saved probe outputs as evidence.
- Updated studies/auto_workflow_probe/README.md so the validation section now records the passed sanity and convergence checks and the not-applicable literature comparison.
- Updated studies/auto_workflow_probe/study_state.json to status REPORTING with explicit validation results and report/notebook work left pending.
- Refreshed task_runs/auto_workflow_probe/EXECUTION_SUMMARY.md and RESUME_PROMPT.md so the run state now points cleanly to the next report-phase invocation.

## 2026-04-06T02:28:23.0165884-05:00 - Report phase completed for the lightweight probe
- Read the existing report skeleton, created studies/auto_workflow_probe/report/report.tex.bak, and preserved the configured append-only report behavior by adding a minimal validated extension rather than overwriting the plan-phase scaffold.
- Appended the tiny probe report content, added the needed bibliography entries, and compiled studies/auto_workflow_probe/report/report.pdf successfully after simplifying the appendix artifact table for RevTeX compatibility.
- Created studies/auto_workflow_probe/scripts/reproducibility_notebook.ipynb and executed the fast load-first notebook cells to confirm the saved artifact, CSV data, validation checks, and figure regeneration all work from the notebook path.
- Updated TASK_CHECKLIST.md, EXECUTION_SUMMARY.md, and studies/auto_workflow_probe/study_state.json to reflect report completion and set the study status to REVIEW_REQUESTED.
- Wrote task_runs/auto_workflow_probe/REVIEW_REQUEST.md so the Science Director can review the tiny probe handoff.

## 2026-04-06T03:26:00-05:00 - Revision iteration 2 implement phase resumed
- Re-read research_config.json, AGENTS.md, the study README, study_state.json, FOLLOWUP_PROMPT.md, REVIEW_DIRECTIVE.md, TASK_CHECKLIST.md, PROGRESS_LOG.md, BLOCKERS.md, and the current implementation script before editing.
- Checked the cqed_sim surface again and confirmed that `DispersiveTransmonCavityModel.static_hamiltonian(frame=None)` is available, which provides a lightweight independent path for the required evidence addition.
- Chose a minimal revision scope for this phase only: extend the existing probe script with a static-Hamiltonian evolution comparison, regenerate the figure text without code-style labels, and update the study framing to an explicit model-level validation note.

## 2026-04-06T03:26:00-05:00 - Revision iteration 2 implement phase completed
- Extended `studies/auto_workflow_probe/scripts/free_dispersive_pi_probe.py` so it now produces both the original helper-backed outputs and a second figure-artifact pair based on explicit static-Hamiltonian evolution.
- Updated the primary figure labels to remove code-style text and added `data/phase_difference_hamiltonian_cross_check.csv`, `artifacts/free_dispersive_hamiltonian_cross_check.json`, and `figures/phase_difference_hamiltonian_cross_check.{png,pdf}` to the planned output set.
- Narrowed the README framing so the study is explicitly described as a model-level workflow-validation note rather than a hardware-regime claim, and logged the remaining report-phase gap in `IMPROVEMENTS.md`.
- Corrected the study-state iteration drift by moving the run to a validation-ready checkpoint after the revision implementation tasks completed.

## 2026-04-06T03:34:00-05:00 - Resolved the rotating-frame mismatch in the Hamiltonian cross-check
- The first static-Hamiltonian rerun used `frame=None` and therefore retained the bare cavity carrier, which produced an apparent wrapped disagreement of about `3.05 rad` against the intended dispersive phase law.
- Verified by direct basis-energy inspection that the correct model-level comparison must use the rotating frame defined by the bare qubit and cavity frequencies, where the branch shift reduces to `-2.84 MHz` as expected.
- Patched the probe script to evaluate the static Hamiltonian in that rotating frame, re-ran the script, and recovered machine-precision agreement: the maximum wrapped analytic-versus-Hamiltonian mismatch is now `8.882e-16 rad`, the helper-versus-Hamiltonian mismatch is `0.0 rad`, and the full script runtime is `1.076 s`.

## 2026-04-06T02:59:08.8211203-05:00 - Revision iteration 2 validate phase completed
- Re-read the saved helper and Hamiltonian artifacts, the cross-check CSV, and the current revision state before closing the validation gate.
- Ran a non-destructive Python validation check against the saved outputs and confirmed 51 sampled points, the exact zero-phase initial condition, and the closest sampled $\pi$ crossing at `176.056338028 ns` with the analytic, helper, and Hamiltonian phases all equal to `3.14159265359 rad`.
- Confirmed the saved quantitative validation results remain intact: the maximum wrapped analytic-versus-helper mismatch is `8.882e-16 rad`, the maximum wrapped analytic-versus-Hamiltonian mismatch is `8.882e-16 rad`, the helper-versus-Hamiltonian wrapped difference is `0.0 rad`, and the effective Hamiltonian manifold shift is `-2.8400000000000005 MHz`, matching the helper-derived shift within floating-point precision.
- Marked V3.4 complete in `TASK_CHECKLIST.md`, updated the README validation language to include the independent cross-check, refreshed `EXECUTION_SUMMARY.md`, and advanced `study_state.json` to the watcher-facing `REPORTING` status for the next report-phase invocation.

## 2026-04-06T03:15:25.2382992-05:00 - Revision iteration 2 report phase completed
- Refreshed `studies/auto_workflow_probe/report/report.tex.bak`, then rewrote `studies/auto_workflow_probe/report/report.tex` from the old append-only scaffold into a complete five-page model-validation note.
- Integrated the independent static-Hamiltonian evidence into the main text, moved the cutoff-stability check to the appendix, and removed the unsupported hardware-regime language so the report now matches the actual implemented scope.
- Rebuilt `studies/auto_workflow_probe/report/report.pdf` successfully after fixing two patch-serialization backslash losses and simplifying the appendix artifact table to a RevTeX-safe layout.
- Marked R4.5 complete in `TASK_CHECKLIST.md`, refreshed `EXECUTION_SUMMARY.md` and `REVIEW_REQUEST.md` for iteration 2, and advanced `studies/auto_workflow_probe/study_state.json` to `REVIEW_REQUESTED` with the prior `NEEDS_REWORK` decision cleared.
