# Progress Log

## 2026-04-02 - Unconditional displacement integration

- Executed the new unconditional-displacement driver covering naive, fast, two-tone, echoed, and hardware-aware optimal-control protocols.
- Generated new artifacts: `unconditional_single_tone_summary.json`, `unconditional_two_tone_summary.json`, `unconditional_echo_summary.json`, `unconditional_optimal_control_summary.json`, `unconditional_protocol_comparison.json`, `unconditional_wigner_cases.json`, and `unconditional_displacement_summary.json`.
- Generated new figure pairs for branch mismatch, superposition entanglement, filter tradeoff, `chi` scaling, protocol comparison, and Wigner comparison.
- Headline outcome: the best simple interpretable protocol is the `20 ns` two-tone compensated pulse (`delta_alpha = 8.15e-4`, `F_{+x} = 0.99768` on vacuum), while the best overall mean fidelity on the broad state-test set comes from the constrained `40 ns` optimal-control pulse (`0.9575`).
- Updated the study README, improvement log, run-state files, report, and reproducibility notebook to align with the new unconditional-displacement scope and artifact chain.

## 2026-04-01 - Verification remediation

- Ran a full study-validator-style verification pass after the echoed displacement extension.
- Confirmed that the study structure, artifact metadata, figure pairs, notebook execution, and report section requirements all pass.
- Found one consistency issue: `report.tex`, `EXECUTION_SUMMARY.md`, and `REVIEW_REQUEST.md` still described the pre-extension study state and omitted the echoed displacement artifacts/results.
- Updated the report Limitations/Future Work section, Saved Artifacts list, and Reproduction Procedure so the report now reflects the completed echoed-displacement extension.
- Updated `EXECUTION_SUMMARY.md` and `REVIEW_REQUEST.md` with the echoed-displacement negative result and corrected artifact/figure counts.
- Next state: ready for Science Director review; no active technical blocker remains.

## 2026-04-01 - Phase 10: Echoed Displacement Extension

- Implemented echoed displacement study: D(α/2)→Xπ→D(α/2)→Xπ echo scheme to refocus chi-induced displacement error.
- Analytic pilot computed toggling-frame error scaling: echo should reduce O(χ²T²) → O(χ⁴T⁴/16), predicting 5–2000x improvement depending on duration.
- Four-way numerical comparison: bare vs echo, for qubit in |g⟩, |e⟩, and |+x⟩, across α∈{0.5,1.0,2.0} × T∈{10,20,50,100,200}ns.
- Tested both Gaussian and DRAG π-pulses at T_pi∈{20,40}ns, with both fixed-displacement-duration and fixed-total-duration fairness conventions.
- **Key negative result:** Echo does NOT improve fidelity for superposition states. The vacuum-calibrated π-pulse introduces Fock-dependent errors when the cavity is populated during the displacement. For |+x⟩, the echo consistently underperforms bare displacement (improvement ratio 0.1–1.0x). For |e⟩ alone, the echo provides modest 1.2–3.3x improvement at long durations.
- Generated 4 figure pairs (PNG+PDF): fidelity comparison bars, improvement heatmap, infidelity vs duration, π-pulse variant comparison.
- Saved echoed_displacement_results.json and echo_analytic_estimate.json artifacts.
- Updated IMPROVEMENTS.md with detailed negative-result entry.
- Updated README.md Known Limitations with echo finding.
- This confirms the P1 warning from the ideal-SQR echoed study: echoed constructions fail when inserted π-pulses are not uniform across active Fock manifolds.

## 2026-04-01T00:00:00Z - Planning completed and run initialized

- Source task arrived as an inline prompt rather than a task file; the planner inferred the new study slug `waveform_level_gate_realization_dispersive_cqed` and matching run directory.
- Reviewed AGENTS.md, the task-run state rules, and `research_config.json`. The run will follow the full initialize -> plan -> implement -> validate -> report -> notebook -> review loop with a maximum of 14 iterations.
- Audited reference study helpers before writing the plan. The main reusable patterns come from the corrected SQR metric study, the ideal-SQR direct-vs-echoed multitone follow-up, the literature-informed selective-primitives study, and the unified holographic study.
- Confirmed from the local `cqed_sim` API and physics-conventions docs that the active runtime uses rad/s and seconds, transmon-first tensor ordering, bare-frequency rotating frames by default, and an internal pulse-carrier convention where the public drive frequency is converted to the raw carrier at the pulse boundary.
- Confirmed the repo already has reusable helpers for Gaussian/DRAG pulses, Fock-resolved SU(2) block extraction, conditioned Bloch vectors, Wigner diagnostics, and pulse-sequence simulation.
- Front-loaded the execution plan around a source-level convention audit because the study prompt explicitly requires implementation truth over documentation truth.
- Next checkpoint 1: initialize the study folder and write the mandatory README and IMPROVEMENTS scaffold.
- Next checkpoint 2: write the convention-audit artifact and freeze the ideal gate and metric baseline before any broad sweep.

## Session 2 - Phases 0-1: Study Bootstrap and Convention Audit

- Created study directory structure: scripts/, data/, figures/, artifacts/, report/
- Wrote README.md with all mandatory sections including analytic preliminary, gap analysis, assumptions
- Wrote IMPROVEMENTS.md scaffold with priority headers
- Created common.py (~515 lines) with shared helpers: model/frame builders, pulse constructors, simulation wrappers, operator extraction, SU(2) analysis, displacement analysis, serialization
- Created runtime_compat.py for QuTiP 5.x compatibility
- Conducted convention audit: verified tensor ordering (transmon-first, idx = q*n_cav + n), Hamiltonian terms (chi, chi', K correctly signed), carrier convention (carrier_for_transition_frequency returns negative frequency per exp(+iwt) convention), drive coupling (H_drive = eps(t)*adag + eps*(t)*a)
- Saved convention_audit.json to artifacts/
- Fixed displacement pulse phase bug: epsilon = i*alpha/T requires phase = arg(epsilon), not zero

## Session 3 - Phase 4: Displacement Study

- Created and executed displacement_study.py (~340 lines)
- Key findings:
  - |g> fidelity: >0.999 for T <= 100 ns across all alpha values tested
  - |e> fidelity: severely degraded; drops to 0.06 at T=100ns, alpha=2.0; 0.0022 at T=200ns
  - Ablation: chi accounts for >99% of the displacement error (fidelity 0.999990 -> 0.506 for |e>)
  - chi' and K contribute at the 10^-4 level
  - Entanglement entropy = 0 for all eigenstate initial conditions (physically correct: dispersive H commutes with qubit number operator)
  - Phase-space trajectory: |g> final fid=0.9998, |e> final fid=0.2143
- Generated 5 figure pairs (PNG+PDF): infidelity heatmap, entanglement heatmap, ablation, phase-space, fidelity vs alpha
- Saved displacement_sweep.npz, displacement_ablation.json, displacement_trajectory.json

## Session 3 - Phase 5: Qubit Rotation Study

- Created and executed qubit_rotation_study.py (~350 lines)
- Fixed figure generation bug: removed unused `theta_target = float("pi")` line
- Key findings:
  - Fock-resolved X_pi at T=40ns: fid(n=0)=0.998759, fid(n=5)=0.483370
  - X and Y rotations have identical fidelity profiles (rotational symmetry)
  - DRAG optimal at 0.5 ns: fidelity 0.999841 (vs 0.998759 without DRAG, ~8x improvement in infidelity)
  - DRAG overtuning: fidelity degrades for coefficients > 1.5 ns
  - Entanglement at alpha=2.0: entropy = 0.54 bits, fidelity = 0.36
- Generated 5 figure pairs: Fock-resolved angle, Fock-resolved infidelity, DRAG sweep, entanglement vs cavity, error budget
- Saved qubit_rotation_fock_resolved.json, qubit_rotation_drag_sweep.json, qubit_rotation_entanglement_vs_alpha.json

## Session 3 - Phase 6-7: Cross-Regime Synthesis and Validation

- Created and executed cross_regime_synthesis.py (~350 lines)
- Fixed sanity check bug: session variable reuse (alpha=2.0 instead of 1.0 for |g>/|e> comparison)
- Regime map results:
  - X_pi max Fock for 99% fidelity: n=1 at 20ns, n=0 at 40ns+
  - Displacement: |g> safe at all tested durations; |e> only safe for T < 10ns
- Validation results:
  - Sanity checks: ALL PASS (ideal limit, eigenstate symmetry, qubit-state independence at short T)
  - Convergence: ALL PASS (n_cav delta 3.2e-8, dt converged, n_tr delta 3e-6)
- Generated 3 figure pairs: combined regime map, qubit rotation regime map, convergence
- Saved cross_regime_summary.json, displacement_regime_map.json, qubit_rotation_regime_map.json, validation_sanity.json, validation_convergence.json

## Session 3-4 - Phase 8: Report and State Files

- Created references.bib (8 entries)
- Wrote report.tex (~430 lines) with all AGENTS.md-mandated sections
- Resolved compilation issues:
  - Removed `\usepackage{multline}` (already in amsmath)
  - Fixed "Extra \or" error: siunitx+revtex4-2 conflict with p{} columns; changed artifacts table to l l l with \footnotesize
  - Report compiles cleanly: 8 pages, no fatal errors, no overfull hbox
- Created reproducibility_notebook.ipynb with:
  - User-tunable parameter cell (all knobs in one place)
  - Derived objects cell (model/frame from parameters)
  - Load-saved-results + commented rerun paths for each major result
  - Validation section, convention audit, summary table
- Updated IMPROVEMENTS.md with all findings, bug reports, and compute notes
- Updated README.md: validation checkboxes checked, status = COMPLETE
- Updated TASK_CHECKLIST.md: all items through R8.5 checked
- Updated PROGRESS_LOG.md with full session history
- Next: Write EXECUTION_SUMMARY.md and REVIEW_REQUEST.md (R8.6)

## Session 5 - Verification pass and remediation

- Re-validated the completed study against the study-validator and validate-results requirements.
- Fixed workspace Python analysis so the external `cqed_sim` checkout resolves via `.vscode/settings.json -> python.analysis.extraPaths`.
- Upgraded `common.save_json()` so study JSON artifacts now carry the required metadata fields: `study_name`, `date_created`, `description`, and `load_instructions`, while `common.load_json()` remains backward compatible by returning the payload.
- Re-ran `displacement_study.py`, `qubit_rotation_study.py`, and `cross_regime_synthesis.py` to regenerate all machine-readable JSON artifacts with metadata.
- Backfilled `convention_audit.json`, which is not produced by a study script, to the same metadata schema.
- Repaired the reproducibility notebook to match the actual saved artifact schemas and executed all substantive code sections successfully: setup, parameters, derived objects, displacement, ablation, Fock-resolved rotation, DRAG sweep, entanglement, cross-regime summary, validation, and convention audit.
- Cleaned the report body so the main text no longer mentions internal framework identifiers, updated the Reproducibility appendix to list the actual artifact filenames, and rebuilt `report.pdf` successfully.
- Final verification status:
  - editor diagnostics for the study folder: no errors
  - JSON artifact metadata: present for all 11 JSON artifacts
  - notebook execution: successful through all substantive result sections
  - report compile: successful; only non-fatal underfull hbox warnings remain in the PDF build log
- Next: wait for Science Director review (`phase=review`). Phase 9 remains unstarted because no `REVIEW_DIRECTIVE.md` exists yet.

## 2026-04-02 - Final report integration and closeout

- Backed up `report.tex` and integrated the echoed-displacement negative result into the paper itself rather than leaving it only in the execution summary, notebook, and side artifacts.
- Added a dedicated Results subsection plus appendix evidence (`echo_improvement_heatmap` and a representative-value table) so the echoed-displacement limitation is now directly supported in the report.
- Cleaned the report preamble and minor LaTeX hygiene issues, enabled the RevTeX `floatfix` option, and rebuilt `report.pdf` successfully using explicit `pdflatex` and `bibtex` passes because `latexmk` is unavailable on this machine.
- Created the missing `study_state.json`, wrote `REVIEW_DIRECTIVE.md` with decision `APPROVE`, and wrote `POLISH_COMPLETE.md`.
- Final state: study marked complete, Phase 9 checklist items resolved, and no active blockers remain.

## 2026-04-02 - Multiplex-drive discussion extension

- Extended the completed unconditional-displacement report with an appendix-style discussion that interprets the successful two-tone pulse as the minimal multiplex cavity drive for the two-branch dispersive problem.
- Anchored the new discussion to existing study numbers rather than new simulations: the short two-tone result still solves the leading branch split almost exactly, while the broader-state gap to the constrained waveform now motivates richer multitone or multiplex follow-up work explicitly.
- Updated the study README, improvement log, checklist, and `study_state.json` so the same framing appears consistently in the archival study metadata.
- Next state: rebuild the PDF once and verify that the report still compiles cleanly after the documentation-only extension.

## 2026-04-02 - Explicit multiplex benchmark

- Added and executed `scripts/multiplex_displacement_followup.py`, which compresses the best `40 ns` constrained waveform into explicit full-duration multicarrier drives with `K = 2, 3, 4, 5, 6, 8` tones and evaluates them on the same broad 14-state test set.
- Generated `artifacts/unconditional_multiplex_followup.json` and the figure pair `figures/unconditional_multiplex_followup.{png,pdf}`.
- Main science result: the direct full-duration multiplex family is a negative result. Even the best `8`-tone case reaches only `0.524` mean fidelity and `1.89e-3` minimum fidelity, far below both the calibrated two-tone pulses and the bounded sampled optimal waveform.
- Secondary but important result: once the calibrated two-tone pulses at `20 ns` and `40 ns` are evaluated on the same broad state set, the `20 ns` two-tone case is actually the strongest tested protocol on that metric (`0.9857` mean fidelity, `0.9242` minimum fidelity).
- Next state: update the report and archival summaries so the old claim that the bounded optimal waveform is best overall on the broad state set is superseded by the new explicit two-tone and multiplex benchmark.

## 2026-04-02 - Structured multiplex closure

- Extended `scripts/multiplex_displacement_followup.py` to test two additional interpretable follow-up families: a segmented branch-resonant multiplex family with `1, 2, 4, 8` windows and a four-segment jointly optimized shaped-two-tone family with one complex scale per window.
- Confirmed from the live `cqed_sim` pulse implementation that carrier phase uses global time `exp(i (carrier * t + phase))`, so the poor segmented result is a physical negative result rather than a segment-phase convention bug.
- Final structured-follow-up results:
  - best full-duration multicarrier case (`8` tones): mean fidelity `0.5242`, minimum fidelity `1.89e-3`
  - best segmented branch-resonant case (`8` segments): mean fidelity `0.3612`, minimum fidelity `8.14e-3`
  - best jointly optimized shaped-two-tone case (`4` segments, Powell optimization): mean fidelity `0.1576`, minimum fidelity `0.0459`
- Reconciled the report main text, README, IMPROVEMENTS log, execution summary, review request, and checklist with the final ranking: the short `20 ns` two-tone pulse is the best tested protocol overall on the explicit broad state set, while the bounded `40 ns` optimal-control waveform remains the strongest sampled-waveform benchmark.
- Updated the reproducibility notebook to load the multiplex follow-up artifact, expose a rerun flag for the follow-up script, print the new structured-follow-up metrics, and display the combined multiplex figure. The edited notebook cells were executed successfully in the Python 3.12.10 kernel.
