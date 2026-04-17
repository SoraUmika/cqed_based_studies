# Review Directive -- Iteration 2
Reviewer: Codex 5.4 xHigh Science Director
Study: studies/fast_robust_storage_vacuum_reset_comparison
Run: task_runs/fast_robust_storage_vacuum_reset_comparison
Date: 2026-04-14

## Decision
NEEDS_REWORK

## Journal Review Score
| Dimension | Score (1-5) | Blocking issue? |
|-----------|------------|----------------|
| E. Novelty & scientific significance | 3 | Yes |
| C. Technical soundness & methodology | 2 | Yes |
| A. Clarity & presentation | 2 | Yes |
| D. Reproducibility & completeness | 2 | Yes |

**Equivalent journal verdict:** Reject
**Scores:** 5 = publication-ready, 4 = minor fixes, 3 = significant improvement needed, 2 = major rework, 1 = fundamental flaw

## Summary Verdict
This revision successfully addressed most of the iteration-1 manuscript-structure problems: the report now has literature context, a real methods section, quantitative sensitivity tables, and a substantially improved appendix. The remaining blocker is more serious: the quantitative story is still internally inconsistent across the abstract, summary table, higher-Fock table, convergence table, and machine-readable artifacts. In the same manuscript, the pulsed ladder is reported as a single-photon cooler with $\tau_e = \SI{11.0}{ns}$ in the abstract and summary table, but as $\tau_e = \SI{995}{ns}$ for $|1\rangle_s$ in the higher-Fock table; similarly, the convergence prose claims at-most-factor-of-two changes for $\Delta t \leq \SI{1.0}{ns}$ while the printed convergence table shows changes of more than an order of magnitude for the bright-state and Raman-like schemes. Because the manuscript never explains whether these conflicting values are end-of-pulse, end-of-window, or steady-state metrics, the evidence-to-claim chain is still broken and the study is not yet scientifically defensible.

## A. Writing Quality Assessment
### Strengths
- The revised report now has a usable Introduction and System and Methods section, with appropriate cQED citations and a clearly stated robustness-score definition.
- The auxiliary autonomous benchmark remains clearly labeled as non-native and is not conflated with the physical protocol ranking.
- The main text generally keeps its claims scoped to the present device abstraction, which is the right scientific stance.

### Required Fixes (blocking approval)
- Results section: the same symbol $\tau_e$ is used for numerically incompatible quantities without explanation. The summary table gives the pulsed ladder $\tau_e = \SI{11.0}{ns}$, while Table~\ref{tab:initial_states} gives $\tau_e = \SI{995}{ns}$ for the pulsed ladder starting from $|1\rangle_s$.
- Results and Discussion: the narrative is internally contradictory. The summary table says the pulsed ladder is faster than the bright-state protocol (11 ns versus 13 ns), but the prose later states, "The bright-state protocol is the fastest single-photon cooler."
- Validation section: the prose claims time-step stability "by at most a factor of ${\sim}2$" for $\Delta t = 0.5$ and $\SI{1.0}{ns}$, but the table printed directly below it shows bright-state final storage occupation changing from $0.1925$ to $0.008558$ and Raman-like from $0.09898$ to $0.1091$.
- Reproducibility appendix: the report still does not tell the reader which artifact is canonical when the saved files disagree. For example, data/scheme_summary.csv reports the pulsed-ladder baseline residual as $0.004063$, whereas data/study_results.json reports selected_schemes.pulsed_ladder.final_storage_n as $0.008278$.
- Validation bookkeeping: TASK_CHECKLIST.md still leaves analytic-versus-numeric agreement and matched-noise-model stability unchecked, while README.md marks validation complete. The study state is therefore overstating completion.

### Suggestions (non-blocking)
- Replace the \texttt{cqed\_sim} code-style identifier in main-text prose with ordinary scientific wording.
- After the metric reconciliation, add one short sentence near the start of Results explaining which numbers are extracted at end-of-window, which are steady-state, and which are per-step or per-protocol quantities.

## B. Evidence-Claim Audit
| Claim (exact text, section) | Supporting evidence | Verdict | Action required |
|-----------------------------|---------------------|---------|----------------|
| "the pulsed ladder is the fastest practical protocol (e-fold time $\tau_e = \SI{11.0}{ns}$)" (Abstract) | Table~\ref{tab:summary}; data/scheme_summary.csv | UNSUPPORTED | Reconcile this with Table~\ref{tab:initial_states} and data/initial_state_summary.csv, which report $\tau_e = \SI{995}{ns}$ for the pulsed ladder starting from $|1\rangle_s$, or explicitly define the metrics as different quantities. |
| "the resonant bright-state path matches the pulsed ladder in speed ($\tau_e = \SI{13.0}{ns}$)" (Abstract) | Table~\ref{tab:summary}; data/scheme_summary.csv | WEAK | Clarify whether "matches" means same metric and same extraction window, and remove the later contradictory statement that the bright-state protocol is the fastest single-photon cooler. |
| "The bright-state protocol is the fastest single-photon cooler" (Scheme Comparison discussion) | None that is consistent with the summary table | UNSUPPORTED | Rewrite or remove. As printed, the report itself gives the pulsed ladder as faster (11 ns versus 13 ns). |
| "The pulsed ladder's e-fold time degrades sharply for $|3\rangle_s$ (from \SI{995}{ns} for $|1\rangle_s$ to \SI{1990}{ns})" (Higher-Fock and Mixed Initial States) | Table~\ref{tab:initial_states}; data/initial_state_summary.csv | SUPPORTED | Keep only after defining why this $\tau_e$ is not the same quantity as the 11 ns headline number. |
| "Between $\Delta t = 0.5$ and $\SI{1.0}{ns}$ ... final storage occupations change by at most a factor of ${\sim}2$" (Validation) | Table~\ref{tab:convergence}; data/convergence_summary.csv | UNSUPPORTED | Correct the prose or regenerate the table. The printed bright-state and Raman-like values contradict this statement by factors much larger than 2. |
| "the scheme ranking is stable for $\Delta t \leq \SI{1.0}{ns}$" (Abstract and Validation) | Table~\ref{tab:convergence}; data/convergence_summary.csv | WEAK | Demonstrate this only after the convergence table itself is reconciled with the headline summary values and metric definitions. |
| "the continuous Raman-like protocol achieves the highest robustness score ($0.982$) among the physical schemes" (Abstract) | Table~\ref{tab:summary}; data/scheme_summary.csv | SUPPORTED | none |
| "strengthening the readout dissipation from $1\times$ to $4\times$ reduces the Raman-like residual occupation by a factor of ${\sim}25$" (Abstract / Discussion) | Fig.~3; data/readout_kappa_tradeoff.csv | SUPPORTED | none |
| "The fastest practical scheme is the pulsed ladder ($\tau_e = \SI{11}{ns}$), which is ${\sim}20\times$ faster than the Raman-like protocol for single-photon cooling" (Conclusion) | Table~\ref{tab:summary}; data/scheme_summary.csv | WEAK | Keep only after reconciling the contradictory 995 ns single-photon value in Table~\ref{tab:initial_states}. |

Every UNSUPPORTED verdict is a required fix. Every WEAK verdict must be either strengthened or explicitly qualified in the revised report.

## C. Physics and Methodology Assessment
### What is correct
- The device parameters remain physically plausible for the stated cQED regime and are presented much more clearly than in iteration 1.
- The Raman-like protocol still appears to be the robustness leader among the physical schemes, and the readout-linewidth sweep still qualitatively supports the claim that a stronger readout bath helps it most.
- The revised report now documents the robustness-score construction and the operating-point scans sufficiently well for a reader to understand the intended workflow.

### Required Fixes (blocking)
- Define a canonical metric dictionary for the study. The current manuscript appears to mix distinct quantities such as end-of-window values, steady-state values, and possibly per-step or per-protocol times under the same labels.
- Reconcile the saved artifacts themselves. At minimum, data/scheme_summary.csv, data/initial_state_summary.csv, data/convergence_summary.csv, and data/study_results.json must agree on which quantity is being reported, or the report must explicitly separate them with distinct names.
- Rewrite the convergence section from validated numbers. Right now the printed table and the surrounding prose cannot both be correct.
- Resolve the validation-status mismatch between TASK_CHECKLIST.md and README.md. Either complete the unchecked validation items or reduce the scope of the claimed validation.

### Advisory Issues (non-blocking)
- The dispersive-model section still does not state a quantitative validity bound such as $g/\Delta \ll 1$ for the chosen operating regime.
- The manuscript still does not benchmark its headline performance against published active-reset results, so scientific significance remains suggestive rather than demonstrated.

### Convergence and Uncertainty Audit
- Hilbert space convergence: A table is now present, but it is not trustworthy in its current form because the baseline values do not align with the main summary and the prose misstates what the table shows.
- Optimizer convergence: Not central here. The operating-point scans are described, but the report should avoid stronger-than-supported near-optimality language.
- Uncertainty/error bars on key results: Sensitivity tables are present, but the underlying headline metrics are not yet consistently defined, so the uncertainty discussion is not currently interpretable.
- Multiple restarts / global optimum evidence (OPT/DES): No explicit global-optimum argument; acceptable only if the report stays with "selected operating point" language.
- Parameter sensitivity (+/-10%): Reported.
- Approximation validity bounds: Not stated quantitatively.

## D. Completeness Check
- Reproducibility appendix: Improved but still incomplete, because it does not map each reported headline number to a single canonical artifact field when multiple artifacts disagree.
- Saved artifacts in artifacts/: Present, but the quantitative cross-file inconsistencies currently prevent them from serving as a clean reproduction source.
- IMPROVEMENTS.md: Current with honest limitations.
- Notebook runs end-to-end: Not verified. The notebook file exists, but none of its cells are marked executed in the saved file.
- All figures referenced in text: Yes, based on the revised manuscript structure.
- All claims in abstract supported in body: No. The speed claims and the convergence claim are not self-consistent across the body.

## E. Novelty and Scientific Significance Assessment
- New insight delivered: The study still points toward a meaningful device-level split between the fastest experimentally simple protocol and the most robustness-favorable continuous protocol.
- Competitive with state-of-the-art: Not established quantitatively.
- Contribution delineated from prior work: Partly, but still weakened by the lack of a quantitative external benchmark.
- Scope accurately stated (system-specific vs. general): Mostly yes.
- Missing prior work that must be cited: No new citation omission is blocking in this iteration.

## Required Actions for Next Iteration
1. **[METRICS]** Define one canonical metric dictionary and propagate it everywhere.
   - What: For every reported quantity, explicitly distinguish end-of-window, steady-state, per-step, and per-protocol values. If the 11 ns and 995 ns pulsed-ladder numbers are both valid, they must be given different names and justified.
   - Where: report/report.tex plus a new or revised machine-readable crosswalk artifact in data/ or artifacts/.
   - Success criterion: Every number quoted in the abstract, summary table, higher-Fock table, and appendix maps unambiguously to one artifact field and one metric definition.
2. **[ARTIFACTS]** Reconcile or regenerate the inconsistent summary artifacts.
   - What: Bring data/scheme_summary.csv, data/initial_state_summary.csv, data/convergence_summary.csv, and data/study_results.json into quantitative agreement with the chosen metric definitions.
   - Where: studies/fast_robust_storage_vacuum_reset_comparison/data/ and studies/fast_robust_storage_vacuum_reset_comparison/artifacts/.
   - Success criterion: The single-photon baseline and convergence values no longer disagree across files, or any intentionally distinct quantities are clearly labeled as such.
3. **[REPORT]** Rewrite the contradictory Results, Validation, and Conclusion prose after reconciliation.
   - What: Remove statements that conflict with each other, including the "bright-state protocol is the fastest single-photon cooler" sentence if the canonical numbers still show the pulsed ladder faster.
   - Where: Abstract, Results, Validation, Discussion, and Conclusion.
   - Success criterion: No sentence in the report is contradicted by another sentence or table elsewhere in the same report.
4. **[VALIDATE]** Rebuild the convergence narrative from verified values.
   - What: Either regenerate the convergence table or correct the prose so it accurately describes the actual dt and truncation changes for all three physical schemes.
   - Where: Validation section and appendix.
   - Success criterion: The "factor of ${\sim}2$" and ranking-stability statements are numerically correct when checked directly against the printed table and underlying CSV.
5. **[STATE]** Bring the run-state files back into sync with the true validation status.
   - What: Update TASK_CHECKLIST.md, README.md, and the study state so unfinished validation items and notebook verification status are represented honestly.
   - Where: task_runs/fast_robust_storage_vacuum_reset_comparison/TASK_CHECKLIST.md, studies/fast_robust_storage_vacuum_reset_comparison/README.md, studies/fast_robust_storage_vacuum_reset_comparison/study_state.json.
   - Success criterion: The checklist, README, and study state all tell the same story about what has and has not yet been validated.

## Open Concerns (non-blocking)
- The report would benefit from a quantitative comparison against published active-reset figures of merit once the internal metric definitions are stable.
- A short quantitative dispersive-validity check would strengthen the methods section, but it is not the blocker in this iteration.
