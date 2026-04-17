# Review Directive -- Iteration 1
Reviewer: Codex 5.4 xHigh Science Director
Study: studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed
Run: task_runs/multitone_sqr_no_go_decoupled_echo_dispersive_cqed
Date: 2026-04-06

## Decision
APPROVE

## Journal Review Score
| Dimension | Score (1-5) | Blocking issue? |
|-----------|------------|----------------|
| E. Novelty & scientific significance | 4 | No |
| C. Technical soundness & methodology | 4 | No |
| A. Clarity & presentation | 4 | No |
| D. Reproducibility & completeness | 4 | No |

**Equivalent journal verdict:** Minor Revision
**Scores:** 5 = publication-ready, 4 = minor fixes, 3 = significant improvement needed, 2 = major rework, 1 = fundamental flaw

## Summary Verdict
This report now makes the extension pass do exactly what it needed to do: it turns the previously abstract accidental-cancellation set into an explicit tuned-set map, and it corrects the earlier overstrong echo language to the more defensible "partial rescue only" verdict for symmetry-aligned checkpoints. The main technical claims are supported by the figures, appendix tables, and machine-readable summaries, and the extension materially strengthens the baseline conclusion rather than merely repeating it. The report also keeps the scope disciplined: it does not mistake the decoupled-block success for evidence against the shared-line no-go, and it does not oversell the manifold-aware finite echo improvement. The main remaining weakness is reproducibility polish rather than physics: the saved notebook is present but remains baseline-weighted and was not verified during this review as the extension-facing reproduction path.

## A. Writing Quality Assessment
### Strengths
- The abstract is self-contained and numerically specific; it states the strict control question, the analytic obstruction, the sweep-scale falsification attempt, the tuned-set extension result, and the refined echo verdict in one coherent narrative.
- The report consistently distinguishes the strict shared-line model, the stronger decoupled-block approximation, and the echoed construction instead of blurring them together.
- The extension section improves the paper's honesty: the tuned-set mapping is explicit, and the echo discussion now narrows the claim rather than stretching it.
- Figure captions are readable and mostly self-contained, with the extension figures correctly tied to the limited claims they support.

### Required Fixes (blocking approval)
- None.

### Suggestions (non-blocking)
- Add a one-sentence roadmap at the end of the Introduction so the section flow is explicit on first read.
- Consider adding a compact table of the effective model parameters used in the numerical study to make the Methods section easier to scan.
- During the polish pass, clean up the appendix path-heavy paragraphs that produced severe underfull-box warnings in the LaTeX log.

## B. Evidence-Claim Audit
| Claim (exact text, section) | Supporting evidence | Verdict | Action required |
|-----------------------------|---------------------|---------|----------------|
| "The strict simultaneous shared-line problem fails broadly: the mean restricted average gate fidelity is only 0.6094, the best case reaches 0.8058, and the worst falls to 0.3011." (Results, Headline results) | Table 1 and data/study_summary.json | SUPPORTED | none |
| "The exact reduced blockwise replay matches the full result to machine precision... The failure is therefore not caused by population leaving the addressed blocks." (Results, Headline results) | Table 1, Validation section, data/validation_summary.json, and reduced-vs-full metrics in saved case artifacts | SUPPORTED | none |
| "The strict shared-line fidelity decreases as the pulse grows longer... from 0.6500 at |chi|T/2pi=1 to 0.5689 at |chi|T/2pi=5." (Results, Headline results) | Fig. 1, Appendix Table I, and data/study_results.csv | SUPPORTED | none |
| "The residual phase grows with both duration and the size of the addressed manifold window." (Fig. 2 caption and adjacent text) | Fig. 2 and Appendix Table I | SUPPORTED | none |
| "The second-order Magnus approximation... is not numerically accurate enough to stand alone at the tested drive strengths." (Results, Why the numerical study supports, rather than replaces, the no-go) | Table 1 and saved Magnus-vs-full comparison metrics in data/study_results.json | SUPPORTED | none |
| "The decoupled-block construction reproduces the ideal target with unit fidelity in every case." (Decoupled-block section) | Fig. 3, Table 1, and data/study_summary.json | SUPPORTED | none |
| "The ideal instantaneous echo suppresses the mean maximum residual-Z error... but the mean fidelity collapses to 0.2018." (Echo section, Numerical comparison) | Table 1, Fig. 4, and data/study_summary.json | SUPPORTED | none |
| "The first nontrivial equal-angle aligned-x root occurs at |chi|T/(2pi)=0.7151483265621014." (Extension section) | Fig. 5 and data/extension_summary.json plus data/extension_tuned_set_map.json | SUPPORTED | none |
| "The exact shared-line checkpoint sitting on it is still poor: the restricted average gate fidelity is only 0.426986 and the maximum residual-Z error is 1.3824 rad." (Extension section) | Appendix Table III, Fig. 7, and data/extension_summary.json | SUPPORTED | none |
| "A manifold-aware finite multitone refocusing pulse improves both fidelity and residual Z relative to the plain strict pulse for the tuned and off-plus checkpoints... [but] remains far from ideal." (Extension section) | Fig. 6 and data/extension_validation_summary.json | SUPPORTED | none |

Every UNSUPPORTED verdict is a required fix. Every WEAK verdict must be either strengthened or explicitly qualified in the text.

## C. Physics and Methodology Assessment
### What is correct
- The report treats the no-go as a controlled result inside the block-resolved dispersive model rather than as an all-regime theorem, which is the right level of claim for the analytic machinery actually used.
- The two-block Magnus derivation isolates the obstruction in the correct place: spectator-induced blockwise Z terms that survive after the available amplitude and azimuth knobs are already consumed by the target transverse rotations.
- The exact reduced replay is used appropriately as a sanity discriminator. That comparison directly supports the statement that the strict failure is already present in the shared-line block-resolved dynamics rather than being a leakage artifact.
- The extension pass materially improves the physics narrative. Mapping the tuned set explicitly removes ambiguity about whether the earlier "generic" claim hid a missed exceptional locus, and the echo follow-up correctly retracts the older universal-failure wording without overshooting into a false rescue claim.
- The echo discussion handles symmetry carefully: generic XY targets are not symmetry matched to the toggling algebra, and the aligned-x improvement is framed as partial and special-case only.

### Required Fixes (blocking)
- None.

### Advisory Issues (non-blocking)
- The report could be clearer that the extension checkpoint study is a focused local probe around the tuned locus, not an exhaustive neighborhood search. The current wording is acceptable, but one additional sentence would prevent over-interpretation by a hurried reader.
- The numerical Methods section would benefit from a compact statement of optimizer restart policy so readers can immediately see that the study is a failed falsification attempt rather than a claim of global optimality.

### Convergence and Uncertainty Audit
- Hilbert space convergence: Reported with numbers for the representative baseline case; increasing the transmon truncation from two to three levels changes the cited fidelity only from 0.696774 to 0.696958.
- Optimizer convergence: Reported indirectly through the representative higher-budget rerun and through the report's non-optimality framing; this is adequate here because the headline claim is not "optimal gate found" but "analytic obstruction survives exact shared-line falsification attempts."
- Uncertainty/error bars on key results: No stochastic error bars are given, but the deterministic simulation uncertainty is bounded by explicit timestep and truncation checks, and those changes are small compared with the large gap from the strict-model fidelities to unity.
- Multiple restarts / global optimum evidence (OPT/DES): Not required for approval in the strict sense here because the report does not claim global optimality; the analytic no-go, not the optimizer alone, carries the central argument.
- Parameter sensitivity (+/-10%): Partially reported. The extension explicitly checks the tuned point with chi-prime included, but there is no broader +/-10% sensitivity sweep. This does not block approval because the report is scoped as a structural control-ansatz study rather than a hardware-tolerance study.
- Approximation validity bounds: Stated and respected at the model level. The report consistently limits its strongest claim to the dispersive block-resolved shared-line model and does not advertise it as a full strong-drive theorem.

## D. Completeness Check
- Reproducibility appendix: Complete for the report and script-based artifact trail; the extension files and reproduction steps are documented explicitly.
- Saved artifacts in artifacts/: Present and documented.
- IMPROVEMENTS.md: Current with honest limitations.
- Notebook runs end-to-end: Not verified in this review. The notebook exists and is structurally coherent, but the saved notebook has no executed cells and its content remains baseline-focused rather than directly reproducing the extension figures and extension summaries.
- All figures referenced in text: Yes.
- All claims in abstract supported in body: Yes.

## E. Novelty and Scientific Significance Assessment
- New insight delivered: The report now makes the exceptional two-block tuned set explicit and shows that this accidental cancellation locus still does not rescue the exact shared-line gate, while also sharpening the echo conclusion to "partial rescue only" for symmetry-aligned checkpoints.
- Competitive with state-of-the-art: Not applicable as a performance race. This is a constrained no-go and model-discrimination study, not a report of a new highest-fidelity control protocol.
- Contribution delineated from prior work: Yes. The report explicitly separates this strict no-detuning shared-line question from earlier repository studies that used detuning offsets or richer waveform families.
- Scope accurately stated (system-specific vs. general): Yes. The report repeatedly confines its strongest conclusion to the strict shared-line dispersive model and states what it does not prove.
- Missing prior work that must be cited: None that blocks approval at the chosen scope.

## Open Concerns (non-blocking)
- The reproducibility notebook should be refreshed in a later maintenance pass so the extension artifacts appear in the notebook itself rather than only in the report appendix and script-based reproduction path.
- The LaTeX build log still contains float-placement and underfull-box warnings. They do not undermine the science, but they are worth cleaning up during the polish pass.