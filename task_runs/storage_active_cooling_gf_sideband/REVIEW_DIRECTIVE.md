# Review Directive — Iteration 1
Date: 2026-04-08

## Decision
APPROVE

## Journal Review Score
| Dimension | Score (1–5) | Blocking issue? |
|-----------|------------|----------------|
| E. Novelty & scientific significance | 4 | No |
| C. Technical soundness & methodology | 4 | No |
| A. Clarity & presentation | 4 | No |
| D. Reproducibility & completeness | 4 | No |

Scores: 5=publication-ready, 4=minor fixes, 3=significant improvement needed, 2=major rework, 1=fundamental flaw
Equivalent journal verdict: Accept

## Summary Verdict
The study gives a technically coherent effective-model answer to the target question: the ladder `|g,0_r,n_s> -> |f,0_r,n_s-1> -> |g,1_r,n_s-1> -> |g,0_r,n_s-1>` is viable through `n_s <= 4` in the local `cqed_sim` device abstraction, with explicit spectroscopy, pulse recommendations, and open-system cooling evidence. The strongest aspect is the evidence chain from resonance derivation to pulse scans to repeated-cycle cooling. The most important remaining limitation is also stated honestly: the work validates effective-control feasibility, not microscopic pump generation or hardware power calibration. One evidence-to-claim mismatch was found during review in the Step A pulse-family discussion and has already been corrected in the report, leaving no blocking issues within the stated scope.

## A. Writing Quality Assessment
### Strengths
- The report states the required `|q,n_r,n_s>` basis convention explicitly at the start of the methods section and keeps that notation consistent.
- The comparison to `arXiv:2503.10623v1` is specific enough to distinguish borrowed sideband-control logic from the modified cooling ladder studied here.
- The revised single-cycle cooling table now reports residual transmon and readout excitation in addition to success probability.

### Required Fixes (blocking)
- None.

### Suggestions (non-blocking)
- If this study is turned into a broader manuscript, add a short experiment-facing table that maps each recommended pulse to a suggested calibration sequence in the lab.
- If space permits, add one sentence clarifying why the long-duration phase-modulated Step A solutions are not preferred despite their competitive peak transfer.

## B. Evidence-Claim Audit
| Claim (exact text, section) | Supporting evidence | Verdict | Action |
|-----------------------------|---------------------|---------|--------|
| “The storage-cooling ladder ... is viable in the local device model for `n_s<=4`.” (Conclusion) | `data/recommendation_table.csv`, `data/study_results.json`, `figures/cooling_per_cycle.pdf` | SUPPORTED | none |
| “The best practical Step A controls are bump pulses for `n_s=1,3,4` and a short square pulse for `n_s=2`.” (Conclusion) | `data/recommendation_table.csv`, `data/study_results.json` | SUPPORTED | none |
| “The full open-system primitive achieves `96%` to `97%` single-cycle success ...” (Conclusion) | `data/recommendation_table.csv`, `data/study_results.json` | SUPPORTED | none |
| “The Floquet results therefore support the target-doublet picture and the absence of a large unexpected shift, but they should not be over-interpreted ...” (Floquet section) | `data/floquet_summary.csv`, `artifacts/floquet_summary.json` | SUPPORTED | none |
| “The phase-modulated bump ... performs similarly poorly in Step A.” (older Results wording) | `data/study_results.json` | UNSUPPORTED | fixed in current report: replaced with the accurate statement that Step A can reach high transfer, but only with much longer pulses |

## C. Physics and Methodology Assessment
### What is correct
- The Step A and Step B resonance conditions are derived from the same dispersive Hamiltonian used numerically.
- The bosonic enhancement claim is quantitatively verified through the extracted storage-sideband matrix elements and Floquet doublet splittings.
- The protocol distinguishes coherent transfer from irreversible cooling by explicitly including readout decay and tracking residual readout occupation.
- Convergence is reported numerically for truncation and timestep changes rather than being asserted without evidence.

### Required Fixes (blocking)
- None.

### Convergence and Uncertainty Audit
- Hilbert space convergence: reported quantitatively for transmon, storage, and readout truncation changes; all deltas remain small.
- Optimizer convergence: N/A in the strict optimization sense; the pulse-design study is a scan over discrete families, amplitudes, and durations rather than a stochastic optimizer.
- Uncertainty/error bars: absent, but acceptable within the present deterministic simulator scope because the study does not claim experimental uncertainty bounds.
- Multiple restarts / global optimum evidence (OPT/DES): N/A for the chosen scan-based workflow.
- Parameter sensitivity: reported through local detuning and amplitude scans.
- Approximation validity bounds: stated qualitatively and tied to the explicit limitation that the sidebands are effective control operators.

## D. Completeness Check
- Reproducibility appendix: Complete within the study scope.
- Saved artifacts in artifacts/: Present and documented.
- IMPROVEMENTS.md: Current and specific.
- Notebook runs end-to-end: Verified.
- All abstract claims supported in body: Yes.

## E. Novelty and Significance Assessment
- New insight delivered: the transmon-`f` storage-cooling ladder is shown to be viable in the local effective model, with explicit per-`n_s` spectroscopy and pulse recommendations rather than only a qualitative mechanism sketch.
- Competitive with state-of-the-art: not benchmarked as state-of-the-art control, but technically credible as an experiment-design study.
- Contribution delineated from prior work: yes; the report clearly distinguishes the present cooling ladder from the multimode-memory sideband protocol in `arXiv:2503.10623v1`.
- Missing prior work: none blocking within the report’s stated scope.

## Required Actions for Next Iteration
1. **[NONE]** No blocking follow-up is required for the current effective-model study.
   - What: The verification pass found one evidence-to-prose mismatch in the Step A pulse-family discussion; it has already been corrected.
   - Where: `studies/storage_active_cooling_gf_sideband/report/report.tex`
   - Success criterion: Current report wording matches `data/study_results.json` for the phase-modulated Step A family and the report remains compilable.

## Open Concerns (non-blocking)
- Hardware transfer still depends on a microscopic pump-calibration bridge that is explicitly out of scope here.
- Floquet leakage predictions remain qualitative because the dominant-doublet analysis carries a truncation-boundary warning.
- The passive `4/kappa_r` ringdown leaves a residual readout population near `1.7e-2`, which may matter for rapid experimental repetition.
