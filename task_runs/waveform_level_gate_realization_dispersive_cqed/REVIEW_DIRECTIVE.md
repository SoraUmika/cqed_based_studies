# Review Directive — Iteration 1
Date: 2026-04-02

## Decision
APPROVE

## Summary Verdict
This study gives a scientifically coherent and quantitatively supported waveform-level account of when dispersive cavity displacements and qubit rotations cease to behave like ideal gates. The physics claims are consistent with the reported parameter regime, the validation evidence is explicit, and the report now includes direct paper-level evidence for the echoed-displacement negative result that had previously lived only in side artifacts. The study is suitable for final archival use as a credible internal reference.

## Writing Quality Assessment
### Strengths
- The abstract, methods, and results sections now present the main quantitative findings without relying on external context.
- The regime-map narrative is clear, and the dominant error mechanisms are stated in physically meaningful terms.
- The echoed-displacement extension is now integrated into the paper text and appendix rather than being confined to execution logs and saved artifacts.

### Issues (must fix before approval)
- None.

### Suggestions (optional improvements)
- If the report is expanded further, consider tightening a few long appendix artifact descriptions to reduce underfull-line noise in the LaTeX log.
- A future revision could move the echoed-displacement comparison into a dedicated main-text figure if that extension becomes central rather than supplementary.

## Evidence-Claim Audit
| Claim (location) | Supporting evidence | Verdict |
|------------------|--------------------|---------|
| Excited-state displacement breaks down once |chi| T becomes appreciable (Results: Cavity displacement fidelity) | Fig. 1 and Fig. 2 | SUPPORTED |
| The dispersive shift chi dominates the displacement error budget (Results: Coupling ablation) | Fig. 3 and Appendix Table I | SUPPORTED |
| Qubit rotations become strongly Fock dependent away from the vacuum manifold (Results: Qubit rotation characterization) | Fig. 4 and Fig. 5 | SUPPORTED |
| DRAG improves vacuum-manifold fidelity but does not solve Fock-dependent detuning (Results: DRAG correction) | Fig. 6 | SUPPORTED |
| Populated-cavity qubit pulses generate substantial qubit-cavity entanglement (Results: Entanglement) | Fig. 7 | SUPPORTED |
| The naive echoed displacement fails for |+x> superpositions even though it can help the |e> branch alone (Results: Echoed displacement extension; Limitations) | Appendix Fig. 10 and Appendix Table IV | SUPPORTED |

## Physics and Methodology Assessment
### What is correct
- The reported device parameters are physically reasonable for a dispersive transmon-cavity platform.
- The dispersive and Fock-dependent error interpretations follow directly from the Hamiltonian and are consistent with the reported regime boundaries.
- Convergence is documented numerically rather than asserted qualitatively.
- The echoed-displacement negative result is framed correctly: the failure arises from the inserted qubit pulses becoming manifold dependent once the cavity is populated.

### What is problematic or missing
- No blocking physics issue remains.
- The closed-system scope remains a known limitation, but it is stated clearly and does not invalidate the reported unitary conclusions.

## Completeness Check
| Item | Present? | Notes |
|------|---------|-------|
| Reproducibility appendix | yes | Includes artifact inventory and reproduction procedure |
| All optimized or calibrated values needed for reproduction | yes | Pulse conventions and representative operating points are stated |
| Artifacts in artifacts/ directory | yes | 13 JSON artifacts plus NPZ data |
| scripts/reproducibility_notebook.ipynb | yes | Present and previously verified |
| IMPROVEMENTS.md updated | yes | Includes echoed-displacement limitation and follow-up actions |
| All questions from the Introduction answered | yes | Regime boundaries and dominant mechanisms are addressed |
| Limitations section specific | yes | Includes the echo failure mode and next-step guidance |

## Required Actions for Next Iteration
None. The study is approved.

## Open Concerns (non-blocking)
- Open-system and multi-waveform extensions remain valuable future work, but they are already framed as limitations rather than hidden assumptions.
- The MiKTeX environment on this machine lacks perl, so future automated report builds should continue to use explicit pdflatex and bibtex passes unless latexmk is repaired.
