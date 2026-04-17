# Review Request

Date: 2026-03-31
Study: `studies/fock_resolved_black_box_sqr_inference`
Run: `task_runs/fock_resolved_black_box_sqr_inference`

The execution stage is complete and the study is ready for critical review.

## What to Review
1. Verify that the main conclusion is stated honestly:
   the protocol recovers weighted transverse sector information but does not recover `p_n` or sector-resolved `Z_n`.
2. Check the evidence-claim mapping for the coherence-witness claim:
   the combined protocol should expose coherent cases through a large residual, while wait-only should remain blind.
3. Check the pulse-level interpretation:
   near-diagonal black-box outputs should be reconstructed accurately, while coherent and leakage-prone outputs should be presented as failure modes / diagnostic triggers rather than as successful reconstructions.
4. Confirm that the report keeps filenames and code-style identifiers out of the main scientific prose and confines reproducibility-specific file references to the appendix / reproducibility section.

## Primary Artifacts
- `studies/fock_resolved_black_box_sqr_inference/report/report.pdf`
- `studies/fock_resolved_black_box_sqr_inference/report/report.tex`
- `studies/fock_resolved_black_box_sqr_inference/data/study_results.json`
- `studies/fock_resolved_black_box_sqr_inference/data/validation_summary.json`
- `task_runs/fock_resolved_black_box_sqr_inference/EXECUTION_SUMMARY.md`
