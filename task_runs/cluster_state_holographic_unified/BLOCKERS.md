# Blockers

## Active
- None currently.

## Resolved
- Missing unified-study task-run state files were created so the rerun can be tracked and resumed cleanly.
- Corrected rerun execution initially failed because the generic blocks keyword was passed to builders that do not accept it; the dispatch logic was restricted to the CPSQR builder.
- Active-subspace propagation initially failed due to QuTiP tensor-dimension mismatches; propagation was normalized to flat Hilbert-space operators before cavity reductions were taken.
- The appendix Wigner generator was still reading legacy Strategy-B data; it now regenerates panels from corrected_scope_summary.json in the local figure directory.
