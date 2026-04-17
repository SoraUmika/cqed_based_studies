# Progress Log: Comprehensive Hybrid Control Survey

## 2026-03-24 — Initialization
- Created study folder structure
- Fixed scipy 1.17.1 -> 1.14.1 for cqed_sim compatibility
- Read full cqed_sim API Reference (3846 lines)
- Created README.md, IMPROVEMENTS.md, study_state.json
- Next: Science Director planning phase

## 2026-03-24 — Science Director: Planning
- Classified problem as ANA+OPT+DES
- Designed 8-task benchmark suite (Tasks A-H)
- Wrote SCIENCE_DIRECTIVE.md with ordered action items
- Created benchmark_common.py with shared utilities

## 2026-03-24 — Execution Engineer: Implementation (Iteration 1)
- Created run_all_benchmarks.py with 8 benchmark tasks
- **Task A (Ancilla X_pi)**: 12 results — Gaussian/DRAG at 20/40/80ns for n=0,2. DRAG best: F=0.9998 at 80ns.
- **Task B (Synthesis)**: Restructured after TargetStateMapping failure. 9 results — SQR+SNAP F~0.975, D+SNAP F~0.970.
- **Task B2 (GRAPE state prep)**: Fock|1> 500ns F=0.9999, Fock|2> 700ns F=1.0.
- **Task C (Conditional Phase)**: F=1.0, T=176.1ns (analytic free evolution).
- **Task D (GRAPE sweep)**: 9 results — F>=0.9999 at 200ns, F=1.0 at >=400ns.
- **Task E (GRAPE unitary)**: n_match=2 F=0.813, n_match=3 F=0.976 (limited by 50 iterations).
- **Task F (Chi sweep)**: 18 results — GRAPE robust across all chi; Gaussian degrades at weak chi.
- **Task G (Convergence)**: 12 results — stable across all truncations.
- **Task H (Open system)**: 8 results — decoherence crossover at ~40ns.
- Total: 61 benchmark results in 439s wall-clock.

### Key debugging:
- TargetStateMapping always returns objective=2.0 — used make_target() workaround
- ConditionalDisplacement not available as gateset — ECD benchmarks omitted
- Subspace requires n_tr=2 (not 3) to match dimension constraints
- GrapeResult: use `objective_value` and `nominal_final_unitary` (not `objective`)
- Qobj to numpy: use `.full().flatten()` (not `np.array()`)

## 2026-03-24 — Execution Engineer: Figures & Report
- Generated 7 publication-quality figure pairs (PNG+PDF)
- Wrote report.tex (revtex4-2, 7 pages) with full appendix and reproducibility section
- Created references.bib (8 citations)
- Compiled PDF: pdflatex -> bibtex -> pdflatex -> pdflatex (455KB, 7 pages)
- Created artifacts/benchmark_summary.json

## 2026-03-25 — Execution Engineer: State File Updates
- Updated IMPROVEMENTS.md with P1/P2/P3 items, failed approaches, compute notes
- Updated study_state.json (status: VALIDATE, iteration 2)
- Updated TASK_CHECKLIST.md (marked all completed items)
- Updated BLOCKERS.md with resolved issues
- Next: Consistency review loop

## 2026-03-25 — Consistency Review & Completion
- Ran consistency_check.py: ALL 61 RESULTS PASSED
- Physical checks: DRAG >= Gaussian (PASS), monotonic fidelity (PASS), open <= closed (PASS), CPhase=1.0 (PASS)
- Convergence: fid_spread=0.000000 for all n_tr values (PASS)
- File checks: All 14 figures, report.pdf, convergence_results.json, artifacts present (PASS)
- Report table cross-check: All table values match benchmark data exactly
- Status set to COMPLETE

## 2026-03-24 — Rerun Verification
- Created archive backup at `archive/rerun_backup_20260324_222356` before rerunning the study pipeline
- Reran `python scripts/run_all_benchmarks.py`: completed successfully with 61 results in 331.9s
- Reran `python scripts/generate_figures.py`: regenerated all 7 figure pairs successfully
- Reran `python scripts/consistency_check.py`: ALL CHECKS PASSED on regenerated outputs
- Recompiled report via `pdflatex -> bibtex -> pdflatex -> pdflatex`: report.pdf rebuilt successfully at 455,471 bytes
- Compared archived baseline vs rerun outputs:
	- benchmark_results.json: zero numerical deltas across all 59 keyed result rows
	- convergence_results.json: zero numerical deltas across all 12 rows
	- PNG figures: byte-identical
	- PDF figures and report.pdf: identical sizes, non-identical hashes consistent with regenerated PDF metadata
- Fixed missing `sec:discussion` and `sec:conclusion` labels in report.tex and recompiled; `report.log` now has no undefined-reference warnings
