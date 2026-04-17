# Task Checklist: Comprehensive Hybrid Control Survey

## Phase 1: Planning
- [x] Bootstrap study folder structure
- [x] Fix cqed_sim import issues
- [x] Write SCIENCE_DIRECTIVE.md (planning)

## Phase 2: Implementation — Core Benchmarks
- [x] Create common benchmark framework (model setup, metrics, plotting utilities)
- [x] Benchmark A: Ancilla-only control (X_pi) — Gaussian, DRAG at 20/40/80ns for n=0,2
- [x] Benchmark B: Cavity state preparation — SQR+SNAP vs SQR+D+SNAP vs D+SNAP synthesis comparison
- [x] Benchmark B2: GRAPE state preparation — Fock|1>,|2> at 500/700ns
- [x] Benchmark C: Conditional phase gate — Free dispersive evolution (F=1.0, T=176.1ns)
- [x] Benchmark D: GRAPE optimal control sweep — durations=[200,400,800]ns × targets=[1,2,3]
- [x] Benchmark E: GRAPE subspace unitary — n_match=2,3 (F=0.813, 0.976)
- [ ] Benchmark ECD-style control — BLOCKED: ConditionalDisplacement not available as gateset
- [x] Benchmark F: Chi sweep — 6 chi values × 3 methods (Gaussian, GRAPE, CPhase)

## Phase 3: Parameter Sweeps
- [x] Chi sweep (weak/intermediate/strong dispersive) — 18 results
- [x] Open vs closed system comparison — 8 results at 20/40/80/200ns
- [ ] Control budget sweep (amplitude, bandwidth, duration) — not implemented
- [x] Truncation convergence checks — 12 runs, results stable

## Phase 4: Validation
- [x] Sanity checks (analytical CPhase=1.0, DRAG>Gaussian, convergence stable)
- [x] Convergence analysis (Hilbert space dimension sweep)
- [x] Cross-method comparison (chi-dependent recommendations)
- [x] Full rerun reproducibility audit (fresh rerun matched archived numerical outputs exactly)

## Phase 5: Report
- [x] Write comprehensive survey + numerical results report (report.tex, 7 pages)
- [x] Generate all publication-quality figures (7 PNG+PDF pairs)
- [x] Create recommendation matrix (Table V in report)
- [x] Compile LaTeX report to PDF (455KB)
- [x] Create artifacts/benchmark_summary.json
- [x] Finalize IMPROVEMENTS.md and state files
- [x] Consistency review — ALL CHECKS PASSED (61 results verified)
