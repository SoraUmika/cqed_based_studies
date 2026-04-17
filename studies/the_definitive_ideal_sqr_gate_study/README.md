# The Definitive Ideal SQR Gate Study

## Problem Class
OPT | REP | ANA

## Motivation
This study asks the repository's central selective-control question in one place: can a parameterized multitone qubit-drive waveform realize an ideal x-axis selective qubit rotation (SQR) in dispersive cQED, and if so on which addressed manifolds and durations? The study is intentionally self-contained and incorporates the three earlier SQR reports as explicit baselines rather than treating them as hidden prerequisites.

## Study Composition
| Component | Source study | Role here |
|---|---|---|
| Arbitrary SU(2) block targets | `multitone_sqr_arbitrary_fock_conditional_rotations` | Historical prior on broader block-diagonal control families |
| Residual-Z cancellation | `parameterized_waveform_residual_z_cancellation` | Historical prior on richer envelopes and echoed constructions |
| Ideal-SQR direct vs echoed baseline | `ideal_sqr_direct_vs_echoed_multitone` | The original strict ideal-SQR baseline |
| Native / rich SQR / CPSQR follow-up | `native_rich_multitone_sqr_cpsqr_feasibility` | Current strongest numerical evidence and "new results" layer for this definitive aggregation pass |

## Goals
1. Load and normalize all prior SQR study results into one machine-readable table.
2. Reconstruct per-manifold ideal-SQR error generators for the saved ideal-target artifacts.
3. Quantify whether the repository already contains strict ideal-SQR cases above 0.99 and where that success fails to generalize.
4. Add the missing lightweight diagnostics that were still absent from the first aggregation pass: standalone refocusing-pulse characterization, spectral-crowding analysis, and addressed-manifold scaling fits.
5. Generate a unified report, figures, validation summary, objective summaries, and reproducibility notebook in this study folder.

## Methods
- Load saved source-study CSV/JSON/artifact files with `scripts/load_prior_studies.py`.
- Recompute per-manifold Pauli error-generator components from saved ideal-SQR artifacts with `scripts/diagnostics.py`.
- Run standalone refocusing-pulse and spectral-crowding diagnostics with `scripts/extra_analyses.py`.
- Generate unified comparison figures with `scripts/plotting.py`.
- Aggregate source-study validation references and local sanity checks with `scripts/validate_results.py`.
- Build the report, ranking table, and reproducibility notebook with `scripts/run_full_study.py`.

## Analytic Preliminary
For an addressed manifold `n`, the desired strict target is `R_x(theta_n)`. A direct multitone waveform must therefore achieve the correct transverse rotation angle while also suppressing manifold-dependent coherent `Z` structure. This immediately suggests a split outcome:
- strict ideal SQR is possible only when both transverse and residual-phase errors are simultaneously small;
- echoed constructions can look much better on a relaxed CPSQR objective if they mostly refocus conditional phase while leaving some transverse mismatch unresolved.

That analytic picture matches the saved numerical corpus: the strongest direct solutions can exceed 0.99 on restricted addressed windows, while echoed gains are much more robust for relaxed CPSQR than for strict ideal SQR.

## cqed_sim Gap Analysis
| Functionality | Needed? | Available in cqed_sim? | Plan |
|---|---|---|---|
| Source-study simulations | Already completed | Yes | Reuse saved artifacts from the original studies |
| Cross-study normalization | Yes | No | Implement local study-level loaders |
| Artifact-level ideal-SQR error decomposition | Yes | Partial | Use local analysis on saved propagators |
| Standalone pulse replay across manifolds | Yes | Partial | Reuse patched native-rich helpers built on the public API |
| Full fresh optimization campaign | Not in this definitive pass | Yes, but expensive | Leave as future work tracked in `IMPROVEMENTS.md` |

## Assumptions
- The saved source-study artifacts are treated as the authoritative numerical records for this iteration.
- Strict ideal-SQR diagnostics use the blockwise addressed-manifold ordering already adopted by the source studies.
- Practical ranking in this pass uses a decoherence-only analytic fallback and does not yet include finite-difference parameter-drift sweeps.
- The standalone refocusing-pulse audit is a local compromise-pulse scan, not a full echoed-grid re-optimization.

## Compute & Resource Strategy
- No fresh long multitone optimization campaign is launched in this definitive pass.
- Parallel work is used for artifact-level error decomposition; the runner exposes `--n-workers` and `--sequential`.
- Local diagnostics (refocusing-pulse sweep, spectral crowding, scaling fits) run after aggregation and keep the full rebuild on this machine near 30-45 seconds.

## Expected Outcomes
- Recover the original low-fidelity baseline for the direct-vs-echoed ideal-SQR study.
- Show that the later native-rich study already contains strict ideal-SQR cases above 0.99.
- Demonstrate that strict-SQR success is conditional rather than universal across the full tested grid.
- Identify whether the remaining hard cases look more like finite-refocusing and crowding limits than like a simple missing-baseline issue.

## Known Limitations
- This is a definitive aggregation-and-diagnostics pass, not a fresh robust-`X_pi` / GRAPE / hybrid optimization campaign.
- Open-system robustness and parameter-drift sensitivity are not re-simulated in this pass.
- The report's practical ranking is therefore provisional.

## Validation
- [x] Sanity checks
- [x] Convergence
- [x] Literature comparison (implemented here as explicit repository-baseline comparison across the prior SQR studies)

## Status
COMPLETE
