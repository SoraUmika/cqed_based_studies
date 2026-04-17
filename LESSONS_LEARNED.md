# Lessons Learned

> Cross-study insights for all agents. **Read this file before starting any new study.**
> Updated as studies complete. Never delete entries — mark resolved items with `[RESOLVED]`.

---

## cqed_sim Quirks

| Issue | Detail | Workaround |
|-------|--------|------------|
| Model dimension mismatch in helpers | `run_conditioned_multitone_validation(..., n_tr=3)` fails comparing 2x2 targets against 3x3 reduced states | Use explicit qubit-subspace truncation or stick to `n_tr=2` |
| Phase factor sign conventions | Wait-phase inverse model requires `e^{-i phi_n(t)}` not `e^{+i phi_n(t)}`; sign errors silently break validation | Verify against exact closed-form analytic limits before deploying |
| Missing black-box inverse workflows | Fock-state inference, black-box gate forensics not first-class operations | Implement locally; upstream if reused across >1 study |
| Missing selective pulse builders | SNAP and ConditionalPhaseSQR pulse constructors not in main cqed_sim | Study-local definitions work; P1 upstream candidates |
| No reduced blockwise replay helper | Checking full vs. reduced blockwise dynamics requires custom code | Excellent validation sanity check; upstream candidate |

## Failed Approaches

| Approach | Why It Failed | Lesson |
|----------|---------------|--------|
| Global process-fidelity metric for conditioned phase | Undercounts usable gate quality because target only defined up to branch-local Z-phases | Use **branch-resolved relaxed fidelity** for structured controls |
| Direct Kerr-free correction via inferred chi' | Replaying with updated chi' did not improve recommendation set | Residual error is not dominated by chi' alone; multiple failure modes need separate diagnosis |
| Ideal instantaneous echo as rescue | Reduced Z-error but collapsed fidelity 0.71 -> 0.20 | Echoed refocusing is not universal; test on actual targets, not just residuals |
| Finite Gaussian echo for multitone SQR | Fidelity worse; maximum residual-Z 0.0786 -> 0.4146 rad | Echoes amplify certain block misalignments |
| Amplitude/azimuth optimization fallback | Multi-start Powell stuck at fidelity ~0.6 mean | When optimizer stalls uniformly, ansatz is insufficient; consider structural limits |
| Basis-expanded envelope family as catch-all | Never outperformed single-pulse Gaussian baseline | Simply adding degrees of freedom does not fix deep controllability limits |
| Symmetry-constrained echo variants | Phase-flip echo variant produced near-zero fidelity | Manual symmetry constraints often backfire; allow independent parameterization |

## Parameter Ranges That Work

| Parameter | Tested Range | Optimal/Typical | Notes |
|-----------|-------------|-----------------|-------|
| chi*T / (2pi) duration | 1-5 | 2-3.5 | Longer durations hit control bounds; no monotonic frontier |
| n_cavity (Fock truncation) | 4-12 | 8-12 | n_cav=4 for quick validation; 8+ for production |
| n_tr (transmon levels) | 2-3 | 2 (main), 3 (replay) | n_tr=2 for optimization; validate with n_tr=3 for leakage |
| Optimizer | Powell -> L-BFGS-B | Two-stage | Powell finds neighborhood; L-BFGS-B refines |
| maxiter | 50-200 | 200 | Can hit FACTR stop at 0 iterations; increase if stalled |
| Random ensemble size | 1-5 seeds | 5+ seeds | Single seed is fragile; 5 seeds reveal brittleness |

## Compute Timing Data

| Phase | Setup | Duration |
|-------|-------|----------|
| Single optimization (n_tr=2, Powell+L-BFGS-B) | Pilot SQR/CPSQR | ~7 s |
| Full 192-case SQR grid (single-pulse) | chi only, N_active 1-5 | ~49 min (15.4 s mean, 72.8 s max) |
| Echoed extension (180 cases) | Repeated-half refocusing | 1.5-26.9 s per case |
| Literature-informed pulse-family sweep | Noisy replication | ~369 s |
| Multitone no-go block-limit study | Production sweep | ~529 s |
| Black-box Fock inference | Optimized multitone + sweeps | ~105 s total |
| Primitive GRAPE diagnostics (N_cav 8-12) | Three seeds | ~7 s per case |

**Key notes:** No GPU acceleration in any study. Validation reruns typically 1-5% of main sweep cost. Optimizer convergence is the bottleneck for random targets.

## LaTeX / Build Notes

| Issue | Fix |
|-------|-----|
| `latexmk` hangs on Windows | Use direct: `pdflatex` -> `bibtex` -> `pdflatex` -> `pdflatex` |
| Stale duplicate `\end{document}` | Grep and keep only the first one |
| "Misplaced \noalign" in RevTeX appendix tables | Replace loop-generated rows with direct listings |
| UTF-8 BOM in Windows JSON | Use `encoding="utf-8-sig"` in all JSON loaders |

## Common Pitfalls

| Pitfall | Prevention |
|---------|-----------|
| Optimizer stagnation (0 iterations) | Increase maxiter; multi-start to rule out local minima |
| Process-fidelity definition mismatch | Use branch-resolved relaxed fidelity for structured targets |
| sign/phase convention errors | Cross-validate against hand-derived analytic limits |
| Echoed pulse sequences not universal | Test independently; don't assume echo uniformly improves |
| chi' higher-order corrections are subtle | Treat as separate experimental axis; separate chi-only and chi+chi' sweeps |
| Spectator excitation hypothesis | Measure crosstalk separately; don't assume it dominates |
| Identifiability nullspace in inverse problems | Add auxiliary measurement or impose prior |
| QuTiP import hangs on Windows | Use `scripts/runtime_compat.py` patch before any cqed_sim import |

## Upstream Candidates

- **Reduced blockwise replay helper** (P1, LOW) — Validates full vs. reduced-block equivalence
- **Blockwise residual-generator diagnostics** (P2, LOW) — Per-block X, Y, Z, and gauge analysis
- **Black-box Fock-state inference kernels** (P1, MEDIUM) — D(alpha) wait, recoverable-subspace analysis
- **SNAP and ConditionalPhaseSQR pulse builders** (P1, MEDIUM) — Study-local implementations work but awkward to reuse
