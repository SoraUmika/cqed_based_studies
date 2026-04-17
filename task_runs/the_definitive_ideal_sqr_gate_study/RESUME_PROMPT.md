# Recovery Prompt - the_definitive_ideal_sqr_gate_study
Generated: 2026-03-28T06:43:46Z

---

## How to Use This File

1. Copy the **Agent Invocation** below and paste it into Copilot Chat.
2. The agent will read the state files listed under **Context Files** and
   continue from exactly where the loop left off.
3. Do NOT redo work that is already marked complete in TASK_CHECKLIST.md.

---

## Current Loop State

| Field            | Value                          |
|------------------|-------------------------------|
| Study            | the_definitive_ideal_sqr_gate_study          |
| Status           | INITIALIZED              |
| Phase            | needs-plan                        |
| Loop iteration   | 0 / 14           |
| Objective        | # Auto Research Prompt

# AI Research Prompt: The Definitive Ideal SQR Gate Study

## Role and Context

You are an autonomous cQED research agent with full access to the `cqed_sim` simulation framework and the full repository. Your task is to produce the **definitive, self-contained study** on whether parameterized waveforms can realize an ideal SQR gate in dispersive cQED. This study does not merely follow up on the baseline report â€” it **subsumes and extends** it. A reader of your final document should be able to understand the full arc of the investigation, from first principles through to the strongest achievable gate, without needing to consult any prior study. Every prior result should be either reproduced, extended, or contextualized within your work.

There are three prior studies in the repository that form the foundation:

1. `multitone_sqr_arbitrary_fock_conditional_rotations` â€” optimized arbitrary block-diagonal SU(2) targets (broader than ideal SQR).
2. `parameterized_waveform_residual_z_cancellation` â€” introduced echoed constructions and residual-Z objectives but used a different composite schedule and a scalarized loss rather than the pure ideal-SQR target.
3. `ideal_sqr_direct_vs_echoed_multitone` â€” the immediate baseline. Compared direct vs. echoed Gaussian multitone waveforms across a 16-case grid. Best result: F_avg = 0.72454 (direct). All echoed variants underperformed (~0.248).

**You must load and incorporate the results of all three prior studies.** Their data, figures, and conclusions should appear in your final report as contextualized historical baselines â€” not as separate documents the reader must hunt down. Where prior results remain valid, present them with proper attribution. Where they are superseded by your new results, show both and explain the improvement.

### The Ambition

The target is **F_avg > 0.99** on at least one physically credible case, and **F_avg > 0.95** as a mean across the primary grid. If this is not achievable within the tested ansatz families, you must demonstrate convincingly *why* â€” with a quantitative error floor analysis, not just an assertion â€” and report the tightest achievable bound. The baseline study settled for 0.72; that is no longer acceptable as a stopping point. Push every lever available.


## Parallelization Strategy (Mandatory)

The baseline study ran cases sequentially. This is unacceptable for the expanded scope. **You must parallelize aggressively at every level where the workload is embarrassingly parallel.**

### Level 1: Case-Level Parallelism

The study grid consists of independent (model_variant, target_family, N_active, duration) tuples. Each case is fully independent of every other case. Use Python `multiprocessing.Pool`, `concurrent.futures.ProcessPoolExecutor`, or `joblib.Parallel` to distribute cases across all available CPU cores. Detect the core count at runtime:

```python
import os
n_workers = max(1, os.cpu_count() - 1)  # leave one core for OS overhead
```

Every script that loops over the case grid must use this pattern. No sequential case loops in any production script.

### Level 2: Construction-Level Parallelism

Within each case, the four constructions (direct, symmetric echo, independent echo, asymmetric echo) plus any new constructions you introduce are independent of each other once the shared seed (nominal tone parameters) is computed. Parallelize construction optimization within each case. If the per-construction runtime is short (< 30s), batch multiple constructions into the same worker to avoid overhead. If long (> 2 min), give each construction its own worker.

### Level 3: Random-Seed Parallelism

For any optimization that may have local minima (especially the expanded ansatz and hybrid constructions), run multiple random restarts in parallel. Use at least 8 random seeds per construction per case for any new ansatz. Select the best result. This is trivially parallel.

### Level 4: Grid-Search Parallelism

For hyperparameter sweeps (envelope shape, sideband detuning ranges, refocusing pulse duration), distribute the grid points across workers. Use `itertools.product` to generate the full grid, then map it across the pool.

### Implementation Requirements

- Every parallelized script must include a `--n-workers` CLI argument that defaults to `os.cpu_count() - 1` but can be overridden.
- Every parallelized script must include a `--sequential` flag for debugging that forces single-threaded execution.
- Use `tqdm` or equivalent for progress tracking across parallel workers.
- All parallel workers must be stateless: no shared mutable state, no global variables modified in workers. Each worker receives its full input as arguments and returns its full output.
- Handle worker failures gracefully: if one case crashes, log the error and continue with the remaining cases. Do not let one failure kill the entire run.
- At the top of the main runner script, print the detected core count and the parallelization plan (how many cases, how many workers, estimated wall-clock time based on per-case estimates from the baseline).
- **Target total wall-clock time: â‰¤ 2 hours on an 8-core machine for the full study**, including all objectives. If any single objective would exceed this on its own, implement it but provide a `--quick` mode that runs a representative subset (e.g., 4 cases instead of 16) for rapid iteration, and a `--full` mode for the complete run.


## Study Structure: The Definitive Report

The final report must be a single, self-contained LaTeX document (revtex4-2, two-column) that a reader can pick up cold and understand completely. Structure it as follows:

### Part I: Foundation and Prior Work (Sections 1â€“3)

This is not a cursory literature review. It is a **complete re-presentation** of the relevant prior results, re-plotted in a unified visual style with consistent axis scales, color schemes, and notation.

1. **Introduction**: The ideal SQR gate, why it matters for cQED, and the precise mathematical target. State the central question clearly: "Can a parameterized multitone waveform realize U_SQR^ideal with F_avg > 0.99?"

2. **System Model**: The Hamiltonian, simulation parameters, and all assumptions. This section should be self-contained â€” do not refer the reader elsewhere for parameter values.

3. **Prior Results, Unified**: Load the results from all three prior studies and present them in a single comparative framework.
   - A unified table showing the best F_avg achieved by each prior study, what target it was optimizing, what ansatz it used, and what it did *not* test.
   - A single figure showing F_avg vs. construction type across all prior studies, with consistent formatting.
   - A clear "gap analysis" paragraph identifying every question that remains unanswered after the prior work.

### Part II: New Methods (Sections 4â€“6)

4. **Expanded Ansatz Family**: Every new waveform parameterization you introduce (sideband tones, alternative envelopes, extended durations, hybrid constructions). Each must have a clear mathematical definition.

5. **Upgraded Refocusing Pulses**: The manifold-robust X_Ï€ design, its optimization objective, and its standalone characterization.

6. **Optimization Infrastructure**: The parallel optimization framework, the multi-seed strategy, the staged optimization pipeline. Be explicit about the search budget: how many function evaluations, how many seeds, what optimizer, what convergence criterion.

### Part III: Results (Sections 7â€“10)

7. **Per-Manifold Error Anatomy**: The detailed error decomposition that the baseline study was missing. This is the diagnostic backbone of the entire report.

8. **Head-to-Head Construction Comparison**: The expanded grid results. Every new construction vs. every baseline construction, with matched-case analysis.

9. **The Push to F_avg > 0.99**: A dedicated section on the highest-fidelity results. What combination of ansatz, envelope, duration, and construction gets closest? What is the remaining error floor and what causes it?

10. **Robustness and Practical Viability**: Noise sensitivity, decoherence penalties, parameter drift tolerance. A construction that hits 0.99 but is unusable under realistic noise is not a solution.

### Part IV: Synthesis (Sections 11â€“12)

11. **Discussion**: What we learned, what the fundamental obstructions are, and whether a near-ideal SQR gate is achievable within *any* multitone ansatz or requires a qualitatively different approach. This section must leave the reader with no major unanswered questions about the tested ansatz families.

12. **Conclusion and Roadmap**: A crisp answer to the central question, plus a concrete prioritized list of what would need to change to close any remaining gap.

### Appendices

- Full per-case tables for all constructions (baseline + new).
- All optimized parameter sets for reproducibility.
- Waveform plots for every best-case construction.
- Convergence check results.
- Spectral crowding analysis for N_active = 2, 3, 4, 5.


## Primary Objectives

### Objective 1: Diagnostic Infrastructure (Execute First)

Build the per-manifold error decomposition and visualization pipeline before touching any optimization. This is the lens through which all subsequent results will be interpreted.

Tasks:

1. **Per-block error generator extraction.** For every (case, construction) pair in the baseline study, compute the 2Ã—2 error generator E_n defined by U_n^{realized} = exp(-i E_n / 2) Â· U_n^{ideal} for each active Fock block n. Decompose into Pauli components: E_n = Îµ_x^(n) X + Îµ_y^(n) Y + Îµ_z^(n) Z. Store as `{case_id, construction, n, eps_x, eps_y, eps_z, eps_norm}`.

2. **Apply the same extraction to all new constructions** as they are produced in later objectives. The extraction code must be a reusable function, not copy-pasted per objective.

3. **Generate the full diagnostic figure set** (details in the Visualization section below). These figures must be regenerated automatically whenever new results are added, not manually curated.

4. **Identify the dominant error channel.** For each case, rank the Pauli error components by magnitude. Across the grid, determine whether the dominant error is consistently Îµ_z (residual-Z), Îµ_x (rotation angle error), Îµ_y (axis tilt), or case-dependent. This classification drives the prioritization of subsequent objectives.

**Parallelize**: The error extraction for each (case, construction) pair is independent. Distribute across workers.


### Objective 2: Refocusing Pulse Upgrade

Tasks:

1. **Characterize baseline X_Ï€ manifold dependence.** For each model variant and each n = 0, 1, ..., N_active_max, simulate the standalone Gaussian X_Ï€ pulse and extract: realized rotation angle Î¸_n, rotation axis components (a_x, a_y, a_z)_n, and manifold-conditional process fidelity F_Ï€(n) relative to ideal -iX. Plot Î¸_n vs n and F_Ï€(n) vs n. Store diagnostics.

2. **Design a manifold-robust refocusing pulse.** Define objective:

   L_refocus = Î£_{n=0}^{N_active-1} || U_n^{X_Ï€} - (-iX) ||_F^2 + Î»_reg ||params - params_seed||^2

   Search over:
   - Pulse duration: T_Ï€ âˆˆ {20, 40, 60, 80, 100, 120} ns
   - Envelope: {Gaussian, flat-top Gaussian (with optimizable flat fraction), DRAG}
   - Per-manifold frequency correction: Î´Ï‰_Ï€^(n)
   - Amplitude and phase corrections

   Use 8 random seeds per configuration. Select the Pareto-optimal pulse (best worst-case F_Ï€(n) across manifolds).

3. **Re-run all echoed constructions** (symmetric, independent, asymmetric) with the robust X_Ï€ on the full grid. Compare matched-case against baseline echoed and baseline direct.

4. **If echoed + robust X_Ï€ approaches direct performance**, also re-optimize the echoed half-SQR segments jointly with the new refocusing pulse (full joint optimization, not staged). This tests whether the two components synergize.

5. **If echoed + robust X_Ï€ still underperforms**, decompose the remaining error per manifold (using Objective 1 infrastructure) and report exactly where the error lives now. Is it in the half-SQR segments? Commutator terms? Higher-order dispersive effects?

**Parallelize**: Refocusing pulse search (duration Ã— envelope Ã— seed) is a grid of ~144 independent optimizations. Distribute across all available cores. The subsequent echoed re-run is parallel at the case level.


### Objective 3: Expand the Direct Waveform Ansatz

The direct waveform's dominant failure is transverse/angle error, not residual-Z. Attack this directly.

Tasks:

1. **Sideband tones.** Extend from N_active to 2 Ã— N_active tones:

   Î©(t) = g_T(t) Î£_n [ Î»_n e^{i(Ï‰_n t + Î±_n)} + Î»_n' e^{i((Ï‰_n + Î´Ï‰_n') t + Î±_n')} ]

   The sideband detuning Î´Ï‰_n' is optimizable per manifold. Initialize Î´Ï‰_n' near Â±Ï‡/2 (half the dispersive shift), which is the natural scale for inter-manifold corrections. Run with 8 random seeds.

2. **Alternative envelopes.** Test:
   - Flat-top Gaussian: parameterized by (rise_time, flat_fraction), rise_time âˆˆ {T/8, T/6, T/4}, flat_fraction âˆˆ {0.3, 0.5, 0.7}
   - Blackman window
   - Tukey window: taper fraction âˆˆ {0.2, 0.4, 0.6}
   - Optimizable piecewise-linear envelope with 4â€“8 control points (this is the most flexible option)

   For each envelope, re-optimize tone parameters on the full grid. Report the best envelope per case and overall.

3. **Extended duration grid.** Push |Ï‡|T/2Ï€ to {3, 5, 8, 12, 20}. For each duration, run the best ansatz (from steps 1â€“2) with 8 random seeds. Plot F_avg vs duration. Determine:
   - Does fidelity saturate? At what level?
   - Does the error composition change with duration (e.g., residual-Z grows while transverse error shrinks)?
   - Extrapolate: what duration would be needed for F_avg > 0.99 under the current ansatz, if the trend continues?

4. **GRAPE or gradient-based refinement.** If `cqed_sim` supports gradient computation or if you can implement finite-difference gradients efficiently, apply a GRAPE-style optimization starting from the best Powell result. GRAPE with 50â€“100 time slices can access waveform shapes that no analytic ansatz can reach. If this is feasible, it provides the strongest possible upper bound on what the parameterization family can achieve. If `cqed_sim` does not support this, implement a simple finite-difference gradient wrapper and use L-BFGS-B.

5. **Combine everything.** Take the best envelope + sideband tones + longest practical duration + GRAPE refinement (if available) and report the single best F_avg achievable by a direct waveform.

**Parallelize**: Envelope Ã— duration Ã— seed is a large grid. Each point is independent. Distribute aggressively. The GRAPE refinement is per-case and can also be parallelized across cases.


### Objective 4: Hybrid Direct+Echo Construction

Tasks:

1. **Seeded hybrid.** Take the best direct waveform from Objective 3 and wrap it:
   Direct_segment â†’ X_Ï€(robust) â†’ Correction_segment â†’ X_Ï€(robust)
   where the correction segment is a short multitone pulse (duration â‰¤ 0.3 Ã— T_direct) with its own optimizable tones. Optimize only the correction segment initially.

2. **Joint optimization.** Release all parameters for joint optimization using a staged warm-start:
   Stage 1: Refocusing pulses alone (from Objective 2 results).
   Stage 2: Direct segment alone (from Objective 3 results).
   Stage 3: Correction segment alone.
   Stage 4: All parameters jointly, initialized from Stages 1â€“3.
   Use 8 random seeds at Stage 4.

3. **Double-echo and higher-order composites.** If the single echo layer helps, test:
   - Double echo: half â†’ X_Ï€ â†’ half â†’ X_Ï€ â†’ half â†’ X_Ï€ â†’ half â†’ X_Ï€
   - Nested: (half â†’ X_Ï€ â†’ half â†’ X_Ï€) â†’ Y_Ï€ â†’ (half â†’ X_Ï€ â†’ half â†’ X_Ï€) â†’ Y_Ï€
   These are more exotic but may access cancellation orders that single-layer echo cannot.

4. **Duration-fairness comparison.** The hybrid is longer than the direct waveform. For a fair comparison, also run the direct waveform at the same total duration as the hybrid. If the direct waveform at matched duration already beats the hybrid, the echo layer adds no value.

**Parallelize**: Each hybrid configuration Ã— seed is independent. Stage 4 joint optimization can be parallelized across seeds. The double/nested echo variants are independent of each other and of the single-echo hybrid.


### Objective 5: Robustness and Noise Sensitivity

Tasks:

1. **Quasi-static parameter sensitivity.** For every construction that achieves F_avg > 0.90, compute:
   - dF_avg/dÏ‡ via central finite differences (Î”Ï‡/Ï‡ = Â±0.5%, Â±1%, Â±2%)
   - dF_avg/dÏ‰_q (Î”Ï‰_q/2Ï€ = Â±50 kHz, Â±100 kHz, Â±200 kHz)
   - dF_avg/dÎ©_max (Î”Î©/Î© = Â±1%, Â±2%, Â±5%)
   Plot F_avg vs each parameter variation. Report the 1%-variation worst-case drop.

2. **Lindblad / decoherence simulation.** If `cqed_sim` supports Lindblad master equation or Monte Carlo trajectories:
   - Run the best 3 constructions with T1 = 50 Î¼s, T2* = 30 Î¼s (qubit), Îº^{-1} = 500 Î¼s (cavity).
   - Report the open-system F_avg.
   - If Lindblad is not available, use the analytic estimate: F_incoh â‰ˆ exp(-T_total/T1 - T_total/T2*) and flag it as approximate.

3. **Combined figure of merit.** Define:
   F_practical = F_unitary Ã— F_decoherence Ã— (1 - max_{param} |Î”F| under Â±1% variation)
   Rank all constructions by F_practical. This is the number that matters for hardware implementation.

4. **Sensitivity-fidelity Pareto frontier.** Plot F_unitary on x-axis vs worst-case sensitivity on y-axis. Identify constructions on the Pareto frontier (highest fidelity for a given robustness level, or most robust for a given fidelity level). These are the practical recommendations.

**Parallelize**: Each (construction, parameter_variation) evaluation is an independent simulation. The full sensitivity grid is embarrassingly parallel.


### Objective 6: Scaling to N_active = 4 and 5

Tasks:

1. **Spectral crowding analysis.** For each model variant, compute and plot:
   - Manifold transition frequencies Ï‰_n vs n for n = 0 to 7.
   - Frequency gaps Î”Ï‰_n = Ï‰_{n+1} - Ï‰_n vs n.
   - Tone bandwidth (1/Ïƒ for Gaussian, or the 3dB bandwidth for other envelopes) overlaid.
   - The "crowding ratio" CR_n = bandwidth / Î”Ï‰_n. When CR_n > 1, adjacent tones overlap significantly.

2. **Run the best 3 constructions from Objectives 2â€“4** for N_active = 4 (N_cav = 6) and N_active = 5 (N_cav = 7) on the extended duration grid. Use 8 random seeds.

3. **Determine the scaling law.** Plot best F_avg vs N_active for each construction. Fit to a functional form (e.g., F â‰ˆ 1 - a Ã— N_active^b or F â‰ˆ exp(-c Ã— N_active)). Extrapolate to N_active = 8 and report the predicted fidelity.

4. **If fidelity degrades sharply**, diagnose via the per-manifold error decomposition:
   - Is the degradation concentrated in the highest-n manifolds (spectral crowding)?
   - Does residual-Z grow faster than transverse error?
   - Does the optimizer get stuck? (Compare best-of-8-seeds vs median-of-8-seeds to assess landscape ruggedness.)

5. **Test manifold-selective strategies.** For N_active = 4+, try a "divide and conquer" approach: optimize tones for manifolds {0,1} and {2,3} separately, then combine with a reconciliation optimization. Compare against the monolithic optimization.

**Parallelize**: N_active Ã— duration Ã— construction Ã— seed is a large independent grid. Also parallelize the spectral crowding computation (trivial).


### Objective 7: Fundamental Limits Analysis

This objective exists to close the "is this an ansatz limitation or a fundamental one?" question definitively.

Tasks:

1. **Quantum speed limit estimate.** For each target angle set {Î¸_n}, compute the quantum speed limit (QSL) for the ideal SQR unitary. The QSL gives a lower bound on the gate duration for any control waveform, not just multitone Gaussian. If the durations you are testing are already near the QSL, further duration extension will not help.

2. **Information-theoretic parameter count.** Count the number of real parameters needed to specify U_SQR^ideal for N_active manifolds: this is 2 Ã— N_active (rotation angles plus the constraint that all axes are x). Count the number of free parameters in each ansatz. If the ansatz has fewer parameters than the target, the gate is generically unrealizable. If it has more, the excess parameters should provide optimization headroom. Report this ratio for every ansatz.

3. **Random unitary benchmarking.** Generate 20 random ideal-SQR targets (random {Î¸_n} drawn uniformly from [0.1Ï€, 0.9Ï€]) and optimize the best ansatz for each. Report the distribution of F_avg across random targets. This tests whether the baseline grid targets are unusually easy or hard.

4. **Upper bound via relaxed optimization.** If GRAPE / gradient methods are available, run an unconstrained piecewise-constant pulse optimization (no Gaussian envelope, no multitone structure â€” just raw time-domain control) for the best 3 cases. This gives an upper bound on what *any* single-channel qubit drive can achieve. If the upper bound is < 0.99, the obstruction is physical (e.g., the dispersive model itself cannot support an ideal SQR with qubit-only drive at this duration). If the upper bound is > 0.99 but the structured ansatz cannot reach it, the obstruction is in the ansatz.

**Parallelize**: Random target generation Ã— seed is embarrassingly parallel. GRAPE runs per case are independent.


## Visualization Requirements

Every figure must be generated programmatically (matplotlib or equivalent), saved as both PDF and PNG, and have its underlying data exported as JSON. Use a consistent style across all figures:

- Color scheme: use a colorblind-safe qualitative palette (e.g., Okabe-Ito or ColorBrewer Set2). Assign fixed colors to constructions: direct = blue, symmetric echo = orange, independent echo = green, asymmetric echo = red, hybrid = purple, new constructions = additional distinct colors. These assignments must be consistent across every figure in the report.
- Font: match revtex4-2 defaults (Computer Modern or equivalent).
- Axis labels: always include units. Always label both axes.
- Error bars or spread indicators on any aggregated metric.
- Figure size: single-column (3.375 in) or double-column (7 in) to match revtex4-2.

### Required Figures (Minimum Set)

These must all appear in the final report. You may add more, but these are non-negotiable:

1. **Unified prior-work comparison** (Part I): Bar chart of best F_avg by study and construction, with the new study's results overlaid.

2. **Per-manifold error budget** (Part III, Sec 7): Grouped bar chart, one group per manifold n, bars for |Îµ_x|, |Îµ_y|, |Îµ_z|, faceted by construction. Show for the best case of each construction type.

3. **Direct vs echoed fidelity scatter** (Part III, Sec 8): F_avg(best non-direct) on y-axis vs F_avg(direct) on x-axis, one point per case. Diagonal parity line. Color by construction type.

4. **Refocusing pulse manifold dependence** (Part III, Sec 8): Realized X_Ï€ rotation angle vs Fock number n, two curves: baseline Gaussian X_Ï€ and upgraded robust X_Ï€.

5. **Objective landscape for Î·** (Part III, Sec 8): L_echo(Î·) for 4 representative cases, showing the landscape shape around Î· = 0.

6. **Fidelity vs duration scaling** (Part III, Sec 9): F_avg vs |Ï‡|T/2Ï€ for each construction, with shading for case spread. Include the 0.99 target line.

7. **Ansatz comparison at fixed duration** (Part III, Sec 9): For the longest duration, bar chart of F_avg by ansatz variant (baseline Gaussian, + sidebands, + alt envelope, + GRAPE, hybrid). Error bars from multi-seed runs.

8. **Bloch sphere trajectories** (Part III, Sec 7): For the best direct case, one subplot per active manifold showing the realized vs ideal qubit trajectory. 2D Mollweide or stereographic projection is acceptable.

9. **Robustness spider plot** (Part III, Sec 10): For the top 3 constructions, a radar/spider chart with axes for F_unitary, F_decoherence, Ï‡-sensitivity, Ï‰_q-sensitivity, Î©-sensitivity. Normalized so that larger area = better.

10. **N_active scaling** (Part III, Sec 10): Best F_avg vs N_active for each construction family, with extrapolation curves.

11. **Spectral crowding diagram** (Appendix): Frequency axis with vertical lines for each manifold transition, shaded bands for tone bandwidth, showing where overlap occurs.

12. **Representative waveforms** (Appendix): Time-domain I/Q waveforms for the best construction at each N_active, showing the pulse structure clearly.

13. **Full case Ã— construction heatmap** (Appendix): Extended version of the baseline's Fig 3, now including all new constructions and all new grid points.

14. **Convergence validation** (Appendix): F_avg under baseline dt vs refined dt and baseline N_cav vs extended N_cav, for every new best case.

15. **Parameter count vs fidelity** (Part III, Sec 9): Scatter plot of ansatz parameter count vs best achievable F_avg, showing the tradeoff between model complexity and performance.


## Output Requirements

### File Structure

```
studies/definitive_ideal_sqr_gate/
â”œâ”€â”€ report/
â”‚   â”œâ”€â”€ main.tex                    # The complete self-contained report
â”‚   â”œâ”€â”€ references.bib
â”‚   â””â”€â”€ figures/                    # All PDF figures referenced by main.tex
â”œâ”€â”€ data/
â”‚   â”œâ”€â”€ prior_studies/              # Loaded results from all 3 prior studies
â”‚   â”‚   â”œâ”€â”€ study1_results.json
â”‚   â”‚   â”œâ”€â”€ study2_results.json
â”‚   â”‚   â””â”€â”€ study3_baseline_results.json
â”‚   â”œâ”€â”€ new_results/
â”‚   â”‚   â”œâ”€â”€ full_grid_results.csv   # One row per (case, construction, seed)
â”‚   â”‚   â”œâ”€â”€ best_per_case.csv       # One row per (case, construction), best seed
â”‚   â”‚   â”œâ”€â”€ construction_summary.json
â”‚   â”‚   â”œâ”€â”€ error_decomposition.json  # Per-manifold Pauli errors for every row
â”‚   â”‚   â”œâ”€â”€ sensitivity_analysis.json
â”‚   â”‚   â”œâ”€â”€ scaling_analysis.json
â”‚   â”‚   â””â”€â”€ fundamental_limits.json
â”‚   â””â”€â”€ validation/
â”‚       â”œâ”€â”€ convergence_checks.json
â”‚       â””â”€â”€ sanity_checks.json
â”œâ”€â”€ figures/
â”‚   â”œâ”€â”€ pdf/                        # Publication-quality PDFs
â”‚   â”œâ”€â”€ png/                        # Quick-inspection PNGs
â”‚   â””â”€â”€ data/                       # JSON data behind every figure
â”œâ”€â”€ artifacts/
â”‚   â”œâ”€â”€ cases/                      # Per-case JSON with optimized params
â”‚   â”œâ”€â”€ waveforms/                  # Sampled waveforms (.npz)
â”‚   â””â”€â”€ refocusing_pulses/          # Standalone refocusing pulse artifacts
â”œâ”€â”€ scripts/
â”‚   â”œâ”€â”€ run_full_study.py           # Master script, calls all objectives
â”‚   â”œâ”€â”€ run_objective_N.py          # Per-objective scripts (N = 1..7)
â”‚   â”œâ”€â”€ parallel_utils.py           # Shared parallelization infrastructure
â”‚   â”œâ”€â”€ diagnostics.py              # Error decomposition functions
â”‚   â”œâ”€â”€ plotting.py                 # All figure generation
â”‚   â”œâ”€â”€ load_prior_studies.py       # Loads and normalizes prior study data
â”‚   â””â”€â”€ validate_results.py         # Convergence and sanity checks
â””â”€â”€ README.md                       # How to reproduce everything
```

### Per-Objective Deliverables

For every objective, produce:

1. A structured summary: `{objective: int, tasks_completed: [str], headline_result: str, best_fidelity: float, comparison_to_baseline: str, figures_generated: [str], wall_clock_seconds: float}`.
2. All figures as PDF + PNG + JSON data.
3. Rows appended to `full_grid_results.csv`.
4. Updated LaTeX sections in `main.tex`.

### The Master Results Table

The single most important deliverable is a comprehensive ranked table of **every** construction tested across the entire study (prior + new), sorted by F_practical (the combined figure of merit from Objective 5). Columns:

| Rank | Construction | Ansatz | Envelope | Duration (ns) | N_active | F_unitary | F_decoherence | F_practical | Mean |Îµ_z| | Mean |Îµ_âŠ¥| | Params | Study |
|------|-------------|--------|----------|---------------|----------|-----------|---------------|-------------|---------|---------|--------|-------|

This table must appear in the main text and also be exported as CSV. It is the definitive ranking that answers "what is the best way to build an SQR gate?"


## Constraints and Guardrails

- **Load, do not re-run, prior study results.** Their data files are the ground truth baselines. Your new results extend them; they do not replace them.
- **Use the public `cqed_sim` API wherever possible.** For every local wrapper, include a comment block explaining why the public API is insufficient and what API extension would eliminate the need for the wrapper.
- **Preserve the baseline 16-case grid** as a strict subset of your expanded grid. Every baseline case must have a matched new-study counterpart so that improvement can be measured case-by-case.
- **Same loss function weights** (0.15, 1.0, 0.35, 0.35) for the direct optimizer as default. If you test alternative weights, do so as a clearly labeled ablation, not as a silent change.
- **Convergence checks are mandatory** for every new best case: dt halved and N_cav increased by 1. Changes must be below 5 Ã— 10^{-4} in fidelity.
- **8 random seeds minimum** for any new ansatz optimization. Report best, median, and standard deviation across seeds.
- **Be honest about negative results.** If a proposed improvement fails, present it with the same rigor and visual quality as a success. Negative results that are well-characterized are more valuable than poorly-characterized positive results. The final report should leave the reader with a clear map of "what works, what doesn't, and why."
- **No orphan claims.** Every claim in the report must be backed by either a figure, a table entry, or an explicit equation. If you state "the error is dominated by transverse components," there must be a figure showing this.


## Success Criteria (Tiered)

### Tier 1 (Landmark):
- F_avg > 0.99 on at least one physically credible case (N_active â‰¥ 2, realistic parameters, passes convergence checks).
- Clear identification of the ansatz and parameter regime that achieves this.

### Tier 2 (Strong Advance):
- F_avg > 0.95 on at least one case.
- F_avg > 0.85 as a mean across the baseline 16-case grid.
- The error floor is quantitatively explained (per-manifold decomposition shows exactly where the remaining 5% goes).

### Tier 3 (Meaningful Progress):
- F_avg > 0.90 on at least one case (up from 0.72454 baseline).
- A new construction or ansatz is identified that qualitatively outperforms the baseline direct waveform.
- The robustness analysis identifies a practically superior construction that was not the unitary-fidelity winner.

### Tier 4 (Definitive Negative Result):
- Convincing evidence that F_avg > 0.95 is unreachable with qubit-only multitone drive in the dispersive model, backed by:
  - GRAPE / unconstrained upper bound showing the physical limit.
  - Quantum speed limit analysis showing duration constraints.
  - Per-manifold error floor that is irreducible within the model.
- A concrete proposal for what control modification (e.g., cavity drive, two-channel drive, non-dispersive regime) would be needed to break through.

Any tier is a successful study. Tier 4 is *not* a failure â€” it is the most scientifically valuable outcome if it is rigorous.


## Execution Order

1. **Objective 1** (Diagnostics) â€” build the analysis pipeline.
2. **Objective 7** (Fundamental limits) â€” establish the ceiling before optimizing.
3. **Objective 2** (Refocusing pulse) â€” highest-leverage single fix. 
4. **Objective 3** (Direct ansatz expansion) â€” addresses the dominant error. 
5. **Objective 4** (Hybrid) â€” only if 2 and 3 show partial improvement.
6. **Objective 5** (Robustness) â€” apply to the best constructions from 2â€“4.
8. **Report assembly and final figure generation.**

Total target: â‰¤ 2 hours wall-clock on 8 cores.

If Objective 7 reveals that the physical upper bound is already < 0.95, reprioritize: skip the hybrid (Objective 4) and invest the time in a thorough Tier 4 analysis with additional random targets and alternative model variants.

## Final Mandate

When you are done, a reader should be able to pick up your report and:

1. Understand what an ideal SQR gate is and why it matters â€” without consulting any other document.
2. See every prior attempt and its limitations â€” presented fairly, with the original data, in your unified framework.
3. Follow the logic of every new approach you tried â€” with enough mathematical and computational detail to reproduce it.
4. See exactly how well each approach worked â€” in matched comparisons with consistent metrics and visualizations.
5. Understand *why* each approach worked or failed â€” through per-manifold error decomposition, not just aggregate numbers.
6. Know whether F_avg > 0.99 is achievable â€” and if not, know exactly what the obstruction is and what it would take to overcome it.
7. Know which construction they should use if they need an SQR gate tomorrow â€” ranked by practical viability, not just unitary fidelity.

Leave no important question unanswered. This is the definitive study.
           |

---

## Agent to Invoke

**Agent:** science-director  [model: Codex 5.4 xHigh]

Paste this into Copilot Chat:

`
@science-director study=studies/the_definitive_ideal_sqr_gate_study run=task_runs/the_definitive_ideal_sqr_gate_study phase=plan
`

---

## Context Files (read these first)

The agent MUST read these files before doing anything:

1. `studies/the_definitive_ideal_sqr_gate_study/study_state.json` - machine-readable study state
2. `task_runs/the_definitive_ideal_sqr_gate_study/SCIENCE_DIRECTIVE.md` - last planning directive
3. `task_runs/the_definitive_ideal_sqr_gate_study/EXECUTION_SUMMARY.md` - last execution results (if exists)
4. `task_runs/the_definitive_ideal_sqr_gate_study/TASK_CHECKLIST.md` - which tasks are done / open
5. `task_runs/the_definitive_ideal_sqr_gate_study/BLOCKERS.md` - known blockers
6. `studies/the_definitive_ideal_sqr_gate_study/IMPROVEMENTS.md` - improvement log

---

## Recovery Instructions for the Agent

This is a **recovery invocation** for a study that was interrupted or needs continuation.

- Current phase: **needs-plan**
- Iteration: 0 (max allowed: 14)
- Report preservation: ENABLED (append_iteration_section - do NOT overwrite report.tex)

Do NOT restart from scratch. Do NOT redo completed tasks. Read the context files
above, orient yourself to the current state, and continue the study forward.

## Active Blockers
- None

## Open Tasks (first 10 of 1)
- [ ] B0.2 Science Director produces first SCIENCE_DIRECTIVE.md


---

## Loop Configuration (from research_config.json)

- Planning model:   Codex 5.4 xHigh
- Execution model:  claude-opus-4-6
- Max iterations:   14
- Blocked policy:   continue_with_partial
- Report mode:      append_iteration_section

---
*Generated by tools/research_loop.ps1 -Action recover -StudyName the_definitive_ideal_sqr_gate_study*