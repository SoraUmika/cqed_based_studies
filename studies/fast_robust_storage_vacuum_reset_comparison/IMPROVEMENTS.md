# Improvement Log: Fast and Robust Active Cooling / Vacuum Reset in a Transmon-Storage-Readout cQED System

> Written for future agents. Be specific, honest, and actionable.

## Paper Review Findings (2026-04-21) — cleanness audit of `report/report.tex`

### P0 — BLOCKING: the paper is materially wrong in places, do not circulate

- **[P0 | CRITICAL] The stated robustness-score formula does not match the code.**
  - Eq. (2) of the report declares `R = 0.5 * (mean P_g00^coh + mean P_g00^cal)`.
  - Code at `scripts/run_study.py:1349–1378` computes `R = 0.5 * [1/(1 + mean(n_s + P_tr,exc))_coh + 1/(1 + mean(n_s + P_tr,exc))_cal]`.
  - Verified numerically against `data/coherence_robustness.csv` and `data/calibration_robustness.csv`:
    - Reported Table II: Pulsed 0.969, Bright 0.965, Raman 0.982.
    - Code formula reproduces: 0.9693, 0.9652, 0.9816 (matches).
    - Stated formula gives: 0.9510, 0.9630, 0.9809 (does NOT match, and **flips the pulsed-vs-bright ordering**).
  - Consequence: the paper's "Pulsed beats Bright on robustness (0.969 > 0.965)" is an artifact of a formula that is not the one stated. Under Eq. (2) as written, Bright beats Pulsed.
  - **Action:** either (a) rewrite Eq. (2) to the code formula with a physical justification for summing `n_s + P_tr,exc` (note: this mixes two qualitatively different failure modes — residual storage photons vs intermediate transmon population — and compresses them through a Lorentzian-like `1/(1+x)`), or (b) rerun all tables/figures under the stated P_g00-mean formula and accept the reordering. Cannot leave both.

- **[P0 | CRITICAL] "Within factor of ~40 of the autonomous benchmark" (report.tex:353) is numerically and directionally wrong.**
  - 5.6e-4 / 2.1e-5 = 26.7, not ~40.
  - More importantly, 2.1e-5 < 5.6e-4 means the 4× κ_r Raman residual is **below** the "autonomous benchmark," contradicting the paper's claim (§2.2) that the benchmark sets an upper bound.
  - **Action:** either replace the benchmark Γ_cool with a physically motivated value (e.g., weak-drive Purcell limit `Γ_cool = 4g_eff²/κ_r`), or drop the "upper bound / limiting case" framing and call it a "hand-tuned autonomous reference" whose rate is an input, not a limit.

- **[P0 | HIGH] The autonomous-benchmark Γ_cool is never stated and is not physically motivated.**
  - Scheme 4 in §2.2 describes `L ∝ a_s` but the rate is hidden. From Table II one infers `Γ_cool ≈ 1/535 ns ≈ 1.87 MHz`. This is smaller than κ_r/4 ≈ 1.04 MHz × π / 4 or the achievable two-photon coupling; the "benchmark" is weaker than what a real reservoir-engineered protocol could reach.
  - **Action:** state Γ_cool explicitly, and either justify it from device constraints or replace with a parameter scan over Γ_cool.

### P1 — Significant inconsistencies and misleading framing

- **[P1 | HIGH] "τ_e = 11 ns" mislabels the pulsed ladder's practical speed.**
  - `baseline_protocol_duration_ns` in `data/scheme_summary.csv` is 1000 ns for the pulsed ladder. The 11 ns number is the transfer-pulse population-decay timescale during one sub-segment, not the time to remove the photon.
  - The abstract's claim "fastest practical protocol (τ_e = 11.0 ns)" overstates by two orders of magnitude. The bar chart in `figures/scheme_comparison_summary.pdf` amplifies this by plotting 11 against 243.5 and 535 ns side-by-side without showing the 1000 ns ringdown.
  - **Action:** either add a column for `protocol_duration_ns` in Table II (so reader sees 1000 vs 4000 ns), or redefine τ_e globally as time for `⟨n_s(t)⟩` to drop by 1/e of its initial value (which will be much larger than 11 ns for the pulsed protocol).

- **[P1 | HIGH] Robustness grids are run at Δt = 1.0 ns; headline numbers at Δt = 0.5 ns.**
  - The paper admits this (end of §4.2) but then reports R to three decimal places. Comparing coh-grid baseline (T1×1, T2×1) P_g00 ≈ 0.96 (pulsed) to headline 0.969 shows the time-step bias is 0.5–1 pp — **the same size as the inter-scheme R gaps**.
  - The entire 0.013 Raman-vs-Pulsed margin is within numerical error. Either re-run the grids at Δt = 0.5 ns before quoting 3-digit scores, or publish error bars and demote the ranking to "within numerical uncertainty."

- **[P1 | MEDIUM] κ_r reduction factor quoted inconsistently: abstract says ~25, §3.3 says ~770.**
  - Abstract describes "1× to 4×" (actual ratio 0.000528 / 0.000021 ≈ 25).
  - §3.3 and Fig. 3 caption describe "0.5× to 4×" (actual ratio 0.016246 / 0.000021 ≈ 774).
  - Both factually correct but attached to different sweep ranges; reader cannot tell which is referenced. Figure shows 0.5×–4× data.
  - **Action:** pick one baseline (1× makes more physical sense as "the device") and use it everywhere.

- **[P1 | MEDIUM] Orphan figure.** `report.tex:511` defines `\label{fig:summary}` for `scheme_comparison_summary.pdf` with a full caption, but `\ref{fig:summary}` never appears in the body. Cite it in §3.2 or delete.

- **[P1 | MEDIUM] Bright-state transmon population: Table II `P_{tr,exc} = 0.036` vs §4.1 `max P_f ≈ 0.50`.**
  - Table column is end-of-window; §4.1 figure quotes peak. Off by a factor of 14. Column symbol `P_{tr,exc}` is never defined in text. Add `^{(\text{final})}` or a separate `P_f^{\max}` column.

- **[P1 | MEDIUM] Raman-like baseline τ_e drifts across the paper: 243.5 ns (§3.1, Table II), 246 ns (Table IV, Δt=1), 242 ns (§7 open questions).** Pick one.

- **[P1 | MEDIUM] Coherence grid's best-case P_g00 is below the headline P_g00.**
  - For pulsed, max coh-grid P_g00 = 0.9639 at (T1×2, T2×2); headline `baseline_ground_vacuum_final` = 0.9694. A reader expects that relaxing noise can only help, so the grid best-case should equal or exceed the baseline. The gap is a Δt-bias artifact (grid at 1 ns vs baseline at 0.5 ns). Flag this in the caption or eliminate by matching time steps.

### P2 — Missing equations and figures that a reader will expect

- **[P2 | HIGH] No drive Hamiltonians written.** Eq. (1) ends with `H_drive(t)` as a token. At minimum add:
  - Pulsed: `Ω_{gf}^{(n)}(t) [σ_{fg} e^{-iω_d t} + h.c.]` with n-resolved carrier ω_d = ω_{gf} - n·χ_s.
  - Continuous sideband: `Ω_sb [a_s^† a_r σ_{fg}^† e^{-iω_sb t} + h.c.]` (effective; flag as such).
- **[P2 | HIGH] No Raman-effective-coupling or Purcell-damping equations.** Adiabatic elimination: `g_eff ≈ Ω_1 Ω_2 / (2Δ)` and `Γ_cool ≈ 4 g_eff² / κ_r` in the bad-cavity limit. **This directly answers §7's "open question"** about κ_r vs τ_e (Γ_cool ∝ 1/κ_r means bigger κ_r → slower τ_e, exactly what `data/readout_kappa_tradeoff.csv` shows: 242 → 327 ns as κ_r goes 0.5× → 4×).
- **[P2 | HIGH] Lindblad master equation is never written.** Appendix B just names it. Writing it explicitly gives every rate in Table I a home.
- **[P2 | MEDIUM] τ_e and P_g00 are never defined quantitatively.** Current text assumes the reader knows. Add: `τ_e` = exponential-fit 1/e time on `⟨n_s(t)⟩` over `[0, t_fit]`; `P_g00 = ⟨g,0_s,0_r | ρ | g,0_s,0_r⟩`.
- **[P2 | MEDIUM] No device schematic or level diagram.** Standard in circuit-QED papers; prose-only description of four mechanisms is unnecessarily dense. A two-panel schematic (modes + couplings; level ladders with bright vs detuned arrows) would shorten §2.2 significantly.
- **[P2 | MEDIUM] No Pareto / τ_e-vs-R scatter plot.** The central speed-robustness message of §6.2 lives only in prose. Four labeled points on a log-τ_e vs R plot carries the paper.
- **[P2 | MEDIUM] No κ_r-vs-τ_e overlay.** Fig. 3 apparently plots only fidelity/residual. Adding τ_e(κ_r) on a twin axis visualizes the speed-fidelity tradeoff that §7 flags as an open question.
- **[P2 | MEDIUM] No convergence figure.** Table IV's `\makebox` formatting is fragile. Replace with a log–log plot of `⟨n_s⟩_f` vs Δt per scheme; the 2 ns breakdown becomes visually obvious.
- **[P2 | LOW] No transient vs steady transmon trace.** A zoomed P_f(t) panel for the bright-state would resolve the 0.50-vs-0.036 confusion at a glance.

### P3 — Methodology and scope concerns

- **[P3 | HIGH] The code's R formula mixes `n_s` and `P_tr,exc` additively.** A stray transmon e/f population (which will decay back through g with χ_s imprinting storage phase) is not the same failure mode as a residual storage photon. They should not be summed 1-to-1 unless justified. If using a composite, at least argue the relative weight.
- **[P3 | MEDIUM] Convergence tested only for `|1⟩` initial state.** Table III reports |3⟩, |α=1⟩, thermal results at Δt = 0.5 ns without re-verifying convergence for the larger Fock sector. Fock-5 truncation at |3⟩ is plausibly adequate but untested.
- **[P3 | MEDIUM] "Larger" truncation adds only +1 Fock level per mode.** Weak convergence lever; should test +2 and +3 to see a trend, not a single step.
- **[P3 | MEDIUM] Equal weighting of coherence and calibration in R.** Not motivated. Calibration drift is recoverable (re-tune); hardware coherence is not. A realistic composite would weight them differently, or report separately.
- **[P3 | MEDIUM] The "25 grid points per scheme" wording is ambiguous.** "T1 and T2 **jointly** scaled by factors {0.5,...,2.0}" plus "5×5 grid" is contradictory on a first reading (joint → diagonal → 5 points; 5×5 → 25 points). Clarify: "T1 scaled by {...} independently of T2 scaled by {...}."
- **[P3 | LOW] Table IV's `\makebox[2.1cm][r]{$x$; $τ$}` trickery is fragile across journals.** Replace with two sub-columns per truncation group.
- **[P3 | LOW] Notation drift.** `T_{2,R}` in Table I vs bare `T_2` in body; `\eqref{eq:robustness}` in one caption vs `(Eq.~\ref{eq:robustness})` in another; `Fig.~\ref` vs `Figure~\ref`. Pick one style each.

## Critical Gaps (P1)
- **[P1 | HIGH]** The continuous schemes are still driven with constant square overlaps only: the present comparison shows the right qualitative speed-versus-robustness ordering, but it does not yet test STIRAP-like timing or open-system optimal control.
- **[P1 | HIGH]** The autonomous `L \propto a_s` result remains an auxiliary benchmark rather than a native `cqed_sim` replay.
- **[P1 | MEDIUM]** The continuous-scheme end-of-window residual occupations remain sensitive to the time step: the manuscript now quotes the final-state numbers from the `0.5 ns` baseline and uses the `1.0 ns` sweeps only for robustness ranking, but a tighter convergence campaign is still needed before claiming fully converged absolute residuals.

## Recommended Improvements (P2)
- **[P2 | MEDIUM]** Add pump-aware Stark-shift and parasitic-channel modeling on top of the present effective sideband layer.
- **[P2 | MEDIUM]** Extend the continuous schemes to multi-tone or shaped photon-number-aware driving so higher-Fock cooling is not limited by the `n=1`-centered line choice.
- **[P2 | MEDIUM]** Re-run the Raman-like protocol with engineered larger readout linewidth, since the current device linewidth is not yet fully in the autonomous bad-cavity regime.

## Nice-to-Haves (P3)
- **[P3 | LOW]** Add explicit measurement-backaction and readout-heating models if those become relevant on hardware.

## Open Questions
- How close can a shaped counter-intuitive two-tone protocol get to the autonomous benchmark while staying faster than the present readout lifetime?
- Does the present best pulsed protocol remain best once the sideband controls are embedded in a microscopic pump model?
- What readout linewidth increase is required before the Raman-like protocol becomes clearly superior on both speed and robustness?

## What Was Tried and Did Not Work
- **Constant resonant continuous driving as a generic multi-photon solution**: it remains fast for `n=1` but cools higher-Fock support much less cleanly because the real transmon path is heavily occupied and the fixed carrier is still centered on the lowest manifold.
- **Assuming the most virtual detuned protocol is automatically best on the current device**: on the present readout linewidth it is more robust to transmon coherence, but still slower and less complete than the pulsed ladder within the same wall-clock window.

## Compute & Resource Notes
- Main comparative run: `2510.6 s`
- Continuous candidate scans reused earlier validated pulse winners instead of re-running a global waveform search.

## Resolved
- **Metric-definition mismatch across summary artifacts**: `scheme_summary.csv`, `study_results.json`, and the report now use end-of-run `final_*` values consistently for headline residuals, while tail-averaged `steady_*` values remain available only as diagnostic fields.
- **Pulsed initial-state ladder-depth mismatch**: the pulsed initial-state comparisons now use a matched ladder depth (`n=1` for $|1\rangle$, `n=3` for $|3\rangle$, and the full available ladder for coherent and thermal states) instead of forcing a four-rung sequence for every input state.
- **Single-photon summary versus initial-state timestep mismatch**: the initial-state comparison artifact now uses the same `0.5 ns` baseline step as the headline single-photon summary, so the single-photon rows agree across the report and saved tables.
