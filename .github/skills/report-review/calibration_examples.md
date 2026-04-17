# Evaluator Calibration Examples

> Reference graded examples for each review dimension. Use these to calibrate
> your scoring. A score of 5 is publication-ready; a score of 2 means major rework.

---

## Dimension A: Writing Quality — Score Examples

### Score 5 (Publication-Ready)
> "We investigate the dispersive readout of a transmon qubit coupled to a
> coplanar waveguide resonator. The dispersive shift $\chi/2\pi = -2.84\,\text{MHz}$
> is extracted from steady-state transmission simulations using the Jaynes--Cummings
> Hamiltonian in the dispersive limit ($g/\Delta = 0.04 \ll 1$). Section~II introduces
> the system model and defines all parameters. Section~III presents the optimization
> of readout pulse duration and amplitude. Section~IV validates convergence against
> Hilbert space truncation and compares the extracted $\chi$ to the perturbative
> prediction."

**Why score 5:** Self-contained abstract. Every symbol defined at first use with units. Clear roadmap. No code identifiers. Approximation validity stated quantitatively.

### Score 3 (Significant Improvement Needed)
> "We simulate a transmon-cavity system using `cqed_sim`. The chi shift is computed
> and optimized. Results are shown in Fig. 3. The fidelity is good."

**Why score 3:** Code identifier in prose (`cqed_sim`). No symbol definitions. "Good" is not quantitative. No roadmap. Figure referenced without interpretation. Missing units on chi.

### Score 2 (Major Rework)
> "run_opt.py generates the results. The optimization converged. See data/results.npz
> for the full output."

**Why score 2:** Filenames in prose. No physics content. No interpretation. Reads as a README, not a scientific paper.

---

## Dimension B: Evidence-Claim Mapping — Score Examples

### Score 5 (Every Claim Backed)
> **Claim:** "The optimized GRAPE pulse achieves gate fidelity $\mathcal{F} = 0.9987$
> at $T = 200\,\text{ns}$, a $3\times$ improvement over the naive square pulse
> ($\mathcal{F} = 0.9960$)."
>
> **Evidence:** Figure 4 shows both pulses' fidelity vs. time with labeled curves.
> Table II gives exact values. Appendix A provides the optimized waveform.
>
> **Verdict:** SUPPORTED — quantitative claim with direct figure, table, and artifact support.

### Score 3 (Weakly Backed)
> **Claim:** "The gate performance is significantly improved by GRAPE optimization."
>
> **Evidence:** Figure 3 shows a convergence trace of the cost function.
>
> **Verdict:** WEAK — "significantly improved" requires a comparison baseline. The convergence trace shows the optimizer ran, not that the result is good. Need: fidelity before vs. after, with the comparison baseline shown on the same figure.

### Score 2 (Unsupported)
> **Claim:** "Our approach achieves near-optimal performance for this system."
>
> **Evidence:** None cited. No comparison to theoretical bounds or alternative methods.
>
> **Verdict:** UNSUPPORTED — "near-optimal" requires evidence: either multiple random restarts with consistent convergence, a comparison to an analytic bound, or a landscape analysis. Hard-reject trigger.

---

## Dimension C: Physics and Methodology — Score Examples

### Score 5 (Rigorous)
> "We work in the dispersive regime with $g/\Delta = 0.04$. The second-order
> correction to the dispersive shift is $\delta\chi/\chi \sim (g/\Delta)^2 \approx
> 0.16\%$, negligible compared to our target accuracy of $1\%$. Hilbert space
> convergence is confirmed: increasing the Fock cutoff from $N = 8$ to $N = 15$
> changes the extracted $\chi$ by $< 10^{-4}\,\text{MHz}$ (Table III). Five
> independent GRAPE runs from random initializations converge to fidelities
> within $\pm 2 \times 10^{-4}$ of each other (Fig. 7), providing evidence
> against local-minimum trapping."

**Why score 5:** Approximation validity quantified. Convergence shown with numbers. Multiple restarts address global optimality. Specific tolerances stated.

### Score 3 (Asserted Without Verification)
> "We use the dispersive approximation, which is valid for our parameters.
> Convergence was verified by increasing the Hilbert space dimension.
> The optimization was run once from a physically motivated initial guess."

**Why score 3:** Validity asserted but not quantified (what is $g/\Delta$?). Convergence "verified" but no numbers. Single optimizer run — local minimum not addressed.

### Score 2 (Fundamental Gap)
> "Parameters are chosen to be in the dispersive regime. QuTiP solves the
> master equation."

**Why score 2:** No validation at all. No convergence check. No approximation verification. No uncertainty quantification. Black-box simulation with no sanity checks.

---

## Dimension E: Novelty — Score Examples

### Score 5 (Clear Advance)
> "Prior work [Koch07] established the dispersive shift for the two-level
> transmon approximation. We extend this to the full multi-level transmon
> with $N = 5$ levels and show that the higher-level contribution shifts
> $\chi$ by $-3.2\%$ at $\alpha/2\pi = -255\,\text{MHz}$ — a correction
> that is experimentally resolvable with current readout sensitivity.
> This quantitative correction has not been reported previously."

**Why score 5:** Clearly states what existed before (two-level result), what is new (multi-level correction with specific number), and why it matters (experimentally resolvable).

### Score 3 (Incremental/Unclear)
> "We compute the dispersive shift using cqed_sim and find good agreement
> with the known formula."

**Why score 3:** Reproducing a known result without adding new insight. If the study is REP-class, it needs to explain what specifically is being validated and what new understanding the reproduction provides.

### Score 2 (No Contribution)
> "We simulate the transmon-cavity system and report the parameters."

**Why score 2:** No new insight. No comparison. No advance beyond running existing code with standard parameters.

---

## Using These Examples

When scoring a study:
1. For each dimension, find the calibration example closest to what you see
2. Assign the corresponding score
3. If between two examples, choose the lower score and explain what would raise it
4. Any dimension at score 2 with a blocking issue triggers NEEDS_REWORK
5. Any dimension at score 3 with blocking issues triggers REVISE
6. APPROVE requires all dimensions at 4 or above with no blocking issues
