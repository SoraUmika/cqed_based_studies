# Science Directive

Date: 2026-04-10
Study: `studies/gf_sideband_waveform_optimization`
Run: `task_runs/gf_sideband_waveform_optimization`

## Problem Class
OPT, ANA, DES

## Central Question
What waveform family is best for implementing the effective `|f,n-1> <-> |g,n>` sideband interaction in the local `cqed_sim` device model, when the answer is separated into:

1. selective sideband control, and
2. fastest practical unselective sideband control,

for both the storage mode and the readout mode?

## Scope
1. Use the actual local editable `cqed_sim` sideband-reset device parameters.
2. Study both storage-sideband and readout-sideband red-sideband interactions for `n = 1,2,3`.
3. Compare at least the requested waveform families:
   - square,
   - Gaussian,
   - cosine,
   - flat-top cosine,
   - flat-top Gaussian,
   - smooth compact-support ramp,
   - and one additional physically motivated window family.
4. Optimize reasonably within each family over duration, amplitude, and the relevant shape parameter.
5. Define explicit selective and unselective metrics, then extract the fastest durations satisfying those metrics.
6. Quantify transfer, leakage, projected-subspace phase quality, robustness, and neighboring-manifold selectivity.
7. Include the available open-system channels in a finalist-only follow-up pass.

## Non-Negotiable Modeling Boundaries
- Use `cqed_sim` as the primary simulation engine.
- Treat the sideband interaction as an effective control operator, not as a microscopic pump derivation.
- Base the waveform ranking on the full multilevel driven replay, not on a reduced two-state model alone.
- State clearly which conclusions are robust within the effective model and which require a future pump-aware extension.

## Required Deliverables
- Mode-resolved dressed-frequency tables
- Waveform-comparison sweep results and finalist tables
- Speed-versus-selectivity tradeoff plots
- Robustness maps and representative population traces
- Machine-readable results for every headline conclusion
- Markdown summary, LaTeX report/PDF, notebook, and handoff files
