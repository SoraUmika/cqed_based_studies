# Improvement Log: SQR Gate Design

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **No noise-aware GRAPE optimization**: GRAPE controls were optimized in the closed system and then replayed under noise. Directly optimizing under a noisy objective could close the gap between GRAPE and parametric baselines. Requires `cqed_sim` support for noisy GRAPE objectives.
- **Phenomenological noise model**: The open-system replay uses aggregate T1/T2 parameters rather than a calibrated device-specific noise model. Particularly important for the three-mode concurrent readout results.

## Recommended Improvements (P2)
- **Reuse the literature-informed selective primitive baseline**: before launching new waveform scans, start from `studies/literature_informed_selective_primitives`, which now provides literature-backed Gaussian, cosine-squared, and flat-top Gaussian SQR/SNAP operating points under the repository's typical cQED parameters. This should replace ad hoc re-derivation of selective pulse widths and noisy starting guesses. [LOW difficulty]
- **Echoed/composite SQR sequences**: Can the residual spectator branch-local qubit-Z structure be reduced by echoed or two-layer selective pulses without sacrificing the favorable cavity-phase behavior? [MEDIUM difficulty]
- **Larger logical windows**: The combined cavity-phase compilation + SQR study used N=8. Test at N=12+ to assess scalability of the linear phase profile. [LOW]
- **Non-Fock encodings**: Test SQR families in cat-code or binomial-code logical subspaces where the Fock-selective advantage may change. [HIGH]
- **Realistic χ_sr coupling**: The three-mode study used χ_sr=0. Include physical storage-readout cross-Kerr to get realistic concurrent readout penalties. [MEDIUM]

## Nice-to-Haves (P3)
- **DAC quantization effects**: Map how finite DAC resolution affects the optimized correction vectors (d_λ, d_α, δω). [LOW]
- **Temperature-dependent noise budgets**: Sweep fridge temperature to map the thermal operating window for SQR. [MEDIUM]
- **Waveform bridge for GRAPE exports**: Build a calibrated pulse-export pipeline so GRAPE solutions can be hardware-replayed. [HIGH]

## Open Questions
- Can the corrected-parameterization multitone (Part III of waveform study) close the gap to GRAPE in the all-branch case?
- How do the optimized correction vectors correlate with the linearly-fitted cavity-only phase profile?
- Why does square-family outperform cosine-squared under realistic noise despite being suboptimal in the closed system at longer durations?
- Does the Purcell filter analysis change qualitatively for devices with different qubit-readout detunings?

## What Was Tried and Did Not Work
- **Simultaneous common-Gaussian multitone SQR** (Study 2): Target branch angle response is ~10⁻⁴ even for π-rotation targets. The pulse acts as near-identity on all blocks. Three correction strategies were tried on the hard case ({0,1}, θ=π, χT/2π=3): (1) common amplitude correction — flat response, (2) multistart amplitude/phase/detuning optimization — stuck near identity, (3) naive repeated π/8 step compilation — fails because each step also acts as near-identity. The common Gaussian multitone ansatz is fundamentally insufficient for simultaneous multi-branch SQR.
- **Cavity block-phase compilation for all-branch short-gate** (Study 1, Part II): At χT/2π=0.5, the fitted phase profile is irregular and applying it actually worsens strict logical fidelity (from 0.138 to 0.028). The structured all-branch short-gate failure is not a missing bosonic SNAP-like layer.
- **Active depletion of single-segment heuristic** (Study 1): Did not beat passive ring-down.

## Compute & Resource Notes
- Closed-system single-point evaluation: ~0.1–0.5 seconds
- Open-system multilevel replay: ~2–5 seconds per point (n_tr=3, n_cav=6)
- Three-mode replay: ~10–30 seconds per point
- GRAPE optimization (closed system): ~minutes per duration point; full sweep ~1 hour
- Full production scans: several hours total

## Resolved
- **No repository-level literature baseline for selective SQR/SNAP envelopes**: resolved by `studies/literature_informed_selective_primitives`, which now supplies optimized Gaussian, cosine-squared, and flat-top Gaussian references with noise-aware rankings.
