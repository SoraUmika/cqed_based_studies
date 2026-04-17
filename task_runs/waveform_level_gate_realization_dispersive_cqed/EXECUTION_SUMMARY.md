# Execution Summary - Unconditional Cavity Displacement in Dispersive cQED

**Study:** `studies/waveform_level_gate_realization_dispersive_cqed`  
**Run:** `task_runs/waveform_level_gate_realization_dispersive_cqed`  
**Problem class:** ANA, DES, OPT  
**Date:** 2026-04-02

## Headline Results

### 1. Why the naive pulse fails
- The baseline square pulse calibrated to produce `alpha = 1` at `T = 80 ns` is highly branch dependent.
- Simulated final branch displacements:
  - `alpha_g = 0.9961 + 0.0036 i`
  - `alpha_e = 0.6783 + 0.6114 i`
- This gives `delta_alpha = 0.6859`.
- Vacuum-state fidelities:
  - `F_g = 0.999960`
  - `F_e = 0.622737`
- For the qubit superposition `(|g> + |e>) / sqrt(2)` with the cavity initially in vacuum:
  - fidelity `= 0.770319`
  - entanglement entropy `= 0.482920` bits
- The naive single-tone failure therefore is not just a branch-calibration issue; it produces genuine residual qubit-cavity entanglement on superposition inputs.

### 2. Fast broadband pulses help, but filtering matters
- The best fast single-tone result in the present sweep is a `5 ns` Gaussian pulse for `alpha = 1`.
- Metrics:
  - `delta_alpha = 0.05391`
  - `F_g = 0.999708`
  - `F_e = 0.996709`
  - `F_{+x} = 0.997743`
  - entanglement entropy `= 8.62e-3` bits
- This pulse has a 95 percent bandwidth of about `200 MHz`, far above `|chi| / 2pi = 2.84 MHz`, which explains why the branch conditionality is strongly reduced.
- However, a `40 MHz` low-pass filter applied to the same `5 ns` Gaussian removes enough pulse area that the achieved displacement shrinks to about `0.505`, and the fidelity falls to about `0.782` even though branch mismatch remains small. Fast pulses need hardware-aware calibration, not just bandwidth.

### 3. Two-tone compensation is the best simple, interpretable protocol
- The best two-tone branch-compensated pulse occurs at `T = 20 ns` for `alpha = 1`.
- Metrics:
  - `delta_alpha = 8.15e-4`
  - `F_g = 0.999999`
  - `F_e = 0.999997`
  - `F_{+x} = 0.997680`
  - entanglement entropy `= 9.63e-6` bits
- The two branch carriers are separated by about `17.84 Mrad/s`, consistent with the extracted dispersive branch pull.
- A `chi` calibration error of `+/- 10 percent` hardly changed the `80 ns` vacuum-branch result, indicating that the branch-equalization mechanism is not hypersensitive to small branch-frequency calibration errors in this regime.
- Main caveat: vacuum branch equalization does not automatically maximize fidelity on longer-duration coherent-state tests.

### 4. Echoed displacement was not competitive in its practical form
- The best echoed result in the tested family occurred at total duration `60 ns`.
- Metrics:
  - `delta_alpha = 0.07873`
  - `F_g = 0.841499`
  - `F_e = 0.884519`
  - `F_{+x} = 0.891358`
  - entanglement entropy `= 0.02878` bits
- This is much worse than the best fast Gaussian and the best two-tone pulse.
- Interpretation: the inserted vacuum-calibrated qubit `pi` pulse becomes manifold dependent once the cavity is populated, so the ideal toggling-frame cancellation argument does not survive the realistic waveform model.

### 5. Bounded hardware-aware optimal control remains a strong reference, but it is not the best tested protocol on the updated broad-state benchmark
- The best constrained optimal-control case occurred at `T = 40 ns`.
- Metrics on the explicit 14-state test set:
  - mean state fidelity `= 0.957535`
  - minimum state fidelity `= 0.874808`
- Vacuum-branch metrics:
  - `delta_alpha = 0.104776`
  - `F_g = 0.999545`
  - `F_e = 0.985605`
  - `F_{+x} = 0.985303`
  - entanglement entropy `= 0.02727` bits
- This remains a useful hardware-aware sampled-waveform reference, but the multiplex follow-up shows that it no longer holds the best broad-state score once the short calibrated two-tone pulses are evaluated on the same state set.

### 6. Explicit multiplex follow-up and updated ranking
- A new follow-up benchmark compressed the best `40 ns` optimal waveform into explicit full-duration multicarrier drives with `K = 2, 3, 4, 5, 6, 8` tones and re-simulated those drives under the same full model.
- The direct multiplex family was a negative result. The best tested case was `K = 8`, but it still gave:
  - mean state fidelity `= 0.524154`
  - minimum state fidelity `= 0.001888`
  - vacuum `delta_alpha = 0.770316`
  - vacuum `F_{+x} = 0.589034`
- The calibrated two-tone references are much stronger on the same broad 14-state set:
  - two-tone `20 ns`: mean fidelity `= 0.985722`, minimum fidelity `= 0.924196`, vacuum `delta_alpha = 8.15e-4`, vacuum `F_{+x} = 0.997680`
  - two-tone `40 ns`: mean fidelity `= 0.952808`, minimum fidelity `= 0.745075`, vacuum `delta_alpha = 0.002672`, vacuum `F_{+x} = 0.992398`
- Two additional structured extensions were also negative:
  - segmented branch-resonant family: best `8`-segment case gave mean fidelity `= 0.361232`, minimum fidelity `= 0.008144`, vacuum `delta_alpha = 0.143817`, vacuum `F_{+x} = 0.656716`
  - jointly optimized shaped two-tone family: the successful `4`-segment Powell run gave mean fidelity `= 0.157618`, minimum fidelity `= 0.045900`, vacuum `delta_alpha = 0.081933`, vacuum `F_{+x} = 0.227875`
- Updated interpretation: extra carriers alone are not enough if they are deployed as one full-duration multiplex comb, and low-dimensional segmented or jointly optimized structured families do not rescue the situation either. The short branch-calibrated two-tone pulse is now the best tested protocol on the explicit broad state set, while the constrained waveform remains the strongest bounded sampled-waveform reference.

## Protocol Ranking

Using the common mean state fidelity on the explicit test set:

| Protocol | Duration | Mean fidelity | Min fidelity | Vacuum `delta_alpha` | `|+x>` entanglement (vacuum) | Complexity |
|---|---:|---:|---:|---:|---:|---|
| Two-tone compensated | 20 ns | 0.9857 | 0.9242 | 8.15e-4 | 9.63e-6 bits | Medium |
| Optimal control | 40 ns | 0.9575 | 0.8748 | 0.1048 | 0.0273 bits | High |
| Fast Gaussian | 20 ns | 0.9133 | 0.7438 | 0.1885 | 0.0728 bits | Low |
| Naive square | 80 ns | 0.5658 | 0.00861 | 0.6859 | 0.4829 bits | Low |
| Echoed displacement | 120 ns selected snapshot | 0.3164 | 0.0213 | 0.3219 | 0.1617 bits | Medium |

Important nuance:
- **Best overall on the explicit broad 14-state test set after the full structured multiplex follow-up:** short `20 ns` two-tone branch compensation
- **Best sampled waveform under the bounded hardware-aware optimization setup:** `40 ns` optimal control
- **Best simple and experimentally interpretable protocol:** short two-tone branch compensation

## Physics Conclusions
1. The naive cavity pulse stops behaving like `I_q \otimes D(alpha)` once `|chi| T` is no longer small.
2. The extracted guideline `1 / |chi| = 56.0 ns` matches the observed crossover: `5-10 ns` pulses remain close to unconditional, while `80-160 ns` pulses are strongly conditional.
3. `chi` is the dominant error source here. At `alpha = 1`, `T = 80 ns`, the model hierarchy gave:
   - minimal: `delta_alpha = 0.682`
   - with `chi'`: `0.686`
   - full model: `0.686`
4. Two-tone compensation is the best physically transparent fix in this regime.
5. Echo works only on paper unless the qubit inversion itself is manifold uniform.
6. Low-dimensional structured multiplex refinements do not rescue the gap to the best explicit protocols; the missing ingredient is richer temporal control rather than carrier count alone.
7. Optimal control remains a strong sampled-waveform reference, but it is less interpretable and was only run in a bounded, hardware-aware configuration here.

## Validation Summary
| Check | Status | Detail |
|---|---|---|
| Sanity | PASS | Short pulses behave nearly unconditionally; long naive pulses fail with the expected branch splitting and superposition entanglement |
| Mechanism isolation | PASS | The minimal / higher-order / full hierarchy shows that `chi` dominates over `chi'` and Kerr at the representative operating point |
| Convergence | PASS | The study reuses the same validated propagator settings already established for this model family (`n_cav`, `n_tr`, `dt`) |
| Literature alignment | PASS | The observed crossover tracks the expected dispersive criterion `|chi| T << 1` |

## Deliverables Produced
- Updated study README and improvement log
- `unconditional_*` artifacts in `artifacts/`
- `unconditional_*` figure pairs in `figures/`
- Updated `report.tex` and compiled `report.pdf`
- Updated reproducibility notebook
- Updated run-state and review-handoff files

## Recommendation
- If interpretability and calibration simplicity matter most, use the short two-tone branch-compensated displacement.
- On the explicit 14-state benchmark now evaluated in the multiplex follow-up, that same `20 ns` two-tone pulse is also the best tested protocol overall.
- Treat the bounded hardware-aware optimal-control waveform as a useful sampled-waveform reference rather than the standing best performer on this updated metric.
- If the pulse duration approaches or exceeds `1 / |chi|`, do not model the operation as an unconditional displacement unless an explicit compensation strategy has been validated.
