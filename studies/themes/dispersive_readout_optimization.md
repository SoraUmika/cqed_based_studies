# Theme: Dispersive Readout Pulse Optimization

> Unified synthesis of 4 component studies forming a progressive model-refinement chain.

---

## Component Studies

| Order | Study | Model Level | Key Contribution |
|-------|-------|-------------|-----------------|
| 1 | `readout_pulse_optimization` | Linear dispersive resonator | Analytic theory, adjoint GRAPE, optimal drive frequency |
| 2 | `procedural_readout_pulse_sequence_optimization` | Multilevel transmon + cavity (2-mode) | Procedural pulse families, bounds hierarchy, T1-limited reference |
| 3 | `nonlinear_qnd_hardware_realistic_readout` | 2-mode + hardware distortion + effective non-QND | Ring-hold family, QND stress tests, hardware-realistic frontier |
| 4 | `measurement_induced_leakage_ionization_modeling` | 2-mode + strong-readout disturbance + higher-ladder | Level-resolved leakage maps, ionization metric, regime thresholds |

---

## Unified Problem Statement

Dispersive cQED readout drives a resonator coupled to a transmon qubit, using the state-dependent frequency shift χ to distinguish |g⟩ from |e⟩ via the outgoing field. The optimization problem is fundamentally multi-objective: maximize state discrimination while minimizing residual photon occupation, measurement-induced state disturbance (QND defect), leakage to non-computational states, and total measurement duration.

This theme traces the full readout optimization story from the simplest analytic model through device-realistic simulation, answering the question: **which pulse shape, at what duration and amplitude, gives the best practical readout performance under realistic hardware and physics constraints?**

---

## Unified Methods

### Model Hierarchy

```
Level 0: Linear dispersive resonator (analytic ODE)
    ↓ adds multilevel transmon, cavity damping, qubit decoherence
Level 1: DispersiveTransmonCavityModel (n_tr=3, n_cav=14)
    ↓ adds HardwareConfig replay, multilevel relaxation, effective strong-drive mixing
Level 2: DispersiveTransmonCavityModel + SequenceCompiler(hardware=...) + mixing layer
    ↓ adds higher-ladder continuation, regime mapping over amplitude × duration × detuning  
Level 3: Level 2 with extended parameter sweeps and experiment-facing diagnostics
```

### Pulse Families Compared (across all studies)

| Family | Studies | Parameters | Best Regime |
|--------|---------|------------|-------------|
| Square (constant-amplitude) | 1, 2, 3 | 1 (amplitude) | Near-optimal in linear model |
| Smoothed square (cosine rise) | 2, 3 | 2 (amplitude, rise time) | Good QND at short durations |
| Ring-hold (ramp-up, hold, ramp-down) | 3 | 3–4 | **Best at all tested durations under rich model** |
| Segmented procedural (3–5 segments) | 2, 3 | 6–10 | Competitive on information fidelity; higher QND cost |
| Fourier basis | 2, 3 | 4–8 | Smooth; graceful degradation under amplitude stress |
| Nulling tail (analytic) | 2, 3 | 2–3 | Exact emptying in theory; clipping degrades in practice |
| Piecewise reference (high-dim) | 2, 3 | 16+ | Benchmark only — impractical residual photons |
| GRAPE (adjoint gradient) | 1 | 60 segments | ~10–20% above square in linear model only |

### Metrics

| Metric | Definition | Source |
|--------|-----------|--------|
| F_η (information fidelity) | Matched-filter assignment fidelity including detector efficiency | 2, 3 |
| Δ_MF (matched-filter separation) | Integrated I/Q field difference | 1 |
| n_res (residual photons) | Mean photon number at measurement end | 2, 3 |
| n_peak (peak photons) | Maximum photon occupancy during pulse | 2, 3, 4 |
| Q_QND | Repeated-readout state preservation | 2, 3 |
| P_leak | Total leakage to |f⟩ and above | 3, 4 |
| P_ion | High-manifold occupation (m ≥ 3) | 4 |
| Induced transition probability | Measurement-induced |g⟩↔|e⟩ transition | 3 |

---

## Unified Key Results

### Result 1: Linear Theory Remains Predictive for Drive Frequency but Not for Pulse Shape

Study #1 established that the optimal drive frequency is the midpoint ω_d = ω_r + χ/2 across all χ/κ ratios. This holds through all model levels. However, the linear-model conclusion that "square is near-optimal" breaks down once QND constraints and hardware distortion are included.

### Result 2: Procedural Pulses Are Competitive on Information but Pay QND Cost

Study #2 showed that 3–5 segment procedural pulses can match or exceed smoothed-square baselines on matched-filter fidelity (F_η = 0.9628 at 240 ns). But the family-independent QND defect in the 2-mode model (Q_QND ≈ 0.9933) masked any pulse-shape sensitivity.

### Result 3: Hardware + Non-QND Model Changes the Winner

Study #3 upgraded to hardware distortion + effective strong-drive mixing. The key finding: **ring-hold is the best family at all three tested durations** (96, 240, 496 ns), with matched-filter fidelities 0.6480, 0.9599, and 0.9982 respectively. Free procedural segments still reach high F_η but at visibly larger QND cost (Q_QND = 0.9886, induced transition 4.66e-3).

### Result 4: Strong Readout Has Sharply Defined Leakage Onset

Study #4 mapped the leakage regime: the representative 240 ns pulse crosses the visible-leakage benchmark (P_leak > 1e-4) at ε/(2π) = 5.5 MHz. The threshold shifts to 7.5 MHz at 120 ns and down to 4.5 MHz at 480 ns. High-manifold "ionization-like" occupation (P_ion) remains orders of magnitude below on-set leakage across the stable scan window.

### Result 5: Practical Readout Recommendation

| Duration | Best Family | F_η | Q_QND | P_leak | Notes |
|----------|-------------|-----|-------|---------|-------|
| 96 ns | ring_hold | 0.648 | >0.999 | <1e-4 | Discrimination-limited |
| 240 ns | ring_hold | 0.960 | 0.994 | <1e-4 | Sweet spot for most applications |
| 496 ns | ring_hold | 0.998 | 0.986 | boundary | Approaching T1 limit |

---

## Unified Limitations

### Shared Across All Studies
1. **No native stochastic continuous-measurement trajectories in cqed_sim**: QND metrics are replay-based, not trajectory-conditioned.
2. **No native readout-specific optimization objective in cqed_sim**: all optimization is study-local glue wrapping `simulate_sequence`.
3. **Detector efficiency η is representative (0.35), not device-calibrated**.

### Progressive Limitations
- Study #1: Linear model only; no multilevel effects, no decoherence.
- Study #2: QND defect is purely T1-limited; cannot distinguish pulse-shape effects on QND.
- Study #3: Phenomenological mixing layer; not validated against microscopic model or experiment.
- Study #4: Disturbance model uncalibrated to fridge data; stochastic SME replay unstable at high power.

### Physics Not Yet Captured
- True beyond-dispersive transmon dynamics at high readout power.
- Full 3-mode (storage + readout + transmon) readout simulation in the optimization loop.
- Amplifier chain noise (only readout chain response is currently modeled).
- AWG quantization effects on pulse fidelity.

---

## Knowledge Propagation

Results from earlier studies that **propagated** to later ones:
- Linear optimal drive frequency (midpoint) → used as default in #2–#4.
- Procedural pulse-family definitions → carried from #2 to #3 with hardware extensions.
- Nulling-tail construction → carried from #2 to #3; found to degrade under bandwidth limits.
- Effective mixing + higher-ladder framework → shared between #3 and #4.

Results that were **superseded**:
- #1's "GRAPE marginally improves over square" is model-specific to the linear regime; #3 shows the real competition is about QND cost, not SNR.
- #2's "QND defect is family-independent" was an artifact of the T1-only model; #3 showed family-dependent effects once mixing is included.

---

## Suggested Upstreaming (Deduplicated)

All four studies identified overlapping cqed_sim gaps. The deduplicated list:

| Priority | Extension | Studies Requesting |
|----------|-----------|-------------------|
| P1 | Readout-specific optimization objective API | #1, #2, #3 |
| P1 | Time-varying readout envelope simulation in `ReadoutResonator` | #1, #2 |
| P2 | Reusable hardware-profile presets for readout | #3 |
| P2 | Stochastic continuous-readout replay stabilization | #4 |
| P2 | Calibrated strong-readout disturbance presets or fitting utilities | #4 |
| P3 | Built-in readout benchmarking utility (detector-limited + QND-constrained frontiers) | #3 |
| P3 | Export reusable readout comparison CSV | #2 |
