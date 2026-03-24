# Theme: Hybrid Qubit-Cavity Universal Control

> Unified synthesis of 2 component studies covering gate-library benchmarking and specific hybrid unitary synthesis.

---

## Component Studies

| Order | Study | Scope | Key Contribution |
|-------|-------|-------|-----------------|
| 1 | `hybrid_universal_control_gate_set_comparison` | 6-library head-to-head comparison for 2×2 qubit+cavity control | Ranking of SNAP, SQR, ECD, native, GRAPE, sideband libraries; new cqed_sim primitives |
| 2 | `utarget_cqed_decomposition` | Deep synthesis of U_target = (I⊗H)·CNOT·CNOT | 6-phase variational synthesis; algebraic decomposition; robustness analysis |

---

## Unified Problem Statement

Hybrid qubit-cavity systems in cQED encode quantum information in both the transmon and the storage cavity. Universal control of the joint system requires entangling operations between the two subsystems. The central question: **which physical primitive library provides the best route to universal hybrid control at the device level?**

Study #1 (gate-set comparison) answers this broadly by benchmarking 6 candidate libraries. Study #2 (U_target synthesis) answers it in depth for a specific maximally-entangling target requiring three CNOT-equivalent operations.

---

## Shared Framework

### Logical Subspace

Both studies operate on the **4D logical subspace** {|g,0⟩, |g,1⟩, |e,0⟩, |e,1⟩} with Fock encoding |0⟩_L = |0⟩, |1⟩_L = |1⟩.

### Device Parameters

| Parameter | Value |
|-----------|-------|
| χ/2π | −2.84 MHz |
| χ′/2π | −21 kHz |
| K/2π | −28 kHz |
| α/2π | −255 MHz |
| n_tr | 2 (synthesis), 3 (validation replay) |
| n_cav | 6–10 (convergence tested) |

### Shared cqed_sim API

Both studies use: `DispersiveTransmonCavityModel`, `FrameSpec`, `UnitarySynthesizer`, `GateSequence`, `Subspace`, `TargetUnitary`, `LeakagePenalty`, `MultiObjective`, `subspace_unitary_fidelity`, `leakage_metrics`.

### Gate Primitive Libraries

| Library | Primitives | Study #1 Label | Study #2 Label |
|---------|-----------|----------------|----------------|
| SNAP-based | {R_q, D, SNAP} + native chi-wait | Gate Set A | Library C |
| SQR-based | {R_q, D, SQR, ConditionalPhaseSQR} | Gate Set B | Library B |
| ECD-like | {R_q, D, ConditionalDisplacement} | Gate Set C | — |
| Native chi-wait | {R_q, D, FreeEvolveCondPhase} | Gate Set D | Library A |
| GRAPE | Piecewise-constant waveform | Gate Set E | — |
| Sideband | {R_q, JC Exchange, Blue Sideband} | Gate Set F | — |

---

## Unified Key Results

### Result 1: Native Chi-Wait + Fast Qubit Rotations Is the Best Entangler

Both studies independently confirm: the **native dispersive chi-wait entangler** provides the best practical entangling operation for the 2×2 Fock logical subspace.

- Study #1: Gate Set D achieves F_strict = 0.9999999 with negligible leakage in 256 ns for the CZ entangler.
- Study #2: Library A (D + R_q + SQR + FreeEvolveCondPhase) at depth-11 achieves F = 0.9193 with leakage 0.101 for the more demanding U_target.

The chi-wait entangler is exact (no approximation), fast (one dispersive period), and requires no microwave drive — it is the passive dispersive evolution itself.

### Result 2: SNAP Library Is Best for Local Cavity Control

Study #1 conclusively shows Gate Set A (D-SNAP-D) provides the best local cavity Hadamard with F_strict = 0.9887, leakage = 0.0185, duration = 1260 ns. No other library matched this for local cavity operations in the Fock encoding.

### Result 3: SQR-Based Library Has Structural Phase Mismatch

Both studies found that SQR-based entanglers reach perfect block-gauge fidelity (1.0) but suffer strict-fidelity deficits due to a cavity-phase mismatch:
- Study #1: Single-SQR entangler has F_strict = 1/√2, F_block = 1.0.
- Study #2: Library B achieves the best physical fidelity (F = 0.9193), but the cavity-phase structure of SQR is non-trivial.

A local cavity block-phase correction could potentially resolve this, but study #1's cavity-phase compilation analysis showed this is not always achievable.

### Result 4: ECD and Sideband Libraries Underperform in Fock Encoding

- Gate Set C (ECD-like): F_strict ≈ 0.50 for local cavity and entangler targets.
- Gate Set F (sideband): F_strict = 0.8784, leakage = 0.079, duration = 440 ns for local cavity; fast but leaky.

Both may improve significantly under different cavity encodings (cat, binomial) where displacement-based and exchange-based operations are more natural.

### Result 5: GRAPE Provides a Strong but Not Exhaustive Upper Bound

- Study #1: Best GRAPE local cavity F_strict = 0.9618 at 320 ns; best GRAPE entangler F_strict = 0.9458 at 400 ns. Sensitive to seed and time-grid.
- Study #2: GRAPE was not the primary approach, but variational synthesis explored the same synthesis space.

### Result 6: U_target Synthesis Demonstrates Multi-CNOT Hybrid Operations Are Feasible

Study #2 showed that a 3-CNOT-equivalent hybrid operation U_target = (I⊗H)·CNOT_{c→q}·CNOT_{q→c} can be synthesized with F > 0.99 using ideal gates, with the physical implementation achieving F = 0.9193. The dominant error is leakage (0.101), not coherent mismatch.

Robustness analysis: ±5% χ error → ~0.58% RMS fidelity change; ±5% amplitude → ~0.58% RMS; ±5% phase → ~0.12% RMS. Duration sensitivity is essentially zero (optimizer compensates drift automatically).

---

## Unified Limitations

### Fundamental Constraints
1. **Fock {|0⟩, |1⟩} encoding is restrictive**: conclusions about gate libraries may change for cat, binomial, or larger Fock encodings.
2. **n_tr = 2 in synthesis loop**: |f⟩ leakage numbers are optimistic for sideband/exchange libraries.
3. **No pulse-level export for native primitives**: ConditionalDisplacement, JaynesCummingsExchange, and BlueSidebandExchange are ideal operators, not calibrated waveforms.
4. **No open-system optimization**: synthesis is closed-system; noise replay is validation-only.

### Study-Specific Limitations
- Study #1: GRAPE results sensitive to seed; multistart budget was limited.
- Study #2: Physical fidelity (F = 0.9193) limited by leakage at depth-11; deeper sequences or GRAPE on the full waveform might improve.

---

## Cross-Study Dependencies

```
hybrid_universal_control_gate_set_comparison
    │
    │  (provides library ranking, new cqed_sim primitives)
    │
    ├──── utarget_cqed_decomposition
    │      (applies ranked libraries to specific hard target)
    │
    └──── SQR Theme (sqr_pulse_waveform_design, etc.)
           (SQR quality feeds into Library B performance)
```

---

## Connections to Other Themes

- **SQR theme**: SQR is a primitive in Libraries A and B. SQR gate quality (cosine-squared envelope, χT/(2π) ≥ 2–3) directly determines synthesis fidelity.
- **Readout theme**: Readout-chain Purcell effects and measurement-induced dephasing affect the qubit during long synthesis sequences.
- **Adaptive control** (#10): Gray-box strategies would be needed for the SQR and SNAP primitives in these libraries if χ drifts.
- **RESEARCH_PLAN.md Study A3**: Pulse-level synthesis of U_target is the proposed next step, bridging ideal-gate results to physical waveforms.

---

## Suggested Upstreaming (Deduplicated)

| Priority | Extension | Studies Requesting |
|----------|-----------|-------------------|
| P1 | Waveform bridge / pulse export for native synthesis primitives | #1 (ConditionalDisplacement, JC, BlueSideband) |
| P2 | Cat/binomial encoding subspace definition for synthesis | #1 |
| P2 | Multistart warm-started GRAPE with parameter-ensemble robustness | #1 |
| P2 | Logical CNOT_{q→c} gate in the built-in target library | #2 |
| P3 | Random SU(4) and hybrid Clifford target sets for benchmarking | #1 |
| P3 | GRAPE solution → interpretable circuit decomposition | #1 |
