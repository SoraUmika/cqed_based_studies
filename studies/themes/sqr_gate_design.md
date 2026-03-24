# Theme: SQR Gate Design & Optimization

> Unified synthesis of 3 component studies covering waveform design, multi-branch operation, and open-system noise for selective qubit rotation.

---

## Component Studies

| Order | Study | Scope | Key Contribution |
|-------|-------|-------|-----------------|
| 1 | `sqr_pulse_waveform_design` | Closed-system waveform families + extended targets | Cosine-squared envelope optimal; cphase-SQR is natural primitive; GRAPE upper bound ~1; phase-compilation follow-up |
| 2 | `simultaneous_multitone_sqr_design` | Multi-branch simultaneous logical SQR | Block fidelity high but strict fidelity fails; cavity phase correction insufficient; GRAPE confirms waveform limitation |
| 3 | `sqr_open_system_deep_dive` | 3-mode open-system noise analysis | Multilevel relaxation, thermal photons, Purcell decay; GRAPE vs parametric under noise |

---

## Unified Problem Statement

Selective qubit rotation (SQR) is a foundational control primitive in dispersive cQED: it rotates the transmon qubit conditioned on the cavity Fock number n. A naive single-tone selective pulse achieves only conditional-phase SQR (desired rotation on target branch + unavoidable AC Stark phase on spectators), not true SQR (identity on spectators). This theme answers:

1. **Which waveform family best implements SQR?** (closed-system design)
2. **Can simultaneous multi-branch SQR work as a single common waveform?** (multi-branch viability)
3. **Which noise channels actually limit SQR in a real device?** (open-system constraints)

---

## Unified Methods

### System Model

All three studies share the same dispersive transmon-cavity system:
- ω_q/2π = 6.150 GHz, ω_c/2π = 5.241 GHz
- χ/2π = −2.84 MHz, χ′/2π = −21 kHz, K/2π = −28 kHz
- α/2π = −255 MHz
- Fock subspace: n = 0, 1, 2, 3 (N = 4 logical levels)
- Hilbert space: n_cav = 6, n_tr = 3 → dim = 18

Study #3 extends to 3 modes:
- ω_r/2π = 8.597 GHz, κ_r/2π = 2.4 MHz
- n_storage = 8, n_readout = 6, n_tr = 3

### Target Types

| Target | Definition | Natural primitive? |
|--------|-----------|-------------------|
| True SQR | R(θ,φ) on branch n₀, identity on all others | No — spectator phase unavoidable |
| Conditional-phase SQR (cphase-SQR) | R(θ,φ) on branch n₀, per-branch Z rotations allowed on spectators | **Yes — natural primitive** |
| Multi-branch simultaneous SQR | Same R(θ,φ) on all target branches simultaneously | Only at small angles or with GRAPE |

### Waveform Families Compared

| Family | #Params | Studies | Role |
|--------|---------|---------|------|
| Single-tone Gaussian | 3 | 1, 2 | Baseline |
| Cosine-squared (Hann) | 3 | 1 | **Best single-tone selectivity** |
| Phase-modulated Gaussian | 4 | 1 | Marginal improvement over Gaussian |
| One-segment multitone | 3–8 per tone | 1, 2 | Standard SQR builder |
| Two-segment echoed multitone | 2× multitone | 1 | Echo cancels some spectator phase |
| Optimized multitone (8p, 12p, 16p) | 8–16 | 1 | ≤1% improvement for χT/2π ≥ 3 |
| DE+NM optimized 4-tone (12p) | 12 | 1, 2 | F_block ≥ 0.987 but F_strict ≤ 0.27 |
| GRAPE (piecewise constant) | 30–120 | 1, 2, 3 | Upper bound: F ≥ 0.9999 |

### Key Dimensionless Parameter

The gate performance is primarily controlled by **χT/(2π)**, the number of dispersive cycles during the gate. All studies parameterize their scans in this quantity.

| χT/(2π) | Regime |
|---------|--------|
| < 1 | Severely bandwidth-limited; only GRAPE achieves useful fidelity |
| 1–2 | Transitional; optimized multitone helps significantly |
| **2–3** | **Sweet spot: best fidelity-duration tradeoff considering decoherence** |
| 3–5 | High selectivity; simple envelopes nearly saturate; diminishing returns |
| > 5 | Decoherence-dominated; fidelity degradation from T1 |

---

## Unified Key Results

### Result 1: Conditional-Phase SQR Is the Natural Primitive

Across all waveform families and gate durations, cphase-SQR fidelity substantially exceeds true-SQR fidelity. The spectator AC Stark phase is an inherent property of selective driving, not a pulse-design failure. Practical gate compilation should use cphase-SQR and correct the spectator phases in software.

**Implication:** Gate compilers targeting SQR should natively support per-branch Z corrections.

### Result 2: Cosine-Squared Envelope Best for Single-Branch SQR

Among single-tone envelopes, cosine-squared (Hann) provides the best off-resonant suppression and thus the highest cphase-SQR fidelity at intermediate χT/(2π). It reaches F_block ≥ 0.999 at χT/(2π) ≥ 3 with negligible |f⟩ leakage (< 10⁻⁴).

### Result 3: Simultaneous Multi-Branch SQR Fails at π Rotations

Study #2 showed that a common 4-tone waveform targeting R_X(π) on all branches n = 0..3 simultaneously achieves excellent per-branch rotation fidelity (F_block ≥ 0.987 via DE+NM optimization) but fails on strict logical fidelity (F_strict ≤ 0.27). The inter-branch phase coherence is fundamentally broken for the Gaussian parametric ansatz. GRAPE achieves F ≥ 0.9999, confirming the limitation is the ansatz, not the system.

**Phase-compilation analysis** (study #1 follow-up): cavity-only block-phase correction fails dramatically for the short-gate structured case (compiled fidelity drops from 0.138 to 0.028 at χT/2π = 0.5), confirming this is not a correctable phase problem.

### Result 4: Noise Budget Breakdown (Open System)

Study #3 quantified each noise channel's impact at the recommended operating point:

| Noise Channel | Fidelity Impact | Dominant? |
|--------------|----------------|-----------|
| Multilevel transmon relaxation (T1_ge, T1_fe) | −0.5–2% | **Yes**, at long gates |
| Storage thermal occupation (n_th ≤ 0.05) | −0.1–0.5% per 0.01 n_th | Moderate |
| Purcell decay from readout | Negligible for Δ_qr ≈ 2.4 GHz | No |
| Readout-storage crosstalk (3-mode) | Measurable at χ_sr ~ 50 kHz | Yes, if χ_sr ≠ 0 |

**Key finding:** GRAPE replay under realistic noise does NOT meaningfully outperform calibrated parametric baselines once decoherence is included. The closed-system GRAPE advantage (~1–5%) is partially eroded by decoherence, making the simpler parametric pulse the practical choice.

### Result 5: Optimal Operating Window

For representative targets on the nominal device:
- **Gate duration:** χT/(2π) ≈ 2–3 (i.e., T ≈ 700 ns – 1.1 μs at χ/2π = 2.84 MHz)
- **Envelope:** cosine-squared for single-branch; GRAPE only if F_block > 0.999 insufficient
- **Achievable fidelity:** F_block ≥ 0.99 (closed system), F_block ≈ 0.98 (with T1 = 30 μs)
- **Leakage:** < 10⁻⁴ to |f⟩

---

## Unified Limitations

### Closed-System Studies (#1, #2)
1. No decoherence in the optimization loop; performance may be overoptimistic at long gates.
2. n_tr = 2 for targeted-subspace validation; |f⟩ leakage estimated only via spot-check replay.
3. Phase-compilation analysis restricted to cavity-only corrections; full qubit+cavity gauge search not attempted.

### Open-System Study (#3)
1. Three-mode validation uses conservative baseline χ_sr = 0; real devices may have nonzero storage-readout cross-Kerr.
2. No native open-system SQR optimizer; closed-system pulses are replayed through noise rather than optimized within it.
3. GRAPE replay required careful numerical convergence (substrate-per-slice stabilization).

### Common Limitations
1. All studies use the same device parameter set; generality to other χ/κ regimes is not tested.
2. Simultaneous multi-branch SQR remains unsolved for large angles without GRAPE.
3. No experimental validation of any simulated pulse.

---

## Cross-Study Dependencies

```
sqr_pulse_waveform_design ──────────────┐
  │                                      │
  │  (baseline pulses, calibration)      │  (cavity phase analysis)
  ▼                                      │
sqr_open_system_deep_dive               │
  │                                      │
  │  (noise channels, GRAPE replay)      │
  ▼                                      ▼
[Merged SQR conclusions]    simultaneous_multitone_sqr_design
                                │
                                │  (multi-branch logical validation)
                                ▼
                          [Merged SQR conclusions]
```

---

## Connections to Other Themes

- **Readout theme** (#3, #4): The SQR gate must be fast enough to avoid decoherence from readout-chain Purcell effects quantified in the readout studies.
- **Hybrid control theme** (#8): SQR is a primitive in gate libraries A and B; the SQR gate quality directly limits the quality of hybrid qubit-cavity synthesis.
- **Adaptive control** (#10): GRAPE pulses designed in the SQR theme are the exact controls being re-optimized under chi mismatch in the gray-box study.

---

## Suggested Upstreaming (Deduplicated)

| Priority | Extension | Studies Requesting |
|----------|-----------|-------------------|
| P1 | Qutrit-compatible `targeted_subspace_multitone` diagnostic | #2 |
| P2 | Logical phase-compilation helpers (`global qubit-Z + cavity block-phase`) | #1 |
| P2 | Open-system SQR optimization objective (not just replay) | #3 |
| P3 | `targeted_subspace_multitone` → compiled logical-step operator export | #2 |
