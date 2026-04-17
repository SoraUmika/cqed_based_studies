# Study Summary: Fock-Resolved SQR Inference v2

**Profile:** full  **Runtime:** 105.3s

## Part 1 — Single-Qubit Baseline
Best MLE mean fidelity: **0.9982**

## Part 2A — Identifiability
| Protocol | Transverse rank | Full rank | Cond. # |
|---|---|---|---|
| wait_only | 8 | 9 | 1.21 |
| displacement_only | 2 | 3 | ∞ |
| combined | 8 | 9 | 2.34 |

## Part 2B — Full-State Non-Uniqueness
Objective span across 8 restarts: **220.44**
Population σ by sector: ['0.038', '0.044', '0.057', '0.046']
Gauge family size: **2**

## Part 2C — Recoverable Subspace (MLE)
| Protocol | N=1000 RMSE | N=10000 RMSE |
|---|---|---|
| wait_only | 0.0126 | 0.0036 |
| displacement_only | 0.1739 | 0.1739 |
| combined | 0.0108 | 0.0029 |

## Part 2F — Coherence Sweep
Coherence witness gap (combined, f=1.0 vs f=0.0): **0.1387** vs **2.03e-16**

## Part 2E — Pulse-Level Cases
Best case: **cpsqr_pulse_mix_+** RMSE=0.0028
Worst case: **short_pulse_seed_superpos_g** RMSE=0.0643

## Part 3 — Comparison Answers
Q1 Wait-only recovers transverse: **True**
Q2 Displacement needed for p_n/Z_n: **False**
Q4 Coherence witness works: **True**
Q5 Useful up to coherence fraction: **0.25**
Q6 Recommended protocol: **wait_only**