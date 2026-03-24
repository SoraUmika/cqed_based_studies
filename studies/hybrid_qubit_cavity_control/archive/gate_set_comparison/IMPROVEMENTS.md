# Improvement Log: Hybrid Universal Control Gate-Set Comparison in a 2x2 Qubit+Cavity Logical Subspace

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **Native synthesis primitives still stop at the ideal-operator layer**: `ConditionalDisplacement`, `JaynesCummingsExchange`, and `BlueSidebandExchange` now exist in `cqed_sim.unitary_synthesis`, but they are not yet connected to a waveform bridge or calibrated pulse-export workflow. This matters because the present native-gate comparison is physically motivated yet still cleaner than a full lab-transfer path. Next step: add pulse builders and waveform-bridge coverage, then replay the best native sequences through `SequenceCompiler` and `simulate_sequence`.
- **First-pass encoding remains restricted to Fock `{|0>, |1>}`**: conclusions about ECD- and SWAP-like control are likely encoding-sensitive. The present results are reliable for the 2x2 Fock block, but they may understate the advantages of displacement-native or exchange-native control in cat, binomial, or larger Fock encodings.
- **Direct `|f>` leakage is still omitted from the synthesis loop**: the decomposition benchmark uses `n_tr = 2` for tractability. This is probably fine for the present ranking, but it means the leakage numbers for sideband/exchange libraries are optimistic with respect to real multilevel transmon dynamics.

## Recommended Improvements (P2)
- **Upgrade winning libraries to multilevel noisy replay**: re-run Gate Set A, Gate Set D, and Gate Set F with `n_tr = 3` and realistic noise. This will quantify whether the native exchange library is hurt disproportionately by direct `|f>` leakage compared with the chi-wait entangler.
- **Benchmark a cat-like 2D cavity encoding**: repeat the same target suite on an even/odd cat logical pair once a clean study path is chosen. This is the most important follow-up for judging whether the poor ECD result is a property of the primitive or of the Fock encoding.
- **Add a local cavity-phase primitive to complement SQR**: the single-SQR entangler still reaches block-gauge fidelity `1.0` but strict fidelity only `0.7071` against `CX(c->q)`. A compact cavity block-phase correction could convert that into a competitive exact selective library.
- **Add warm-started and uncertainty-aware GRAPE**: the current GRAPE reference is useful but still sensitive to seed choice and slice count. A multistart warm-start workflow plus parameter-ensemble robustness would make the GRAPE-vs-decomposition gap easier to interpret.
- **Promote native pulse validation beyond the ideal operator model**: once waveform export exists, compare ideal `ConditionalDisplacement`, `JaynesCummingsExchange`, and `BlueSidebandExchange` against pulse-level implementations under the same duration prior to see how much native-gate performance survives realistic compilation.

## Nice-to-Haves (P3)
- **Extend targets beyond `CZ` and `CX`**: add random logical `SU(4)` samples or a small hybrid Clifford set.
- **Compress GRAPE solutions into short interpretable circuits**: test whether the best direct pulses suggest a reusable local cavity primitive that is faster than `D-SNAP-D`.
- **Push the native benchmark toward multimode memory tasks**: the SWAP-/sideband-style route is likely more compelling for transfer, encoding, and multimode state engineering than for exact 2x2 logical `CX`.

## Open Questions
- Does the native exchange family become the practical winner once the cavity encoding is cat-like or binomial instead of strict Fock `{|0>, |1>}`?
- How much of Gate Set F's deficit is true controllability mismatch for logical `CX`, and how much is just missing local phase/gauge correction around an otherwise good entangling primitive?
- Can an ECD- or CNOD-inspired local cavity primitive replace `D-SNAP-D` while keeping the exact native chi-wait entangler?
- Is the current native winner (`chi`-wait plus fast qubit rotations) still dominant for less diagonal entanglers, or does a more exchange-like library catch up once the target is not controlled-phase-derived?

## What Was Tried and Did Not Work
- **ECD-like `R-CD-R-CD-R` ansatz for the local cavity Hadamard**: after promoting conditional displacement to a first-class synthesis primitive, the best strict fidelity still stayed near `0.500` with leakage `0.454` at `440 ns`. The problem appears to be encoding mismatch rather than missing API support.
- **ECD-like `R-CD-R-CD-R` entangler ansatz**: remained near strict fidelity `0.500` and block fidelity `0.508` with negligible logical leakage. The primitive combination was unitary on the block but did not match the target entangler cleanly.
- **Minimal native library local control `D-W-D`**: reached only strict fidelity `0.554` with leakage `0.475`, and its worst perturbed replay dropped to `0.062`. Native dispersive waiting alone is an excellent entangler but a poor standalone local cavity-control resource.
- **Selective hybrid local library with two SQRs**: improved over the ECD-like and minimal-native local probes, but stalled at strict fidelity `0.861` and required `2440 ns`, so the added selectivity did not compensate for its time cost in the 2x2 benchmark.
- **Single-SQR exact entangler attempt**: produced strict fidelity `1/sqrt{2}` while achieving block-gauge fidelity `1.0`. This was not a failure of entangling power, but a mismatch between the native selective phase structure and the exact logical gauge of the chosen `CX(c->q)` target.
- **Native SWAP-/sideband-style local ansatz (`BS-JC-BS` with interleaved qubit rotations)**: this became a serious benchmark rather than a dead end. It reached strict fidelity `0.8784` with leakage `0.0788` in `440 ns`, outperforming selective Gate Set B on both fidelity and time, but still not matching the SNAP local gate.
- **Native SWAP-/sideband-style entangler ansatz**: reached strict fidelity `0.7510` with leakage `0.1834` in `320 ns`. It is better than the selective SQR entangler on strict fidelity and speed, but it still falls well short of the exact native chi-wait entangler.

## Compute & Resource Notes
- **Expanded study run**: `scripts/run_study.py` completed in approximately `73.3 s` for 10 decomposition cases, 6 GRAPE reference cases, robustness replays, and 5 summary figures.
- **Validation run**: `scripts/validate_results.py` completed in approximately `237.7 s` (`4.0 min`) for truncation reruns, deeper optimizer spot checks, GRAPE time-grid checks, and native-vs-selective sanity assertions.
- **Dominant cost**: GRAPE and validation remain the bottleneck; the native-gate decomposition cases are comparatively cheap.

## Resolved
- **No first-class ECD primitive in `unitary_synthesis`**: resolved by adding `ConditionalDisplacement` to the local `cqed_sim` copy, wiring it into the synthesis API, fast-evaluation backend, API reference, documentation, and dedicated unit tests.
- **No first-class SWAP-/sideband-style native primitives in `unitary_synthesis`**: resolved by adding `JaynesCummingsExchange` and `BlueSidebandExchange` to the local `cqed_sim` copy, together with tests and API/docs updates.
- **Native dispersive waiting was only an ad hoc study trick**: resolved at the study level by elevating chi-wait control to its own benchmark family (Gate Set D) and by making the native-vs-selective comparison explicit in the report.
