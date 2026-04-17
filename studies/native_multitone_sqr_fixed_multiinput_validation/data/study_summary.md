# Native Multitone SQR Fixed Multi-Input Validation

## Residual-Z study baseline multitone
- Goal: Compare richer waveform families for residual-Z suppression while still attempting the target conditional rotations.
- Target: Random target-D block-diagonal conditional qubit rotations on 4 addressed Fock levels.
- Addressed levels: `4`
- Duration: `1056.3 ns`
- Fixed restricted process fidelity: `0.008405`
- Saved study restricted process fidelity: `0.008405`
- Shift after fixed reevaluation: `+0.000000`
- Single-input mean / min: `0.574001` / `0.131306`
- Pair-input mean / min: `0.486894` / `0.003145`
- Quartet mean / min: `0.443112` / `0.003145`

## Arbitrary-rotation study direct multitone
- Goal: Fit arbitrary block-diagonal conditional qubit rotations with a native multitone SQR waveform.
- Target: Structured family-C arbitrary conditional SU(2) target on 4 addressed Fock levels.
- Addressed levels: `4`
- Duration: `1056.3 ns`
- Fixed restricted process fidelity: `0.600671`
- Saved study restricted process fidelity: `0.600671`
- Shift after fixed reevaluation: `-0.000000`
- Single-input mean / min: `0.649514` / `0.345480`
- Pair-input mean / min: `0.670009` / `0.345480`
- Quartet mean / min: `0.734353` / `0.345479`

## Ideal-SQR study direct multitone
- Goal: Test whether direct or echoed multitone constructions can realize an ideal x-axis SQR profile.
- Target: Smooth x-axis ideal-SQR target on 3 addressed Fock levels.
- Addressed levels: `3`
- Duration: `1056.3 ns`
- Fixed restricted process fidelity: `0.824191`
- Saved study restricted process fidelity: `0.494918`
- Shift after fixed reevaluation: `+0.329273`
- Single-input mean / min: `0.832981` / `0.719799`
- Pair-input mean / min: `0.908966` / `0.719799`
- Quartet mean / min: `0.874682` / `0.719799`

## Corrected-metric study unitary-optimized
- Goal: Re-optimize a corrected direct multitone SQR under the fixed phase convention using a reduced effective-unitary metric.
- Target: Smooth corrected-SQR profile with explicit per-level (theta_n, phi_n) on 4 addressed Fock levels.
- Addressed levels: `4`
- Duration: `1056.3 ns`
- Fixed restricted process fidelity: `0.998648`
- Single-input mean / min: `0.999992` / `0.999978`
- Pair-input mean / min: `0.999990` / `0.999978`
- Quartet mean / min: `0.999989` / `0.999961`
