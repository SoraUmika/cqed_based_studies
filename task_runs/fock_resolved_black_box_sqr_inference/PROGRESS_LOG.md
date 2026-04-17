# Progress Log

## 2026-03-31T05:46:41Z - Study initialized
- Objective: Validate whether qubit-state tomography combined with known dispersive wait evolution and calibrated cavity displacement can recover effective Fock-resolved qubit output states of a black-box SQR gate using forward simulation, constrained least squares, and MLE ground-truth validation.
- Study path: studies/fock_resolved_black_box_sqr_inference
- Run path:   task_runs/fock_resolved_black_box_sqr_inference

## 2026-03-31T05:58:00Z - API review and analytic preliminary completed
- Confirmed local `cqed_sim` installation at `C:\Users\jl82323\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation`.
- Reviewed `cqed_sim` support for dispersive simulation, displacement operators, conditioned-qubit state extraction, ideal SQR operators, and multitone waveform construction.
- Derived the core identifiability result before implementation:
  - displacement-only is exactly uninformative after tracing out the cavity,
  - wait-only can recover weighted transverse sector information,
  - combined `D(alpha) -> wait` can reveal cavity coherences but still cannot recover `p_n` or sector-wise `Z_n` separately.

## 2026-03-31T06:08:00Z - Pulse-level feasibility scratch test completed
- Verified that `cqed_sim.calibration.conditioned_multitone` can synthesize a near-ideal 4-sector SQR-like waveform on an `n_tr = 2`, `n_cav = 4` model in under one minute.
- Logged a resolved helper limitation for `n_tr = 3` full-mode validation and selected a local truncation workaround for leakage-oriented case studies.

## 2026-03-31T09:20:00Z - Full study execution completed
- Implemented the single-qubit baseline, the diagonal-model kernels, recoverable least-squares / MLE solvers, the exact joint-state coherence witness, and the pulse-level black-box case library.
- Ran the full-profile study and validation scripts. All automated validation checks passed.
- Generated machine-readable outputs, eight publication figures, the reproducibility notebook, and the compiled report PDF.
- Key numerical conclusions:
  - single-qubit baseline MLE reached mean fidelity `0.9986` at the high-shot end,
  - wait-only and combined protocols are both full-rank on the recoverable transverse subspace,
  - displacement-only remains rank deficient and non-identifiable,
  - the combined protocol cleanly exposes cavity coherence through a large residual,
  - pulse-level near-diagonal cases are reconstructed accurately, while coherent and leakage-prone cases fail in diagnostically useful ways.
