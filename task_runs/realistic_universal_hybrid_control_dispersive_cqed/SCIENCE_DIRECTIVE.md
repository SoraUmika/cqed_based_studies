# Science Directive

## Objective
Determine whether the ideal universal primitive gate set for a dispersive transmon-cavity system survives the realistic Hamiltonian once `chi`, `chi'`, Kerr, and finite-duration phase accumulation are enforced, while explicitly avoiding the trap of making unconstrained GRAPE the conceptual answer.

## Ordered Actions
1. Audit the local `cqed_sim` API surface and the inherited repository studies relevant to displacement, spectator qubit rotations, selective pulses, strict SQR, relaxed CPSQR, arbitrary conditional control, and sequence-level replay.
2. Derive the analytic timing and phase-budget conflict between unconditional and selective control.
3. Build one machine-readable synthesis dataset plus publication figures from the inherited validated evidence.
4. Write a report that cleanly separates:
   - the ideal primitive picture,
   - the realistic Hamiltonian-level obstacles,
   - the constructive pulse strategies that survive,
   - the remaining universality gap.
5. Produce the reproducibility notebook and final handoff files.

## Success Criteria
- A precise top-level verdict on whether the strict ideal primitive library survives.
- At least one machine-readable summary artifact and one figure for the new synthesis.
- A compiled report and a reproducibility notebook.
