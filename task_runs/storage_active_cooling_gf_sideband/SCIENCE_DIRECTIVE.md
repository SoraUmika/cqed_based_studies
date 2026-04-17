# Science Directive

Date: 2026-04-07
Study: `studies/storage_active_cooling_gf_sideband`
Run: `task_runs/storage_active_cooling_gf_sideband`

## Problem Class
DES, ANA, OPT

## Central Question
Can the local `cqed_sim` device support a practically useful active-cooling primitive for the storage cavity based on the transmon `|f>` manifold, specifically `|g,n\rangle \leftrightarrow |f,n-1\rangle` conversion followed by a readout-assisted dissipative dump, for `n \le 4`?

## Scope
1. Use the actual device parameters exposed by the local editable `cqed_sim` sideband-reset workflow.
2. Derive the storage-conversion and readout-dump resonance conditions in both the lab frame and the chosen rotating frame.
3. Extract and tabulate actual frequencies, detunings, and matrix-element scalings for `n=1,2,3,4`.
4. Compare at least three pulse-envelope families on the effective storage `g-f` sideband channel.
5. Simulate a full cooling primitive with readout decay, including repeated cycles and realistic failure modes.
6. Produce experiment-facing calibration and validation guidance.

## Non-Negotiable Modeling Boundaries
- Use `cqed_sim` as the primary simulation engine.
- Treat the `SidebandDriveSpec` channel as an effective activated interaction, not as a microscopic derivation of the pump.
- Distinguish clearly between:
  - direct transmon carrier drives,
  - effective storage/readout sideband drives,
  - irreversible readout ringdown.
- Do not claim hardware transferability for absolute pump amplitude without calibration evidence.

## Required Deliverables
- Frequency table and viability table for `n=1..4`
- Pulse-scan heatmaps and representative population traces
- Cooling-cycle performance plots and repeated-cycle diagnostics
- Machine-readable artifacts for every headline result
- Technical report, compiled PDF, notebook, and handoff files
