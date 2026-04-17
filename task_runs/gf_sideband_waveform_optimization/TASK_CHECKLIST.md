# Task Checklist

## Status Summary
- Study: `studies/gf_sideband_waveform_optimization`
- Run: `task_runs/gf_sideband_waveform_optimization`
- Problem class: `OPT`, `ANA`, `DES`
- Current state: extended with simultaneous two-tone transfer analysis and ready for review

## Initialize And Plan
- [x] Read AGENTS instructions and repo conventions
- [x] Confirm that `cqed_sim` supports the relevant storage and readout sideband controls
- [x] Create the study folder, state file, README, and improvement log
- [x] Write the study-scoped science directive
- [x] Record the exact local device-parameter source in the exported artifacts

## Implement
- [x] Implement shared study helpers and waveform-family builders
- [x] Generate dressed frequency tables for storage and readout sidebands
- [x] Run waveform sweeps for both modes and all requested families
- [x] Export machine-readable comparison tables and representative traces
- [x] Export an explicit device-parameter provenance table
- [x] Add transmon-decoherence noise scenarios and rerun the open-system follow-up
- [x] Add a gate-oriented ranking artifact
- [x] Implement the simultaneous two-tone storage-to-readout transfer extension
- [x] Integrate the two-tone extension into the main study runner and export the new artifacts

## Validate
- [x] Demonstrate analytic-versus-numeric agreement for the resonance conditions
- [x] Demonstrate truncation and timestep convergence
- [x] Quantify selective versus unselective operating windows
- [x] Quantify robustness to amplitude and detuning error
- [x] Check whether the finalist rankings change under the available open-system channels
- [x] Check whether the shortlist winners remain threshold-valid under the matched transmon-reference reranking
- [x] Validate the single-photon two-tone reduced-ladder model against the full simulator
- [x] Check timestep and truncation convergence for the selected two-tone cases

## Report And Reproducibility
- [x] Update the README validation section and improvement log
- [x] Write `report/report.md`
- [x] Write `report/report.tex`
- [x] Compile `report/report.pdf`
- [x] Create `scripts/reproducibility_notebook.ipynb`
- [x] Update `EXECUTION_SUMMARY.md` and `REVIEW_REQUEST.md`
- [x] Extend the README, report, notebook, and handoff files with the simultaneous two-tone conclusions
