# Task Checklist

## Status Summary
- Study: `studies/storage_active_cooling_gf_sideband`
- Run: `task_runs/storage_active_cooling_gf_sideband`
- Problem class: `DES`, `ANA`, `OPT`
- Current state: complete and ready for review

## Initialize And Plan
- [x] Read AGENTS instructions and repo conventions
- [x] Confirm that `cqed_sim` already supports an effective multilevel sideband drive and readout decay
- [x] Create the study folder, state file, README, and improvement log
- [x] Write the study-scoped science directive
- [x] Confirm the exact device-parameter source and record it in the analysis artifacts

## Implement
- [x] Implement shared study helpers and device loading
- [x] Extract the dressed spectrum and label the relevant states up to `n=4`
- [x] Generate storage-sideband and readout-dump frequency tables
- [x] Run pulse-family amplitude-duration scans for the storage conversion step
- [x] Run readout-dump scans and assemble a full cooling primitive
- [x] Export machine-readable artifacts and report-quality figures

## Validate
- [x] Demonstrate analytic-versus-numeric agreement for the key resonance conditions
- [x] Demonstrate truncation and timestep convergence
- [x] Quantify leakage, residual transmon excitation, and spectral crowding
- [x] Compare the model-predicted `n=1` storage sideband with the locally stored nominal calibration value

## Report And Reproducibility
- [x] Update the README validation section and improvement log
- [x] Write `report/report.tex`
- [x] Compile `report/report.pdf`
- [x] Create `scripts/reproducibility_notebook.ipynb`
- [x] Update `EXECUTION_SUMMARY.md` and `REVIEW_REQUEST.md`
