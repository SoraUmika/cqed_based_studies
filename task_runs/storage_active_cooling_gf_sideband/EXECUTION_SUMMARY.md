# Execution Summary

Date: 2026-04-08
Study: `studies/storage_active_cooling_gf_sideband`
Run: `task_runs/storage_active_cooling_gf_sideband`

## Main Result
The local device supports a practically useful active-cooling primitive based on the user-requested ladder
`|g,0_r,n_s> -> |f,0_r,n_s-1> -> |g,1_r,n_s-1> -> |g,0_r,n_s-1>`
through `n_s <= 4`, provided the transitions are activated as effective storage and readout red sidebands rather than as direct transmon carrier drives.

The study now explicitly uses **Fast Sideband Control of a Weakly Coupled Multimode Bosonic Memory** (`arXiv:2503.10623v1`) as a technical guide for:
- pulse-family selection,
- photon-number-resolved calibration logic,
- bosonic enhancement checks,
- and strong-drive Stark/Floquet screening.

## Frequency Digest
- Step A frequencies (`|g,0_r,n_s> <-> |f,0_r,n_s-1>`, `GHz`):
  - `6.804115034`
  - `6.798434192`
  - `6.792753350`
  - `6.787072508`
- Step B frequencies (`|f,0_r,n_s-1> <-> |g,1_r,n_s-1>`, `GHz`):
  - `3.448825278`
  - `3.443144436`
  - `3.437463594`
  - `3.431782752`
- Adjacent line spacing for both ladders: `5.680842 MHz`
- Dressed-state overlaps with the intended labels remain `1.0` for every tracked state in the present dispersive model

## Pulse Recommendations
- Best Step A pulses:
  - `n=1`: `bump`, `12 MHz`, `20 ns`, transfer `0.996051`
  - `n=2`: `square`, `6 MHz`, `30 ns`, transfer `0.999421`
  - `n=3`: `bump`, `14 MHz`, `10 ns`, transfer `0.997767`
  - `n=4`: `bump`, `12 MHz`, `10 ns`, transfer `0.996048`
- Best Step B pulses:
  - all `n=1..4`: `cosine_squared`, `12 MHz`, `20 ns`, transfer about `0.99605`
- Strong-drive screening:
  - local detuning scans over `+-2 MHz` find the optimum at `0.0 MHz` for every recommended Step A and Step B pulse in the effective model
  - Floquet doublet splittings track the expected bosonic scaling for Step A and the near-constant readout scale for Step B

## Cooling Performance
- Single-cycle success probabilities:
  - `n=1`: `0.973509`
  - `n=2`: `0.972875`
  - `n=3`: `0.967468`
  - `n=4`: `0.962020`
- Dominant single-cycle leakage channel:
  - residual readout population in `|g,1_r,n_s-1>`
- Repeated ladder cooling:
  - initial `|g,0_r,4>` -> final mean storage occupation `0.075741`
  - coherent test state (`alpha=1.1`) -> final mean storage occupation `0.031285`
  - thermal-like test state -> final mean storage occupation `0.012508`

## Validation
- Analytic resonance conditions match the numerically extracted simulator frequencies exactly to the printed precision for `n=1..4`
- Step A sideband matrix elements follow the expected bosonic ladder exactly: `1`, `sqrt(2)`, `sqrt(3)`, `2`
- Convergence for the hardest `n=4` primitive:
  - nearby truncation changes move success probability by at most `4.63e-05`
  - halving the timestep changes success probability by `9.22e-05`
  - doubling the timestep to `0.5 ns` changes success probability by `-3.39e-04`
- Follow-up weak-drive spectroscopy scans:
  - Step A and Step B peaks occur at `0.0 MHz` detuning relative to the dressed model predictions for all `n=1..4`
  - direct `g-f` carrier control keeps storage photon number fixed and is therefore a calibration tool, not the cooling transition
  - the recommended first-pass Step A measurement is `P_f`, with a dump-assisted witness as the fallback when direct `f` readout is weak

## Important Caveats
- The native framework treats the storage and readout sidebands as effective control operators. The study therefore validates control-layer feasibility, spectroscopy, and cooling performance, but it does not derive the microscopic pump coupling.
- The Floquet check is useful but still qualitative for leakage: every recommended strong-drive case triggers a truncation-boundary warning.
- The old nominal `n=1` storage-sideband calibration stored in the local environment remains offset from the static dressed model by `581.406 kHz`, so hardware spectroscopy is still required before transfer to experiment.
- During the final verification/improvement pass, `report.pdf` could not be overwritten because the file was locked by an open viewer. The current compiled artifact for the corrected source is `report_verify.pdf`.

## Delivered Files
- Report:
  - `studies/storage_active_cooling_gf_sideband/report/report.tex`
  - `studies/storage_active_cooling_gf_sideband/report/report.pdf` (pre-verification locked build)
- Verification compile artifact:
  - `studies/storage_active_cooling_gf_sideband/report/report_verify.pdf` (current verification-improved build)
- Notebook:
  - `studies/storage_active_cooling_gf_sideband/scripts/reproducibility_notebook.ipynb`
- Main machine-readable summaries:
  - `studies/storage_active_cooling_gf_sideband/data/study_results.json`
  - `studies/storage_active_cooling_gf_sideband/data/frequency_table.csv`
  - `studies/storage_active_cooling_gf_sideband/data/recommendation_table.csv`
  - `studies/storage_active_cooling_gf_sideband/data/frequency_calibration_scan.csv`
  - `studies/storage_active_cooling_gf_sideband/data/floquet_summary.csv`
  - `studies/storage_active_cooling_gf_sideband/data/spectroscopy_followup_summary.csv`
  - `studies/storage_active_cooling_gf_sideband/data/spectroscopy_followup_step_a.csv`
  - `studies/storage_active_cooling_gf_sideband/data/spectroscopy_followup_step_b.csv`
- Experiment-facing follow-up note:
  - `studies/storage_active_cooling_gf_sideband/artifacts/spectroscopy_measurement_followup.md`
- Review:
  - `task_runs/storage_active_cooling_gf_sideband/REVIEW_DIRECTIVE.md`
