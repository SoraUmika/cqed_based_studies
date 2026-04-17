# Execution Summary

Date: 2026-04-10
Study: `studies/gf_sideband_waveform_optimization`
Run: `task_runs/gf_sideband_waveform_optimization`

## Main Findings
- Closed-system selective winner:
  - storage sideband: `gaussian (sigma_fraction = 0.24)` with conservative duration `300 ns`
  - readout sideband: `gaussian (sigma_fraction = 0.24)` with conservative duration `220 ns`
- Closed-system fast/unselective winner:
  - storage sideband: `square` at `12 ns`
  - readout sideband: `square` at `12 ns`
- Mode-noise-only replay:
  - storage recommendations remain essentially unchanged
  - selective readout control collapses because the readout linewidth dominates the `220 ns` pulse
- Transmon-inclusive replay using the matched local reference (`T1 = 9.813 us`, `T2,Ramsey = 6.325 us`):
  - the storage selective Gaussian drops to mean target transfer `0.96075`
  - the readout fast square drops to mean target transfer `0.92300`
  - only the storage fast-square pulse remains threshold-valid under the original metrics
- Simultaneous two-tone extension for `|g,0,1> -> |g,1,0>`:
  - the resonant bright-state protocol with equal `12 MHz` leg couplings reaches `0.9999977` closed-system target transfer in `29.5 ns`
  - the best detuned Raman-like case lowers the peak intermediate-state occupation to `0.121` but slows to `348 ns`
  - under the noisy replay, the resonant single-photon case still peaks near `0.953`, whereas the detuned case falls to about `0.556`
  - for `n > 1`, the same constant simultaneous drive opens a longer ladder and no longer delivers clean near-unit transfer

## Key Numbers
- Storage selective per-`n` durations: `220 ns`, `300 ns`, `300 ns`
- Readout selective per-`n` durations: `220 ns`, `220 ns`, `220 ns`
- Storage unselective per-`n` durations: `12 ns`, `12 ns`, `12 ns`
- Readout unselective per-`n` durations: `12 ns`, `12 ns`, `12 ns`
- Mode-only noisy target means:
  - storage selective Gaussian: `0.99225`
  - storage fast square: `0.99542`
  - readout selective Gaussian: `0.33126`
  - readout fast square: `0.92441`
- Transmon-reference noisy target means:
  - storage selective Gaussian: `0.96075`
  - storage fast square: `0.99391`
  - readout selective Gaussian: `0.32240`
  - readout fast square: `0.92300`
- Two-tone single-photon selected cases:
  - resonant closed system: peak target `0.9999977`, peak time `29.5 ns`, peak intermediate `0.49975`
  - resonant mode-only noise: peak target `0.95542` at `29.25 ns`
  - resonant transmon reference: peak target `0.95303` at `29.25 ns`
  - detuned closed system: peak target `0.99677`, peak time `348.0 ns`, peak intermediate `0.12108`
  - detuned mode-only noise: peak target `0.56063` at `303.25 ns`
  - detuned transmon reference: peak target `0.55599` at `303.25 ns`
- Gate-oriented winners:
  - storage: `gaussian`, mean projected SWAP fidelity `0.716`
  - readout: `blackman`, mean projected SWAP fidelity `0.765`
- Extended main runtime: `534.8 s`
- Two-tone extension runtime: `21.1 s`
- Validation runtime after extension: `3.6 s`

## Validation
- Analytic sanity: exact rotating-frame sideband frequencies agree with the dispersive formulas to numerical precision for both modes and `n = 1,2,3`.
- Timestep convergence: using `0.5 ns` instead of `0.25 ns` changes finalist target transfer by only `1.8e-4` to `1.4e-3`.
- Truncation convergence: increasing from `(n_tr, n_s, n_r) = (4,5,5)` to `(5,6,6)` changes finalist target transfer by at most `5.7e-5`.
- Two-tone reduced-model validation: the selected single-photon resonant case matches the reduced three-state model to within `4.1e-7` in peak target probability and exactly in peak time.
- Two-tone convergence: the selected single-photon cases change by at most `2.5e-3` in peak target probability at `1.0 ns` resolution and by less than `2.3e-6` under the larger `(5,6,6)` truncation.

## Important Caveats
- The waveform ranking is still performed inside the effective sideband-control abstraction rather than a microscopic pump Hamiltonian.
- The transmon-noise extension uses a matched local tomography workflow as a sensitivity anchor because the sideband-reset example itself does not export transmon `T1/T2`.
- The winning pulses are strong state-transfer primitives but not phase-clean coherent SWAP gates in the projected two-state metric.
- The simultaneous two-tone extension only studies constant square drives, so it does not yet test shaped counter-intuitive overlap, STIRAP-like timing, or direct open-system optimal control.
- The clean two-tone reduced-ladder mechanism is specific to the single-photon case; for `n > 1`, the drive couples into a longer chain and the protocol is no longer a clean one-step transfer.

## Key Outputs
- `studies/gf_sideband_waveform_optimization/report/report.md`
- `studies/gf_sideband_waveform_optimization/report/report.tex`
- `studies/gf_sideband_waveform_optimization/report/report.pdf`
- `studies/gf_sideband_waveform_optimization/report/report_updated.pdf`
- `studies/gf_sideband_waveform_optimization/scripts/reproducibility_notebook.ipynb`
- `studies/gf_sideband_waveform_optimization/data/device_parameter_table.csv`
- `studies/gf_sideband_waveform_optimization/data/noise_scenarios.csv`
- `studies/gf_sideband_waveform_optimization/data/gate_winner_table.csv`
- `studies/gf_sideband_waveform_optimization/data/open_system_summary.csv`
- `studies/gf_sideband_waveform_optimization/data/open_system_transmon_reference_winner_table.csv`
- `studies/gf_sideband_waveform_optimization/data/two_tone_selected_cases.csv`
- `studies/gf_sideband_waveform_optimization/data/two_tone_open_system.csv`
- `studies/gf_sideband_waveform_optimization/data/two_tone_validation.csv`
- `studies/gf_sideband_waveform_optimization/data/validation_convergence.csv`
- `studies/gf_sideband_waveform_optimization/figures/two_tone_scan_summary.png`
- `studies/gf_sideband_waveform_optimization/figures/two_tone_population_dynamics.png`
- `studies/gf_sideband_waveform_optimization/figures/two_tone_open_system_summary.png`
