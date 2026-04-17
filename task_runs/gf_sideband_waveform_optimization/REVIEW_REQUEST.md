# Review Request

Date: 2026-04-10
Study: `studies/gf_sideband_waveform_optimization`
Run: `task_runs/gf_sideband_waveform_optimization`

## Review Focus
Please review the study as a control-and-physics note answering which waveform families are best for the effective
`|f,n-1> <-> |g,n>` sideband interaction in storage and readout modes, with special attention to the new parameter-provenance and transmon-decoherence extensions.
Please also review the new simultaneous two-tone extension that asks whether a continuous storage-sideband plus readout-sideband drive can move a storage photon into the readout mode efficiently under the local cQED parameters.

The main claims to audit are:
1. A wide Gaussian pulse is the best closed-system selective family for both mode types.
2. A square pulse is the best closed-system fastest-transfer family for both mode types.
3. The exact cQED parameter tuple used for the Hamiltonian and mode noise is now fully explicit in the saved artifacts.
4. Selective readout-sideband control is not practically useful under either the mode-only noisy replay or the transmon-inclusive replay.
5. Under the matched transmon-reference reranking, only the storage fast-square pulse remains threshold-valid under the original metrics.
6. The added gate-oriented ranking shows that the analytic family set still does not produce a phase-clean coherent SWAP across `n = 1,2,3`.
7. In the single-photon manifold, a resonant simultaneous two-tone drive realizes a bright-state transfer that reaches essentially unit closed-system target population in about `29.5 ns` and still peaks near `0.953` under the noisy replay.
8. The detuned Raman-like simultaneous two-tone variant lowers the intermediate-state occupation but is too slow to beat the measured readout decay, so it falls to about `0.556` under the noisy replay.
9. The clean reduced three-state simultaneous-transfer picture is specific to `n = 1`; the same constant continuous drive does not generalize cleanly to `n > 1` because it opens a longer conversion chain.

## Primary Artifacts
- `studies/gf_sideband_waveform_optimization/report/report.pdf`
- `studies/gf_sideband_waveform_optimization/report/report_updated.pdf` (`report.pdf` could not be overwritten because it is locked by another process)
- `studies/gf_sideband_waveform_optimization/report/report.tex`
- `studies/gf_sideband_waveform_optimization/report/report.md`
- `studies/gf_sideband_waveform_optimization/scripts/reproducibility_notebook.ipynb`
- `studies/gf_sideband_waveform_optimization/data/device_parameter_table.csv`
- `studies/gf_sideband_waveform_optimization/data/noise_scenarios.csv`
- `studies/gf_sideband_waveform_optimization/data/gate_winner_table.csv`
- `studies/gf_sideband_waveform_optimization/data/open_system_summary.csv`
- `studies/gf_sideband_waveform_optimization/data/open_system_transmon_reference_winner_table.csv`
- `studies/gf_sideband_waveform_optimization/data/validation_convergence.csv`
- `studies/gf_sideband_waveform_optimization/data/two_tone_frequency_table.csv`
- `studies/gf_sideband_waveform_optimization/data/two_tone_selected_cases.csv`
- `studies/gf_sideband_waveform_optimization/data/two_tone_open_system.csv`
- `studies/gf_sideband_waveform_optimization/data/two_tone_validation.csv`
- `studies/gf_sideband_waveform_optimization/figures/two_tone_scan_summary.png`
- `studies/gf_sideband_waveform_optimization/figures/two_tone_population_dynamics.png`
- `studies/gf_sideband_waveform_optimization/figures/two_tone_open_system_summary.png`
