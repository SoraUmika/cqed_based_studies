# Progress Log

## 2026-04-02T00:00:00Z - Study initialized and planned
- Objective: quantify how hot microwave components map onto qubit and cavity observables, coherence loss, multimode heating, and transient response in a dispersive cQED platform.
- Verified that `cqed_sim` supports dispersive single-mode and multimode models plus thermal Lindblad baths through `NoiseSpec`.
- Identified one small workflow gap: no dedicated public steady-state helper, so the study uses `cqed_sim` Hamiltonians and noise operators with a documented QuTiP steady-state bridge.
- Next: implement the shared study library and execute thrust A-D sweeps.

## 2026-04-02T06:45:00Z - Core simulations executed
- Implemented `common.py`, `run_study.py`, and `validate_results.py` for the four research thrusts.
- Completed steady-state thermometer, dephasing, multimode, and transient sweeps and saved JSON artifacts plus four figure sets.
- Realized runtimes were about 11 s for thermometry, 20 s for dephasing, 49 s for the multimode grid, and 3 s for the transient step response.

## 2026-04-02T07:10:00Z - Validation passed
- Ran seven validation checks spanning zero-temperature behavior, Bose scaling, truncation convergence, analytic low-occupation dephasing scaling, and multimode weak-coupling behavior.
- All validation checks passed and were saved to `studies/microwave_component_thermalization_cqed/artifacts/validation_summary.json`.

## 2026-04-02T07:30:00Z - Reports compiled
- Wrote the markdown memo and LaTeX report for the study.
- Compiled `report/report.pdf` successfully after the standard `pdflatex` and `bibtex` passes.
- No blocking LaTeX layout issues remained; only a benign float-placement warning was emitted.

## 2026-04-02T08:00:57.0332120Z - Reproducibility package completed
- Installed `joblib` with `python -m pip install --user joblib` and used threaded parallelism for the multimode map so the notebook stays stable on Windows/Jupyter.
- Added a notebook builder script and generated `scripts/reproducibility_notebook.ipynb`.
- Executed the reproducibility notebook end to end with `nbconvert`; full notebook runtime was about 112 s and the study remained reproducible.
- Wrote the execution summary and review-request handoff files. The run is now ready for independent review.
