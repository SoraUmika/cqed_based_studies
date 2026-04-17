# Progress Log

## 2026-04-14 - Initialization
- Created the study and task-run folder structure.
- Audited the local `cqed_sim` API reference and confirmed the relevant dispersive-model and multitone-calibration surfaces.

## 2026-04-14 - Evidence Gathering
- Read the completed repository studies on waveform-level gate realism, literature-informed selective primitives, strict SQR versus relaxed CPSQR, arbitrary conditional control, and runtime hybrid-unitary replay.
- Identified the central inherited quantitative results needed for a universality verdict.

## 2026-04-14 - Synthesis Build
- Implemented `build_synthesis_dataset.py`.
- Generated `data/synthesis_summary.json`, the primitive-verdict artifact, the analytic phase-budget artifact, and the timing-hierarchy / phase-budget figures.

## 2026-04-14 - Reporting
- Implemented the reproducibility-notebook builder and generated `scripts/reproducibility_notebook.ipynb`.
- Wrote `report/report.tex` and `report/references.bib`.
- Compiled the PDF successfully with `pdflatex -> bibtex -> pdflatex -> pdflatex`.

## 2026-04-14 - Closeout
- Wrote README, IMPROVEMENTS, and `study_state.json`.
- Wrote the task-run handoff files.
- Executed the reproducibility notebook successfully through `python -m jupyter nbconvert --execute`.
- Completed the reviewer pass with an APPROVE decision and final polish record.
- Final synthesis verdict: the strict ideal primitive gate set does not survive literally; a weaker phase-aware constructive library does survive, but a fully pulse-backed non-GRAPE universal stack is still open.
