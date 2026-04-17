---
description: "LaTeX report conventions for cQED study reports. Enforces prose standards, mandatory sections, and formatting rules from AGENTS.md."
applyTo: "**/report/*.tex"
---

# LaTeX Report Conventions

## Template
- Use `revtex4-2` two-column format: `\documentclass[aps,pra,twocolumn,reprint,amsmath,amssymb]{revtex4-2}`
- Required packages: `graphicx`, `booktabs`, `siunitx`, `hyperref`, `float`, `xcolor`

## Mandatory Sections (all required)
1. Abstract (150-300 words, self-contained)
2. Introduction (context, prior work, objectives, roadmap)
3. System and Methods (Hamiltonian, parameters table, analytic preliminary, computational approach)
4. Results (one subsection per study goal)
5. Validation (sanity checks, convergence, literature comparison)
6. Discussion (physical interpretation, limitations)
7. Conclusion
8. Limitations and Future Work (known limitations, suggested improvements with P1-P3 tags, open questions)
9. References (BibTeX via `references.bib`)
10. Appendix: Detailed Results and Data
11. Appendix: Reproducibility (optimized parameters, waveforms, gate sequences, assumptions, procedure, saved artifacts)

## Main Text Prose Rules
- **No filenames or script references** in the main text. File references belong in the Reproducibility appendix only.
- **No code-style identifiers** (`snake_case`, `camelCase`, backtick-quoted names) in prose. Write "the dispersive shift" not "`chi_shift`".
- Every sentence should read as standard scientific writing (Physical Review style).
- Every equation must be numbered and referenced in the text.
- Every figure and table must be referenced in the running text.

## Equations
- Long equations must fit within column margins. Use `multline`, `split`, or `align` to break at binary operators.
- Fix any `Overfull \hbox` warnings related to equations.

## Figures
- Use `\includegraphics` with relative paths from the `report/` directory.
- Prefer `.pdf` (vector) format in compiled reports.
- Every figure needs: descriptive caption, labeled axes with units, legible font (8pt minimum), colorblind-friendly palette.

## Parameters
- Use `siunitx` for all values with units: `\SI{6.150}{GHz}`, `\SI{-2.84}{MHz}`.
- Tabulate system parameters with `booktabs` rules (`\toprule`, `\midrule`, `\bottomrule`).
