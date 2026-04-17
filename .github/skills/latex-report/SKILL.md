---
name: latex-report
description: "Generate and compile a LaTeX study report. Use when: writing the final report, creating report.tex, compiling PDF, finishing a study, generating publication-quality documentation. Enforces the full AGENTS.md report structure, appendix requirements, and self-review checks."
argument-hint: "Path to study folder, e.g. studies/transmon_chi_shift_optimization"
---

# Generate LaTeX Study Report

## When to Use

- Writing the final report for a completed study (AGENTS.md Step 5)
- User says "write the report", "generate PDF", "compile report"
- All validation checks have passed

## Prerequisites

Before generating the report, verify:

1. All goals in the study README are met
2. Validation checks are complete (use the `validate-results` skill if not)
3. Figures exist in `studies/<study_name>/figures/`
4. Data exists in `studies/<study_name>/data/`

## Procedure

### 1. Generate report.tex

Use the [report template](./assets/report_template.tex) as the starting point.

Fill in the full report structure:

| Section | Content |
|---------|---------|
| **Abstract** | Self-contained summary of motivation, methods, key result, and conclusion. |
| **Introduction** | Context, prior work, objective, and roadmap. |
| **System and Methods** | Hamiltonian, parameter table, analytic preliminary, computational approach, cqed_sim classes used. |
| **Results** | Figures, tables, quantitative findings, and interpretation. Reference all figures in `../figures/`. |
| **Validation** | Sanity checks, convergence, and literature comparison where applicable. |
| **Discussion** | Physical interpretation, comparison to literature, and limitations. |
| **Conclusion** | Summary of findings, potential follow-up work. |
| **Limitations and Future Work** | Known limitations, suggested improvements with `[P1|P2|P3]` and `[LOW|MEDIUM|HIGH]`, and open questions. |
| **Appendices** | `Detailed Results and Data` plus `Reproducibility`. |

### 2. Include Figures

All figures must use relative paths from the report directory:

```latex
\includegraphics[width=\columnwidth]{../figures/fidelity_vs_drive_amplitude.pdf}
```

- Prefer `.pdf` figures for vector quality in the compiled report.
- Every figure must have a `\caption` and `\label`.
- Reference all figures in the text with `\ref`.
- For OPT/DES studies, include time-domain and frequency-domain waveform plots in the appendix.

### 3. Add References

Use BibTeX for citations:

1. Create `studies/<study_name>/report/references.bib` with all cited papers.
2. Use `\cite{key}` in the text.
3. The template includes `\bibliographystyle{apsrev4-2}` and `\bibliography{references}`.

### 4. Write Reproducibility Appendix (MANDATORY)

Every report **must** include a `\section{Reproducibility}` appendix. This is required by AGENTS.md and is checked before marking a study COMPLETE.

Include these subsections:

| Subsection | Content |
|------------|---------|
| Optimized Parameters | Full table of every parameter that produced the final result |
| Waveform and Pulse Information | Time slices, dt, amplitude bounds; reference to artifact files |
| Gate Sequence and Decomposition | Exact ordered gate sequence with types, parameters, durations |
| Modeling and Simulation Assumptions | Hilbert space dims, solver tolerances, approximations |
| Reproduction Procedure | Step-by-step: which scripts to run, expected outputs, how to verify |
| Saved Artifacts | List of files in `artifacts/` and `data/` with format, contents, and load code |

Also ensure `studies/<study_name>/artifacts/` contains machine-readable files (JSON/NPZ/CSV) for key optimized results.

The `Saved Artifacts` subsection should use a table with `Filename`, `Format`, `Contents`, and `Load Example`. If a file is central to the study but is absent from that table, the appendix is incomplete.

### 5. Compile to PDF

Run from the `studies/<study_name>/report/` directory:

```
pdflatex report.tex
bibtex report
pdflatex report.tex
pdflatex report.tex
```

Or use `latexmk`:

```
latexmk -pdf report.tex
```

### 6. Finalize

After successful compilation:

1. Verify the PDF renders correctly (figures, equations, references).
2. Verify `artifacts/` directory contains saved results.
3. **Create `scripts/reproducibility_notebook.ipynb`** if it does not already exist (see AGENTS.md Step 6). The notebook is mandatory for every study.
4. Update the README `## Validation` section so it matches the actual report evidence.
5. Update the study README status to `COMPLETE`.

### 7. Mandatory Self-Review

Before declaring the report complete, check all of the following:

- No filenames or script names appear in the main text. Keep them in the Reproducibility appendix only.
- No code-style identifiers (`snake_case`, `camelCase`, backtick-quoted symbols) appear in prose.
- `Validation`, `Limitations and Future Work`, `Detailed Results and Data`, and `Reproducibility` sections are present.
- Every major claim in Results/Discussion points to a figure, table, or saved artifact.
- Every important equation is numbered and referenced.
- The PDF and log are free of unresolved layout problems.

## Rules

- Never submit a report without completed validation.
- All figures must be referenced in the text — no orphan figures.
- Do not put filenames, script names, or artifact paths in the main text.
- Parameter values in the report must match those in the simulation scripts exactly.
- Every report must include a Reproducibility appendix with artifacts.
