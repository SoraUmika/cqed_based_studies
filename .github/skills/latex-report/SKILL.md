---
name: latex-report
description: "Generate and compile a LaTeX study report. Use when: writing the final report, creating report.tex, compiling PDF, finishing a study, generating publication-quality documentation. Provides a 5-section LaTeX template with BibTeX support."
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

Fill in the 5 mandatory sections:

| Section | Content |
|---------|---------|
| **Motivation** | What was studied and why. For REP class, cite the original paper. |
| **System & Methods** | Hamiltonian (use LaTeX equations), parameters table, simulation approach, cqed_sim classes used. |
| **Results** | Figures (use `\includegraphics`), tables, quantitative findings. Reference all figures in `../figures/`. |
| **Discussion** | Physical interpretation, comparison to literature, limitations, error analysis. |
| **Conclusion** | Summary of findings, potential follow-up work. |

### 2. Include Figures

All figures must use relative paths from the report directory:

```latex
\includegraphics[width=\columnwidth]{../figures/fidelity_vs_drive_amplitude.pdf}
```

- Prefer `.pdf` figures for vector quality in the compiled report.
- Every figure must have a `\caption` and `\label`.
- Reference all figures in the text with `\ref`.

### 3. Add References

Use BibTeX for citations:

1. Create `studies/<study_name>/report/references.bib` with all cited papers.
2. Use `\cite{key}` in the text.
3. The template includes `\bibliographystyle{apsrev4-2}` and `\bibliography{references}`.

### 4. Compile to PDF

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

### 5. Finalize

After successful compilation:

1. Verify the PDF renders correctly (figures, equations, references).
2. Update the study README status to `COMPLETE`.

## Rules

- Never submit a report without completed validation.
- All figures must be referenced in the text — no orphan figures.
- Parameter values in the report must match those in the simulation scripts exactly.
