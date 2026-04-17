---
name: report-preflight
description: "Automated pre-compilation scan of report.tex for common violations: filenames in main text, snake_case identifiers in prose, missing label/ref pairs, unreferenced figures, undefined citations, and overfull hbox patterns. Use when: before compiling the report, after writing report.tex, or as part of the self-review pass."
argument-hint: "Path to study folder, e.g. studies/transmon_chi_shift_optimization"
---

# Report Preflight Check

## When to Use

- After writing or editing `report/report.tex`
- Before compiling the report to PDF
- As part of the mandatory self-review pass (AGENTS.md Step 5)
- Before writing REVIEW_REQUEST.md

## Procedure

Read the full `report/report.tex` file and run each check below. Report results as PASS or FAIL.

### Check 1: No Filenames in Main Text

Search the text **before** `\appendix` (or `\section{Detailed Results}` if no `\appendix` command) for patterns that indicate filename references:

- File extensions: `.json`, `.npz`, `.csv`, `.py`, `.npy`, `.pdf`, `.png`, `.tex`, `.bib`
- Path separators: `data/`, `artifacts/`, `scripts/`, `figures/`
- Script names: any word ending in `.py`

**Exception:** `\includegraphics` paths and `\bibliography` commands are allowed.
**Exception:** Content inside `\begin{verbatim}...\end{verbatim}` or `\begin{lstlisting}...\end{lstlisting}` is allowed.

For each violation found, report the line number and the offending text.

### Check 2: No Code-Style Identifiers in Prose

Search main text (before `\appendix`) for:

- `snake_case` patterns: words containing underscores that are NOT inside math mode (`$...$`, `\(...\)`, equation environments) and NOT LaTeX commands (starting with `\`)
- `camelCase` patterns: words with internal capitals that are not proper nouns, acronyms, or standard physics terms
- Backtick-quoted identifiers: `` `anything` ``

For each violation, suggest the scientific-language replacement.

### Check 3: Label/Reference Consistency

- Extract all `\label{...}` definitions
- Extract all `\ref{...}`, `\eqref{...}`, `\cref{...}` references
- Report any reference to an undefined label
- Report any label that is never referenced (warning, not error)

### Check 4: Figure/Table References

- Extract all `\begin{figure}` and `\begin{table}` environments
- Check each has a `\label{...}` inside it
- Check each label is referenced somewhere in the text with `\ref` or equivalent
- Check each `\includegraphics{...}` path points to a file that exists in `figures/`

### Check 5: Citation Completeness

- Extract all `\cite{...}` and `\citep{...}` and `\citet{...}` keys
- Read `references.bib` and extract all BibTeX entry keys
- Report any cited key not defined in the `.bib` file
- Report any `.bib` entry never cited (warning)

### Check 6: Mandatory Sections Present

Check the report contains these section commands (case-insensitive on the title):
- `\begin{abstract}`
- `\section{Introduction}`
- `\section{System and Methods}` (or `\section{Methods}`)
- `\section{Results}`
- `\section{Validation}`
- `\section{Discussion}`
- `\section{Conclusion}`
- `\section{Limitations and Future Work}` (or `\section{Limitations}`)
- `\appendix`
- A section after `\appendix` containing "Reproducibility"

### Check 7: Equation Numbering

- Every `\begin{equation}` should have a `\label{...}`
- Every equation label should be referenced with `\eqref{...}` or `\ref{...}`
- Warn on unnumbered equation environments (`equation*`, `align*`) — these are acceptable but should be the minority

## Output Format

```markdown
## Report Preflight: <study_name>

**Overall:** PASS / FAIL (N issues found)

### Errors (must fix)
| # | Check | Line | Issue | Suggested Fix |
|---|-------|------|-------|---------------|
| 1 | Filenames | L42 | `grape_result.npz` in main text | Move to Reproducibility appendix |
| 2 | Code IDs | L78 | `chi_shift` in prose | Replace with "the dispersive shift $\chi$" |

### Warnings (should review)
| # | Check | Issue |
|---|-------|-------|
| 1 | Unused label | `fig:sweep_full` defined but never referenced |

### Summary
- Filenames in main text: N violations
- Code identifiers in prose: N violations
- Undefined references: N
- Missing figure files: N
- Undefined citations: N
- Missing sections: N
```
