# cQED Autonomous Research Platform

Agents use the **cqed_sim** framework to simulate, optimize, analyze, and validate circuit quantum electrodynamics (cQED) problems — then produce publication-quality reports.

---

## Quick Reference (Start Here)

> **Read this section first.** It summarizes the entire document. Refer to later sections for details.

### 5 Non-Negotiable Rules

| # | Rule | Consequence of Violation |
|---|------|--------------------------|
| 1 | Always use `cqed_sim` — no ad-hoc simulation code unless a gap is documented | Results are not reproducible; study is invalid |
| 2 | Never skip workflow steps (Initialize → Plan → Implement → Validate → Report) | Incomplete studies; missing validation or documentation |
| 3 | Install packages as needed — use `pip install --user`, log in IMPROVEMENTS.md | Keep the environment lean; document what was added and why |
| 4 | No virtual environments — use system Python 3.12.10 directly | Breaks reproducibility across the platform |
| 5 | No notebooks unless the user explicitly requests them | Violates script-based workflow |

### Workflow at a Glance

```
Step 1: Initialize  →  Create study folder + README + IMPROVEMENTS.md
Step 2: Plan        →  Check API Reference, identify gaps, state assumptions
Step 3: Implement   →  Write scripts, save data, generate figures, update IMPROVEMENTS.md
Step 4: Validate    →  Sanity checks ✓  Convergence ✓  Literature comparison ✓
Step 5: Report      →  Write report.tex (with MANDATORY appendices) → compile PDF → COMPLETE
```

### Key Paths

| What | Where |
|------|-------|
| Study root | `studies/<study_name>/` |
| Scripts | `studies/<study_name>/scripts/` |
| Data | `studies/<study_name>/data/` |
| Figures | `studies/<study_name>/figures/` |
| Report | `studies/<study_name>/report/report.tex` |
| Improvement log | `studies/<study_name>/IMPROVEMENTS.md` |
| API Reference | [API_REFERENCE.md](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/API_REFERENCE.md) |

### DO / DON'T Checklist

**DO:**
- ✅ Consult the API Reference before writing any simulation code
- ✅ Classify the problem (OPT/REP/DES/ANA) before starting
- ✅ Update `IMPROVEMENTS.md` in real time during implementation
- ✅ Put findings and discussion in the **main text**; put raw data, optimal parameters, and detailed plots in **appendices**
- ✅ Validate results before writing the report (all 3 checks)
- ✅ Tag every suggested improvement with priority (P1/P2/P3) and difficulty (LOW/MEDIUM/HIGH)
- ✅ Save figures in both `.png` (300 dpi) and `.pdf` (vector) formats
- ✅ Install packages that improve the analysis — use `pip install --user` and log in IMPROVEMENTS.md

**DON'T:**
- ❌ Write ad-hoc simulation code when `cqed_sim` already has the functionality
- ❌ Skip the appendix — it is **required**, not optional
- ❌ Put raw pulse shapes, full parameter tables, or sweep data dumps in the main Results section
- ❌ Delete entries from `IMPROVEMENTS.md` — move resolved items to a `## Resolved` section
- ❌ Write vague limitation entries ("optimization didn't converge" → say **why** and **what to try next**)

---

## Critical Rules

> These rules are **non-negotiable**. Violating any of them constitutes a hard failure.

1. **Always use `cqed_sim`.** Never write ad-hoc simulation code unless a specific technical gap is documented in the study README. Consult the [API Reference](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/API_REFERENCE.md) before writing any simulation code.
2. **Never skip workflow steps.** Every study must follow the full lifecycle: Initialize → Plan → Implement → Validate → Report.
3. **Install packages as needed.** If a package would make the analysis more useful or complete, install it with `pip install <package> --user`. Prefer well-known libraries (numpy, scipy, matplotlib, qutip, lmfit, seaborn, etc.) but any package is acceptable if it genuinely benefits the study. Log every new installation in `IMPROVEMENTS.md` under `## Compute & Resource Notes`.
4. **No virtual environments.** Use system Python (Python 3.12.10) directly. No venvs, conda, or poetry unless the user explicitly requests one.
5. **No notebooks** unless the user explicitly requests notebook-based work.

---

## Problem Classes

Classify every new task into one or more classes before starting work.

| ID    | Class                    | Description                                                        | Typical Deliverables                                      |
|-------|--------------------------|--------------------------------------------------------------------|-----------------------------------------------------------|
| `OPT` | Parameter Optimization   | Optimize control parameters for a target unitary or state transfer | Converged parameters, fidelity plots, landscape scans     |
| `REP` | Result Reproduction      | Reproduce published results — spectra, dynamics, benchmarks        | Comparison plots, quantitative agreement metrics          |
| `DES` | Experiment Design        | Design state preparation, gate implementation, or measurement protocols | Pulse sequences, protocol specs, simulated outcomes  |
| `ANA` | System Analysis          | Extract physical insights, identify optimal operating points       | Parameter sweeps, phase diagrams, operating-point recommendations |

---

## Workflow

Every task follows these five steps in order.

### Step 1 — Initialize Study

Create the study folder and README:

```
studies/<descriptive_name>/
├── README.md          ← created in this step
├── IMPROVEMENTS.md    ← living log of limitations, ideas, and future work
├── scripts/           ← simulation & analysis code
├── data/              ← raw and processed outputs
├── figures/           ← plots for the report
└── report/
    ├── report.tex     ← final LaTeX report
    ├── references.bib ← BibTeX bibliography
    └── report.pdf     ← compiled PDF
```

The README **must** contain these sections:

```markdown
# <Study Title>

## Problem Class
<!-- OPT | REP | DES | ANA — pick one or more -->

## Motivation
<!-- Why this study matters. Link to paper if REP class. -->

## Goals
<!-- Numbered, concrete, falsifiable goals. -->

## Methods
<!-- Which cqed_sim modules/functions will be used. -->

## Expected Outcomes
<!-- What success looks like — quantitative where possible. -->

## Known Limitations
<!-- Updated throughout the study. What approximations are being made?
     What is constrained by compute time or framework capability?
     This section feeds directly into the report and IMPROVEMENTS.md. -->

## Status
<!-- ACTIVE | COMPLETE | BLOCKED — update as work progresses. -->
```

#### IMPROVEMENTS.md — The Improvement Log

Every study must maintain an `IMPROVEMENTS.md` file. This is a **living document** — the agent updates it throughout the study, not just at the end. It serves as the bridge between the current work and all future work on the same topic.

```markdown
# Improvement Log: <Study Title>

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
<!-- Things that could make the current results qualitatively wrong. -->
<!-- Format: - **<What>**: <Why it matters>. <What to do about it>. -->

## Recommended Improvements (P2)
<!-- Things that would meaningfully improve accuracy or scope. -->

## Nice-to-Haves (P3)
<!-- Lower-priority enhancements. -->

## Open Questions
<!-- Unresolved physics or numerical observations worth investigating. -->

## What Was Tried and Did Not Work
<!-- Failed approaches, dead-end parameter ranges, algorithms that
     diverged — anything that saves the next agent from repeating
     the same mistakes. Include enough detail to understand WHY
     it failed, not just that it did. -->

## Compute & Resource Notes
<!-- Wall-clock times for key simulations. Memory usage.
     Which runs were the bottleneck. Helps future agents
     plan their compute budget. -->
```

**Rules for maintaining IMPROVEMENTS.md:**

1. **Start it in Step 1** (Initialize), even if it only has placeholder headings.
2. **Update it in real time** during Steps 3–4 (Implement & Validate). When you hit a limitation, log it immediately — don't wait until the report.
3. **Never delete entries.** If an issue is resolved, move it to a `## Resolved` section at the bottom with a note on how it was fixed.
4. **Be specific about failures.** "Optimization didn't converge" is useless. "Nelder-Mead on 12-parameter pulse stalled at fidelity 0.987 after 500 iterations; cost landscape appears flat near the minimum; GRAPE or gradient-based method likely needed" is useful.

### Step 2 — Plan & Validate Approach

Before writing simulation code:

1. **Check the [API Reference](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/API_REFERENCE.md)** — confirm the required functionality exists in `cqed_sim`.
2. **Identify gaps** — if `cqed_sim` cannot handle the task, document the gap in the README and explain why standalone code is needed.
3. **State assumptions** — list all physical assumptions, parameter ranges, and convergence criteria in the README.

> If new reusable functionality is developed, add a `## Suggested Upstreaming` section to the README.

### Step 3 — Implement & Execute

- Write scripts in `studies/<study_name>/scripts/`.
- Save raw numerical output to `studies/<study_name>/data/`.
- Generate publication-quality figures to `studies/<study_name>/figures/`.
- **Update `IMPROVEMENTS.md` in real time** — log limitations, failed approaches, and compute notes as they arise. Do not defer this to the report phase.
- Update README status as work progresses.

### Step 4 — Validate Results

Before reporting, **all three checks must pass**:

- [ ] **Sanity checks** — Verify limiting cases, conservation laws, or known analytic results.
- [ ] **Convergence** — Confirm results are stable with respect to truncation (Hilbert space dimension), time steps, and optimization iterations.
- [ ] **Literature comparison** (if applicable) — Quantitatively compare to published benchmarks; report percent error or fidelity.

After validation, **finalize `IMPROVEMENTS.md`**: review every limitation discovered during implementation and validation, ensure suggested improvements have priority and difficulty tags, and record any open questions that emerged.

### Step 5 — Report

Write `report/report.tex` following the **Scientific Review Paper Format** defined below. The report **must** include:

1. A `Limitations and Future Work` section — content drawn from `IMPROVEMENTS.md`.
2. A **mandatory `Appendix`** containing detailed results and data (optimal pulses, full parameter tables, sweep data, etc.). The main text presents findings and discussion; the appendix presents the supporting data.

Then: compile to PDF → update README status to `COMPLETE`.

---

## Report Format: Scientific Review Paper

All reports must follow this structure. The format mirrors a peer-reviewed journal article: abstract at the top, structured body sections, and references at the bottom.

### Required LaTeX Template

```latex
\documentclass[11pt, a4paper]{article}

% ── Packages ──────────────────────────────────────────────
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath, amssymb, amsfonts}
\usepackage{graphicx}
\usepackage[margin=1in]{geometry}
\usepackage{hyperref}
\usepackage[numbers, sort&compress]{natbib}
\usepackage{booktabs}
\usepackage{siunitx}
\usepackage{caption}
\usepackage{float}
\usepackage{xcolor}

% ── Metadata ──────────────────────────────────────────────
\title{\textbf{<Study Title>}}
\author{<Author(s)> \\ <Affiliation>}
\date{\today}

\begin{document}

\maketitle

% ══════════════════════════════════════════════════════════
% 1. ABSTRACT
% ══════════════════════════════════════════════════════════
\begin{abstract}
% 150–300 words. Summarize the motivation, methods, key results,
% and primary conclusion. Must be self-contained — a reader
% should understand the scope and outcome without reading further.
\end{abstract}

% ══════════════════════════════════════════════════════════
% 2. INTRODUCTION
% ══════════════════════════════════════════════════════════
\section{Introduction}
% - Context: physical system and its significance in cQED.
% - Problem statement: what question or challenge is addressed.
% - Prior work: brief survey of relevant literature with citations.
% - Objective: what this study aims to achieve.
% - Outline: one-sentence roadmap of remaining sections.

% ══════════════════════════════════════════════════════════
% 3. SYSTEM & METHODS
% ══════════════════════════════════════════════════════════
\section{System and Methods}

\subsection{Hamiltonian}
% Write the full Hamiltonian in equation form.
% Define every symbol and state its numerical value with units.

\subsection{Simulation Parameters}
% Tabulate all parameters:
% \begin{table}[H]
%   \centering
%   \caption{System parameters used in this study.}
%   \begin{tabular}{@{} l S[table-format=+1.3] l @{}}
%     \toprule
%     Parameter & {Value} & Unit \\
%     \midrule
%     ...
%     \bottomrule
%   \end{tabular}
% \end{table}

\subsection{Computational Approach}
% - Which cqed_sim modules/functions were used (with version).
% - Hilbert space truncation dimensions.
% - Time-stepping scheme and step size.
% - Optimization algorithm and convergence criteria (if OPT class).
% - Any standalone code and justification (link to README gap analysis).

% ══════════════════════════════════════════════════════════
% 4. RESULTS
% ══════════════════════════════════════════════════════════
\section{Results}
% Present key findings with figures and tables.
% Focus on INSIGHTS and INTERPRETATION here — what the results
% mean, not just what they are.
%
% IMPORTANT: Detailed raw data (optimal pulse shapes, full parameter
% tables, sweep data dumps, cost landscapes) belong in the APPENDIX,
% not in this section. Include summary figures/tables here and
% reference the appendix for full details.
%
% Every figure must have a descriptive caption and labeled axes with units.
% Every table must have a caption.
% Reference all figures/tables in the running text.
%
% Subsection structure should mirror the study goals:
%   \subsection{Goal 1: <description>}
%   \subsection{Goal 2: <description>}

% ══════════════════════════════════════════════════════════
% 5. VALIDATION
% ══════════════════════════════════════════════════════════
\section{Validation}
% Document each of the three required checks:
%   \subsection{Sanity Checks}
%   \subsection{Convergence Analysis}
%   \subsection{Literature Comparison}  % if applicable

% ══════════════════════════════════════════════════════════
% 6. DISCUSSION
% ══════════════════════════════════════════════════════════
\section{Discussion}
% - Physical interpretation of results.
% - Comparison to literature (quantitative where possible).
% - Limitations and assumptions that affect generality.
% - Unexpected findings or open questions.

% ══════════════════════════════════════════════════════════
% 7. CONCLUSION
% ══════════════════════════════════════════════════════════
\section{Conclusion}
% - Concise summary of key findings.
% - Whether each study goal was met (cross-reference the README).

% ══════════════════════════════════════════════════════════
% 8. LIMITATIONS & FUTURE WORK
% ══════════════════════════════════════════════════════════
\section{Limitations and Future Work}
% THIS SECTION IS MANDATORY. It serves as the primary handoff
% mechanism between the current agent and future agents.
%
% Write this section with the assumption that a future agent
% will read it to decide what to do next. Be specific, honest,
% and actionable.

\subsection{Known Limitations}
% For each limitation, state:
%   (a) What the limitation is.
%   (b) Why it exists (e.g., compute budget, cqed_sim capability,
%       Hilbert space truncation, optimizer choice).
%   (c) How it affects the reliability or generality of the results.
%
% Examples of things to flag:
% - Optimization may be stuck in a local minimum.
% - Hilbert space truncation too aggressive for high-photon-number states.
% - Only one optimizer was tested; others may perform better.
% - Decoherence / noise channels were omitted.
% - Parameter sweep was coarse-grained due to compute time.
% - Rotating-wave or dispersive approximation may break down at
%   the edges of the explored parameter range.

\subsection{Suggested Improvements}
% Concrete, actionable next steps a future agent should take.
% For each suggestion, include:
%   (a) What to do.
%   (b) Why it would improve the study.
%   (c) Estimated difficulty (LOW / MEDIUM / HIGH).
%   (d) Priority (P1 = critical gap, P2 = meaningful improvement,
%       P3 = nice to have).
%
% Format as an itemized list for machine readability:
% \begin{itemize}
%   \item \textbf{[P1 | HIGH]} Re-run optimization with GRAPE in
%         addition to Nelder-Mead to check for local minima.
%         Current fidelity (99.2\%) may not be the global optimum.
%   \item \textbf{[P2 | LOW]} Increase cavity Fock-state truncation
%         from $N=8$ to $N=15$ and verify convergence of $\chi'$ shift.
%   \item \textbf{[P3 | MEDIUM]} Include $T_1$ and $T_\phi$ decay
%         channels to assess gate fidelity under realistic noise.
% \end{itemize}

\subsection{Open Questions}
% Physics questions or unexpected observations that emerged
% during the study but were not resolved. Future agents or
% the user may choose to investigate these.

% ══════════════════════════════════════════════════════════
% 9. REFERENCES
% ══════════════════════════════════════════════════════════
\bibliographystyle{unsrtnat}
\bibliography{references}

% ══════════════════════════════════════════════════════════
% APPENDICES (REQUIRED)
% ══════════════════════════════════════════════════════════
% The appendices are MANDATORY. They contain the detailed data
% that supports the main text. The main body presents findings,
% interpretation, and discussion; the appendices present the
% raw/detailed results that a reader (or future agent) needs
% to reproduce or extend the work.
%
% Include whichever appendix subsections are relevant to the
% problem class. See the "Appendix Content by Problem Class"
% table in AGENTS.md for guidance.

\appendix

\section{Detailed Results and Data}
% This is the primary appendix. Include the raw results that
% back up the main-text findings. Examples:
%
% For OPT studies:
%   - Optimal pulse shapes (I/Q vs time plots)
%   - Full parameter tables for converged solutions
%   - Cost function / fidelity convergence traces
%   - Optimization landscape cross-sections
%
% For REP studies:
%   - Side-by-side comparison plots (simulation vs published)
%   - Full numerical comparison tables with percent errors
%
% For DES studies:
%   - Complete pulse sequence diagrams with timing
%   - Protocol parameter tables
%   - State tomography / process tomography matrices
%
% For ANA studies:
%   - Full parameter sweep heatmaps or line cuts
%   - Phase diagram data
%   - Extracted fit parameters with uncertainties

% \section{Supplementary Derivations}
% Optional — include if the study required non-trivial
% analytic work that supports but does not belong in the
% main text.

% \section{Simulation Configuration}
% Optional — include full simulation config (Hilbert space
% dimensions, solver tolerances, time grids) if they are
% too detailed for the Methods section.

\end{document}
```

### Report Section Rules

| Section | Required? | Guidelines |
|---------|-----------|------------|
| Abstract | **Yes** | 150–300 words. Self-contained summary: motivation, methods, key results, conclusion. No citations, no equations, no acronyms without definition. |
| Introduction | **Yes** | Establish context, cite prior work, state objectives. End with a one-sentence roadmap. |
| System and Methods | **Yes** | Full Hamiltonian with all terms defined. Parameters in a table with units (use `siunitx`). Identify `cqed_sim` modules used. |
| Results | **Yes** | One subsection per study goal. All figures captioned with labeled axes and units. Reference every figure and table in the text. |
| Validation | **Yes** | Report all three checks (sanity, convergence, literature). Include convergence plots where applicable. |
| Discussion | **Yes** | Interpret results physically. Compare quantitatively to literature. State limitations honestly. |
| Conclusion | **Yes** | Summarize findings. State whether each goal was met. |
| Limitations & Future Work | **Yes** | **The agent-to-agent handoff section.** Must include Known Limitations (with cause and impact), Suggested Improvements (with priority P1–P3 and difficulty), and Open Questions. See detailed guidance below. |
| References | **Yes** | Use BibTeX (`references.bib`). Cite all referenced papers, the `cqed_sim` framework, and any external tools. Use `unsrtnat` style (numbered, order of appearance). |
| Appendices | **Yes** | **Required for all studies.** Must include a "Detailed Results and Data" section with raw outputs: optimal pulse shapes, converged parameters, sweep data, comparison tables, etc. Main text presents findings and discussion; appendices present the supporting data. See **Appendix Content Guide** below. |

### Main Text vs. Appendix — Content Split Rule

The main text (Results, Discussion) should present **findings, interpretation, and physical insight**. The appendices should present the **underlying data** that a reader or future agent needs to reproduce, verify, or extend the work.

| Content Type | Goes In | Example |
|--------------|---------|---------|
| Key result (e.g., "fidelity reached 99.5%") | **Main text** (Results) | "The optimized gate achieves $\mathcal{F} = 0.995$." |
| Physical interpretation | **Main text** (Discussion) | "The fidelity plateau is caused by leakage to the $|2\rangle$ state." |
| Optimal pulse shape (I/Q vs. time) | **Appendix** | Figure: `optimal_pulse_iq.pdf` |
| Full converged parameter table | **Appendix** | Table with all 12 optimized pulse coefficients |
| Summary parameter table (key values only) | **Main text** (Results) | Table with gate time, fidelity, leakage |
| Convergence trace (fidelity vs. iteration) | **Main text** (Validation) or **Appendix** | If brief, in Validation; if many traces, in Appendix |
| Parameter sweep heatmap (full) | **Appendix** | 2D grid of $\chi$ vs. $\kappa$ showing fidelity |
| Selected sweep cross-sections | **Main text** (Results) | 1D line cuts at optimal operating point |
| Cost landscape / optimization landscape | **Appendix** | Contour plots of cost function |
| Side-by-side comparison with literature | **Main text** (Validation) | Summary comparison; **Appendix** for full data tables |
| Pulse sequence timing diagrams | **Appendix** | Detailed protocol spec with all gates and delays |

#### Appendix Content by Problem Class

| Problem Class | Required Appendix Content |
|---------------|---------------------------|
| **OPT** (Optimization) | Optimal pulse shapes (I/Q plots), full parameter tables, convergence traces, cost landscape cross-sections |
| **REP** (Reproduction) | Full numerical comparison tables (simulation vs. published), overlay plots, percent-error breakdowns |
| **DES** (Experiment Design) | Complete pulse sequence diagrams, protocol parameter tables, simulated measurement outcomes, state/process tomography data |
| **ANA** (System Analysis) | Full sweep heatmaps, phase diagram data, extracted fit parameters with uncertainties, all line cuts |

> **Rule of thumb:** If a figure or table shows *what the answer is* (detailed data), it belongs in the appendix. If it shows *what the answer means* (insight, comparison, trend), it belongs in the main text. When in doubt, put it in both — a summary in the main text referencing the full version in the appendix.

### Limitations & Future Work — Detailed Guidance

This section exists so that **future agents can pick up where the current agent left off** without re-discovering the same dead ends or missing the same gaps. It is not boilerplate — it is the most practically important section for iterative research.

**Known Limitations** must cover every simplification or constraint that could affect the results. Common categories in cQED studies:

| Category | Examples |
|----------|----------|
| Optimization | Local minima not ruled out; only one algorithm tested; cost landscape not fully explored |
| Truncation | Hilbert space dimension may be too low for high-photon states; Fock basis cutoff not converged |
| Physics omitted | No decoherence channels; no thermal photon population; RWA or dispersive approximation used |
| Compute budget | Parameter sweep was coarse; longer pulse durations not tested; wall-clock time limit reached |
| Framework gaps | `cqed_sim` does not support feature X; workaround Y was used instead |

**Suggested Improvements** must be tagged with priority and difficulty:

| Tag | Meaning |
|-----|---------|
| **P1** | Critical gap — results may be qualitatively wrong without this fix |
| **P2** | Meaningful improvement — results are qualitatively correct but quantitatively limited |
| **P3** | Nice to have — would strengthen the study but not essential |
| **LOW** | Can be done by changing a parameter or re-running a script |
| **MEDIUM** | Requires writing new code or modifying the simulation setup |
| **HIGH** | Requires new `cqed_sim` functionality, significant compute, or new physics |

**Open Questions** should capture anything surprising or unresolved — anomalous data points, unexpected parameter sensitivities, or physical phenomena that deserve their own study.

### Citation and Reference Guidelines

- **Bibliography file:** Every study must include a `references.bib` in the `report/` directory.
- **Citation style:** Numbered references in order of appearance (`unsrtnat` with `natbib`).
- **Minimum citations:** Cite at least (a) the original papers for any reproduced results, (b) the `cqed_sim` framework, and (c) foundational references for the physical system (e.g., transmon, Jaynes–Cummings model).
- **BibTeX entry format:** Use `@article` for journal papers, `@misc` for preprints/software. Always include `doi` when available.

Example `references.bib` entry:

```bibtex
@article{koch2007transmon,
  author  = {Koch, Jens and Yu, Terri M. and Gambetta, Jay and others},
  title   = {Charge-insensitive qubit design derived from the {Cooper} pair box},
  journal = {Physical Review A},
  volume  = {76},
  pages   = {042319},
  year    = {2007},
  doi     = {10.1103/PhysRevA.76.042319}
}

@misc{cqed_sim,
  author = {<Author(s)>},
  title  = {cqed\_sim: Circuit QED Simulation Framework},
  year   = {2025},
  url    = {https://github.com/SoraUmika/qubox_cQEDsim}
}
```

### Figure Standards for Reports

- Save both `.png` (300 dpi, for quick inspection) and `.pdf` (vector, for LaTeX inclusion).
- Use `\includegraphics` with relative paths from the `report/` directory.
- Every figure must have: descriptive caption, labeled axes with units, legible font size (≥ 8 pt in print), and a colorblind-friendly palette.
- Prefer vector formats (`.pdf`) in the compiled report.

### Compilation

```bash
cd studies/<study_name>/report/
pdflatex report.tex
bibtex report
pdflatex report.tex
pdflatex report.tex
```

Or, equivalently: `latexmk -pdf report.tex`.

---

## cqed_sim Framework

| Resource              | Location |
|-----------------------|----------|
| Source code (GitHub)  | [SoraUmika/qubox_cQEDsim](https://github.com/SoraUmika/qubox_cQEDsim) |
| API Reference         | [API_REFERENCE.md](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/API_REFERENCE.md) |
| Physics Conventions   | [physics_conventions_report.tex](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/physics_and_conventions/physics_conventions_report.tex) |
| Local copy            | `C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation` (most up-to-date version) |

**Usage policy:**

- Always prefer `cqed_sim` for simulation, modeling, and experiment reproduction.
- Do not duplicate functionality that already exists in the framework.
- Document any deviation in the study README.

---

## Decision Trees

### Can I use standalone simulation code?

```
Is the required functionality in cqed_sim?
├─ YES → Use cqed_sim. Do NOT write standalone code.
├─ PARTIALLY → Extend using cqed_sim as the foundation.
│              Document what is missing in the README.
└─ NO → Write standalone code.
        Document the gap in the README.
        Add a "Suggested Upstreaming" section if reusable.
```

### Should I install a package?

```
Is the package already installed?
├─ YES → Use it.
└─ NO → Would it genuinely improve the analysis or study?
         ├─ NO → Do not install.
         └─ YES → pip install --user <package>
                  Log the install in IMPROVEMENTS.md (Compute & Resource Notes).
```

### Is the study complete?

```
Are all goals in the README met?
├─ NO → Continue work or document blockers.
└─ YES → Have results been validated (sanity, convergence, literature)?
          ├─ NO → Complete validation first.
          └─ YES → Is IMPROVEMENTS.md populated with limitations, suggestions, and failed approaches?
                    ├─ NO → Complete the improvement log first.
                    └─ YES → Write report.tex:
                              ─ Main text: findings, interpretation, discussion
                              ─ Appendix: detailed data (pulses, parameters, sweeps)
                              ─ Limitations & Future Work section
                             → Compile PDF → README status = COMPLETE
```

---

## Conventions

### Naming

- Study folders: `studies/<lowercase_descriptive_name>/` (e.g., `studies/transmon_chi_shift_optimization/`)
- Scripts: `snake_case.py`
- Figures: `<figure_description>.png` or `.pdf`

### Code Quality

- All scripts must be self-contained and runnable from the study folder.
- Include docstrings explaining the physical setup and parameters.
- Pin key simulation parameters (Hilbert space dims, time steps, etc.) as named constants at the top of each script.

### Figures

- Use matplotlib with consistent styling (labeled axes, units, legends).
- Save both `.png` (for README / quick inspection) and `.pdf` (for LaTeX report).
- Use colorblind-friendly palettes where possible (e.g., `tab10`, `colorblind` from `seaborn`).

---

## Continuous Research Loop

> For detailed architecture, see [RESEARCH_LOOP.md](RESEARCH_LOOP.md).

### Overview

The platform supports a **two-model continuous research loop** for autonomous studies:

| Model | Role | Agent | When Used |
|-------|------|-------|-----------|
| **Science Director** (Codex/GPT) | Physics reasoning, experiment design, result review | `@science-director` | Planning and review phases |
| **Execution Engineer** (Opus) | Implementation, debugging, documentation, reporting | `@execution-engineer` | All execution phases |
| **Research Loop** (combined) | Single-agent mode that switches between both roles | `@research-loop` | Click-and-research mode |

### Quick Start

**New study (click-and-research):**
```
@research-loop study=studies/my_study goal='Optimize chi for 99.5% readout fidelity'
```

**New study (step-by-step):**
```powershell
# 1. Initialize
.\tools\research_loop.ps1 -Action init -StudyName "my_study" -StudyGoal "Optimize chi"

# 2. Plan (Science Director)
@science-director study=studies/my_study run=task_runs/my_study phase=plan

# 3. Execute (Execution Engineer)
@execution-engineer study=studies/my_study run=task_runs/my_study phase=implement

# 4. Review (Science Director)
@science-director study=studies/my_study run=task_runs/my_study phase=review

# 5. Repeat 3-4 until VALIDATE decision, then:
@execution-engineer study=studies/my_study run=task_runs/my_study phase=validate
@execution-engineer study=studies/my_study run=task_runs/my_study phase=report
```

**Resume interrupted study:**
```
@research-loop study=studies/my_study resume
```

**Check status:**
```powershell
.\tools\research_loop.ps1 -Action status -StudyName "my_study"
```

### Loop Phases

```
BOOTSTRAP → PLAN → IMPLEMENT → REVIEW ─┐
                       ↑                │
                       └── CONTINUE ────┘
                           REVISE ──────┘
                       VALIDATE → REPORT → COMPLETE
```

### State Files

| File | Location | Purpose |
|------|----------|---------|
| `study_state.json` | `studies/<name>/` | Machine-readable study state (single source of truth) |
| `SCIENCE_DIRECTIVE.md` | `task_runs/<name>/` | Science Director → Execution Engineer instructions |
| `EXECUTION_SUMMARY.md` | `task_runs/<name>/` | Execution Engineer → Science Director results |
| `TASK_CHECKLIST.md` | `task_runs/<name>/` | Task completion tracking |
| `PROGRESS_LOG.md` | `task_runs/<name>/` | Append-only log of what happened |
| `BLOCKERS.md` | `task_runs/<name>/` | Active and resolved blockers |

### VS Code Tasks

Available from **Terminal → Run Task**:
- **Research: New Study** — Initialize a new study with goal
- **Research: Study Status** — Show current loop state
- **Research: Resume Study** — Detect phase and continue
- **Research: Run Loop Action** — Pick any phase to run

---

## Typical cQED System Parameters

Use these as default starting values unless the study specifies otherwise.

| Parameter                        | Symbol          | Value       | Unit  |
|----------------------------------|-----------------|-------------|-------|
| Dispersive shift                 | χ               | −2.84       | MHz   |
| Second-order dispersive shift    | χ′              | −21         | kHz   |
| Cavity self-Kerr                 | K               | −28         | kHz   |
| Qubit anharmonicity              | α               | −255        | MHz   |
| Qubit frequency                  | ω_q             | 6.150       | GHz   |
| Cavity frequency                 | ω_c             | 5.241       | GHz   |
| Readout resonator frequency      | ω_r             | 8.597       | GHz   |
| Readout resonator linewidth      | κ_r             | 2.4         | MHz   |
