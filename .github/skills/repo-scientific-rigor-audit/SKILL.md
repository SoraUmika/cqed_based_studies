---
name: repo-scientific-rigor-audit
description: "Recurring repo-health and scientific-rigor audit for autonomous cQED studies. Use when: auditing recent studies, checking whether reports and notebooks are trustworthy, reviewing AGENTS.md or skills for drift, verifying reproducibility and artifact traceability, or assessing autonomous-agent readiness after major studies, refactors, or workflow changes."
argument-hint: "Optional scope, for example: 'full repo audit', 'recent OPT studies', or 'post-refactor autonomy audit'"
---

# Repo Scientific-Rigor Audit

## Purpose

This recurring audit checks whether the repository is actually ready for autonomous agent-driven scientific experimentation and reporting.

The standard is not "can an agent run code?" The standard is whether an agent can:

- understand the scientific intent of a study,
- preserve assumptions, conventions, and scope,
- validate results before concluding,
- report uncertainty and failure modes honestly,
- leave behind reproducible artifacts and reports,
- and hand work to another agent or human without guesswork.

## When to Use

- During periodic repo-health reviews
- After any major study is completed or revised
- After workflow refactors or automation changes
- After updating AGENTS.md, RESEARCH_LOOP.md, research_config.json, or any core skill
- When study quality, reproducibility, or agent reliability may have drifted
- Before declaring the repo ready for broader autonomous use

## Audit Scope

### Recent Studies

Inspect the most recent scientifically relevant studies, not only the easiest ones.

At minimum:

1. Inspect the 5-10 most recently modified study folders under studies/.
2. Include any study with recent review-loop activity in task_runs/:
   - REVIEW_REQUEST.md
   - REVIEW_DIRECTIVE.md
   - FOLLOWUP_PROMPT.md
   - EXECUTION_SUMMARY.md
3. Include any study whose state is not fully settled:
   - missing study_state.json
   - README says COMPLETE but machine state does not
   - validation is unchecked or ambiguous
4. Exclude probes or demos only if they are clearly non-scientific. If they affect workflow reliability, assess them separately as workflow-health evidence.

### Repo Guidance Surface

Always inspect:

- AGENTS.md
- RESEARCH_LOOP.md
- research_config.json
- LESSONS_LEARNED.md
- core skills related to planning, validation, reporting, notebooks, and review

At minimum, inspect these skills when present:

- study-init
- validate-results
- study-validator
- report-review
- reproducibility-notebook
- latex-report
- report-preflight
- red-green-validation
- cqed-sim-lookup
- cross-study-analysis

## Required Evidence Sources

For each audited study, inspect as many of the following as exist:

- study_state.json
- README.md
- IMPROVEMENTS.md
- report/report.tex
- report/report.pdf
- scripts/reproducibility_notebook.ipynb
- artifacts/
- data/
- figures/
- task_runs/<study>/SCIENCE_DIRECTIVE.md
- task_runs/<study>/EXECUTION_SUMMARY.md
- task_runs/<study>/REVIEW_REQUEST.md
- task_runs/<study>/REVIEW_DIRECTIVE.md
- task_runs/<study>/FOLLOWUP_PROMPT.md
- task_runs/<study>/BLOCKERS.md

Do not rely on README claims alone. Cross-check claims against reports, artifacts, notebooks, and task-run state.

## Required Checks

### 1. Study Intent Legibility

Can another agent identify the scientific question, problem class, assumptions, approximations, and success criteria without guessing?

Check whether the study clearly states:

- what question is being answered,
- what is exploratory versus concluded,
- what counts as success,
- and what approximations or conventions govern the results.

### 2. Completion and State Consistency

Check whether README.md, study_state.json, and task-run files agree.

Flag any study where:

- README says COMPLETE while review is still pending,
- validation remains unchecked,
- the task-run shows REVISE or NEEDS_REWORK,
- or study_state.json is missing or stale.

In an autonomous-agent setting, stale state is a high-risk defect because it misleads the next agent about whether the study can be trusted or extended.

### 3. Evidence-to-Claim Integrity

Identify the study's headline claims and verify that figures, tables, or artifacts actually support them.

Pay special attention to claims about:

- optimality,
- impossibility or no-go results,
- robustness,
- convergence,
- agreement with literature,
- and device-specific recommendations.

Flag claims that are unsupported, optimizer-dependent, internally inconsistent, or stronger than the evidence allows.

### 4. Validation Coverage

Verify whether the study provides all required validation modes:

- sanity checks,
- convergence checks,
- literature comparison or analytic comparison.

Check whether:

- the reported numbers are explicit,
- tolerances are justified rather than arbitrary,
- failure regimes are described,
- uncertainty or restart variability is reported,
- and the report distinguishes exploratory evidence from validated conclusions.

### 5. Reproducibility and Notebook Quality

Do not treat notebook existence as sufficient.

Check whether the reproducibility notebook is structured for reliable handoff:

- an early tunable-parameter cell,
- a derived-objects cell,
- a default load-saved-results path,
- a rerun path,
- validation reproduction,
- and explanatory markdown before code cells.

Flag notebooks that exist but require hidden assumptions, stale paths, or manual code editing to reproduce the study.

### 6. Artifact Traceability

Check whether headline report numbers map to a single canonical artifact field or table row.

Flag cases where:

- artifacts disagree,
- the report does not say which artifact is canonical,
- or a reader cannot trace a reported claim back to saved machine-readable outputs.

### 7. Separation of Exploration from Conclusions

Check whether the study clearly distinguishes:

- exploratory investigation,
- partial evidence,
- validated conclusions,
- and final recommendations.

Autonomous agents are especially prone to turning exploratory plots into final claims. Treat blurred boundaries here as risky or misleading.

### 8. Honest Uncertainty and Limitations

Check whether assumptions, approximation validity, optimizer limits, parameter provenance, framework gaps, and failure modes are stated explicitly and not buried.

Treat soft language that hides unresolved issues as a risk.

### 9. Prompt-to-Result Provenance

Trace whether the research prompt or directive flows into:

- README goals,
- scripts,
- artifacts,
- report claims,
- and review documents.

Flag missing links that would make continuation difficult for another agent.

### 10. Instruction and Skill Sufficiency

Check whether AGENTS.md and the active skills truly enforce the repo's claimed standards.

Look for rules that are documented but not operationalized, such as:

- required validation with weak or optional checks,
- notebook structure described but not verified,
- self-review required but not gated,
- and state transitions described without enforcement.

### 11. Guardrail Drift

Compare AGENTS.md requirements against what templates, validators, and reviewer skills actually check.

Classify any mismatch as a systemic risk if it could let an agent produce work that appears complete but is not scientifically defensible.

## Judgment Standards

Use these labels consistently:

- ACCEPTABLE: Evidence-backed, reproducible within stated scope, explicit about limits, and consistent across README, state files, artifacts, notebook, and report.
- WEAK: Probably sound, but missing one important control, comparison, or reproducibility element.
- RISKY: Vulnerable to misleading conclusions because of validation gaps, ambiguous artifacts, stale state, or unclear assumptions.
- INCOMPLETE: Missing required outputs, validation, review state, or reproducibility pieces.
- MISLEADING: Claims or statuses overstate what the evidence supports, or repository metadata contradicts the study's real state.
- SYSTEMIC: The issue recurs across multiple studies or is embedded in core agent guidance.

When in doubt, prefer the highest severity that the evidence supports. Do not downgrade a systemic workflow issue just because one study happens to look strong.

## Behavioral Requirements

- Be critical but evidence-based.
- Avoid vague praise and avoid unsupported criticism.
- Point to concrete files, tables, figures, notebooks, or state files.
- Distinguish isolated defects from repo-wide guardrail failures.
- State uncertainty honestly when a conclusion depends on unverified assumptions.
- Prefer actionable recommendations over abstract commentary.
- Do not confuse existence of files with scientific completeness.
- Treat conflicting artifacts, stale state, or unchecked validation as high-risk in an autonomous-agent setting.

## Procedure

1. Define the audit set using recent study modification time plus recent review activity.
2. Build a per-study scorecard covering question, status, validation, reproducibility, traceability, and conclusion quality.
3. Inspect at least one report, one notebook, and one artifact set deeply enough to test whether the standards are real rather than nominal.
4. Compare repo-level instructions against actual skill enforcement and template coverage.
5. Separate findings into:
   - study-specific issues,
   - workflow gaps,
   - instruction and skill gaps,
   - and systemic autonomous-operation risks.
6. Prioritize fixes by how badly they could mislead a future agent or human.
7. Produce the audit report using the template in assets/AUDIT_REPORT_TEMPLATE.md.

## Expected Output

Every audit must produce:

- a recent-study scorecard,
- current strengths,
- current weaknesses and gaps,
- major risks for autonomous-agent scientific work,
- missing guardrails,
- recommended changes to AGENTS.md,
- recommended changes to existing skills,
- proposed process and repo-structure improvements,
- and prioritized action items.

## Common High-Risk Findings

- README says COMPLETE while study_state.json or task-run files show review still pending or revision requested.
- Validation checkboxes are unchecked or generic while conclusions read as final.
- Artifact files disagree about headline metrics.
- A reproducibility notebook exists but does not expose a clear load-vs-rerun path.
- Claims about optimality, impossibility, or robustness lack restarts, uncertainty, or comparison.
- Assumptions or approximation validity are stated but not numerically checked.
- AGENTS.md demands a rule that no skill, template, or validator actually enforces.

## Anti-Patterns

- Treating compiled PDFs as proof of scientific completeness.
- Rewarding file presence over reproducibility.
- Repeating README claims without cross-checking artifacts and task-run state.
- Calling a study complete because the loop produced a notebook and a report.
- Approving strong-sounding conclusions that are only exploratory.