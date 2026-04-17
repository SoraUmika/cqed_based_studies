# Continuous Research Loop -- Architecture

> **System goal:** Start the auto-watcher -> autonomous two-model loop runs fully automatically (no user intervention) -> produces validated, publication-quality results.

> **Canonical reference:** See `AGENTS.md` section "Multi-Agent Research Loop" for the full protocol specification. This document provides the architecture overview and state management details.

---

## Two-Model Division of Labor

### Model 1: Execution Engineer (Opus 4.6 -- Claude Code)

The **primary research executor and communicator**. Handles the full study lifecycle and report writing.

| Responsibility | When Called |
|---------------|------------|
| Self-generate research plan (SCIENCE_DIRECTIVE.md) | Start of each study / iteration |
| Classify problem (OPT/REP/DES/ANA) | Planning phase |
| Write simulation scripts using cqed_sim | Implementation phase |
| Run simulations, save data and figures | Implementation phase |
| Debug failures (environment, syntax, runtime) | Self-debugging phase |
| Validate results (sanity, convergence, literature) | Before reporting |
| Write LaTeX report and compile PDF | Report phase |
| Address reviewer feedback from FOLLOWUP_PROMPT.md | Revision iterations |
| Write REVIEW_REQUEST.md to signal reviewer | End of each iteration |

**Input:** User research prompt (iteration 1) or FOLLOWUP_PROMPT.md (revision iterations)
**Output:** Complete study output + REVIEW_REQUEST.md + EXECUTION_SUMMARY.md

### Model 2: Science Director / Critical Reviewer (Codex 5.4 xHigh -- GitHub Copilot)

The **independent reviewer and science director**. Evaluates quality and gates approval.

| Responsibility | When Called |
|---------------|------------|
| Read the full report (not just the summary) | After each REVIEW_REQUEST.md |
| Evaluate writing quality and readability | Review phase |
| Audit evidence-claim mapping for every non-trivial claim | Review phase |
| Check physics correctness and methodology | Review phase |
| Check completeness (reproducibility, artifacts, notebook) | Review phase |
| Decide: APPROVE, REVISE, or NEEDS_REWORK | End of review |
| Write REVIEW_DIRECTIVE.md with specific required actions | Review phase |
| Write FOLLOWUP_PROMPT.md for next Codex iteration | If REVISE / NEEDS_REWORK |
| Perform final readability polish | After APPROVE (Stage 4) |

**Input:** Full report.tex + EXECUTION_SUMMARY.md + figures + data
**Output:** REVIEW_DIRECTIVE.md + (if needed) FOLLOWUP_PROMPT.md + (after APPROVE) polished report

---

## Research Loop Protocol

```
STAGE 1: EXECUTE (Opus 4.6 -- Claude Code)
  PLAN -> IMPLEMENT -> VALIDATE -> REPORT -> REVIEW_REQUEST.md
         |
         v
STAGE 2: REVIEW (Codex 5.4 xHigh -- GitHub Copilot)
  Read full report -> Evaluate 4 dimensions -> REVIEW_DIRECTIVE.md
         |
         +--> APPROVE  --> STAGE 4: POLISH (Opus) --> COMPLETE
         |
         +--> REVISE / NEEDS_REWORK --> FOLLOWUP_PROMPT.md
                                              |
                                              v
                                   STAGE 3: REFINE (Opus)
                                     Read FOLLOWUP_PROMPT.md
                                     Address required actions
                                     Extend report
                                     REVIEW_REQUEST.md
                                              |
                                              v
                                   Back to STAGE 2 (review again)
```

### Stage Details

**Stage 1 -- Primary Research Execution (Opus)**
1. Read AGENTS.md, cqed_sim API Reference, existing study state
2. Write SCIENCE_DIRECTIVE.md (self-directed planning)
3. Implement scripts, run simulations, save data and figures
4. Validate (sanity checks, convergence, literature comparison)
5. Write report.tex with mandatory appendix, compile PDF
6. Write EXECUTION_SUMMARY.md and REVIEW_REQUEST.md
7. Signal reviewer

**Stage 2 -- Independent Critical Review (Codex)**

> Codex reviews as a **referee for a high-impact physics journal** (PRL / Nature Physics standard). Default stance is REVISE. APPROVE requires the reviewer to be able to say they would recommend this for publication without major revisions.

1. Read full report.tex (not just the execution summary) — summaries may omit weaknesses
2. Evaluate **five dimensions** (all must pass):
   - A. Writing quality and readability
   - B. Evidence-claim mapping (explicit audit table for every non-trivial claim)
   - C. Physics and methodology — including approximation validity bounds, convergence data, uncertainty quantification, multiple restarts for OPT, and sensitivity analysis
   - D. Completeness — reproducibility appendix, artifacts, notebook
   - E. Novelty and scientific significance — new insight vs. prior work, competitive metrics
3. Score each dimension 1–5 and assign an equivalent journal verdict
4. Write REVIEW_DIRECTIVE.md with decision and full journal score table
5. If REVISE or NEEDS_REWORK: write FOLLOWUP_PROMPT.md with specific, executable required actions

**Stage 3 -- Iterative Refinement (Opus)**
1. Read FOLLOWUP_PROMPT.md and REVIEW_DIRECTIVE.md
2. Address every required action from the review
3. Extend (not overwrite) report.tex
4. Write new EXECUTION_SUMMARY.md and REVIEW_REQUEST.md
5. Return to Stage 2

**Stage 4 -- Final Polish (Opus)**
1. After APPROVE: perform readability-only pass on full report
2. Improve sentence clarity, flow, transitions, caption quality
3. Do NOT re-evaluate physics (already approved)
4. Write POLISH_COMPLETE.md
5. Set study status to COMPLETE

### Decision Thresholds

| Decision | Journal equivalent | Criteria | Opus Action |
|----------|--------------------|----------|-------------|
| **APPROVE** | Accept / Minor cosmetic revision | All 5 dimensions pass; study is technically sound, adds verifiable new insight, and is well written; reviewer would recommend for high-impact journal publication | Proceed to Stage 4 polish |
| **REVISE** | Major / Minor revision | Core insight is present and physically sound, but evidence, methodology, scope, or presentation requires significant improvement | Address all required actions in REVIEW_DIRECTIVE.md; extend report |
| **NEEDS_REWORK** | Reject / Resubmit after redesign | Fundamental methodological flaw, unsupported central claims, no new insight over prior work, or any hard-reject trigger present | Revisit experimental design; re-run analyses; may require new approach |

**Hard-reject triggers** (any one → NEEDS_REWORK): unsupported central claim; approximation applied outside validity regime; optimization called "optimal" without restarts or global argument; convergence declared but not shown; no new insight over prior work.

---

## Automation

### Auto-Watcher (tools/auto_loop.ps1)

The watcher polls `task_runs/<study>/` for state-file signals and notifies the user when to paste prompts into GitHub Copilot Chat:

- **Watcher start** -> shows Stage 1 prompt, copies to clipboard
- **REVIEW_REQUEST.md appears** -> shows Stage 2 prompt, copies to clipboard
- **REVIEW_DIRECTIVE.md shows REVISE/NEEDS_REWORK** -> shows Stage 3 prompt, copies to clipboard
- **REVIEW_DIRECTIVE.md shows APPROVE** -> shows Stage 4 prompt, copies to clipboard

Start the watcher:
```powershell
# From VS Code: Terminal -> Run Task -> "Research: Start Auto-Watcher"
# Or from command line:
powershell -ExecutionPolicy Bypass -File tools\auto_loop.ps1 -StudyName <study_name>

# Dry-run (preview without calling CLI):
powershell -ExecutionPolicy Bypass -File tools\auto_loop.ps1 -StudyName <study_name> -DryRun
```

### What Is Automated vs. Manual

| Action | How |
|--------|-----|
| Opus execute (Stage 1) | Manual: Copy prompt to Copilot Chat, select Opus 4.6 model, paste |
| Codex review (Stage 2) | Manual: Copy prompt to Copilot Chat, select Codex 5.4 xHigh model, paste |
| Opus refine (Stage 3) | Manual: Copy prompt to Copilot Chat, select Opus 4.6 model, paste |
| Opus polish (Stage 4) | Manual: Copy prompt to Copilot Chat, select Opus 4.6 model, paste |

---

## State Management

### study_state.json -- Machine-Readable Study State

```json
{
  "study_name": "transmon_chi_optimization",
  "study_path": "studies/transmon_chi_optimization",
  "status": "REVIEW_REQUESTED",
  "problem_class": ["OPT"],
  "created_at": "2026-03-23T10:00:00Z",
  "updated_at": "2026-03-23T12:30:00Z",
  "loop_iteration": 2,
  "objective": "Optimize dispersive shift chi for high-fidelity readout",
  "hypotheses": [
    "Optimal chi is near -3 MHz for our transmon parameters"
  ],
  "assumptions": [
    "Dispersive regime valid (g << Delta)",
    "Cavity Kerr is small compared to chi"
  ],
  "success_criteria": {
    "primary": "Readout fidelity > 99.5%",
    "secondary": "Chi optimized within 5% of analytic prediction"
  },
  "completed_tasks": ["P1.1", "P1.2", "P2.1"],
  "failed_tasks": [
    {
      "id": "P2.2",
      "reason": "Nelder-Mead stalled at 98.7% fidelity",
      "attempted_fixes": ["Increased iterations to 1000", "Tried different initial simplex"],
      "resolution": "Switching to GRAPE optimizer"
    }
  ],
  "pending_tasks": ["P2.3", "P3.1", "P3.2"],
  "blocked_tasks": [],
  "key_results": {
    "best_fidelity": 0.987,
    "optimal_chi_mhz": -2.91,
    "convergence_iterations": 500
  },
  "latest_figures": [
    "figures/fidelity_vs_chi.png",
    "figures/convergence_trace.png"
  ],
  "blockers": [],
  "compute_notes": {
    "total_wall_time_min": 45,
    "bottleneck": "Parameter sweep over 50x50 grid"
  },
  "science_directive_version": 2,
  "review_iterations": 1,
  "reviewer_decision": "REVISE",
  "file_manifest": {
    "scripts": ["scripts/optimize_chi.py", "scripts/sweep_parameters.py"],
    "data": ["data/sweep_results.npz", "data/optimal_params.json"],
    "figures": ["figures/fidelity_vs_chi.png"],
    "report": "report/report.tex"
  }
}
```

### State File Communication Flow

| File | Written By | Read By | Purpose |
|------|-----------|---------|---------|
| `SCIENCE_DIRECTIVE.md` | Opus (self-directed) | Opus (during implement) | Research plan and experiment design |
| `EXECUTION_SUMMARY.md` | Opus | Codex (during review) | Quantitative findings digest |
| `REVIEW_REQUEST.md` | Opus | Codex / auto-watcher | Signal that study is ready for review |
| `REVIEW_DIRECTIVE.md` | Codex | Opus (during revision) | Review assessment and required actions |
| `FOLLOWUP_PROMPT.md` | Codex | Opus (during revision) | Self-contained prompt for next iteration |
| `POLISH_COMPLETE.md` | Opus | User / auto-watcher | Signal that study is finalized |
| `study_state.json` | Both | Both | Machine-readable status (single source of truth) |
| `TASK_CHECKLIST.md` | Opus | Opus (on resume) | Checkpointed task tracking |
| `PROGRESS_LOG.md` | Opus | Both | Append-only log of what happened |
| `BLOCKERS.md` | Opus | Both | Active and resolved blockers |

### SCIENCE_DIRECTIVE.md -- Self-Directed Research Plan (Opus)

```markdown
# Science Directive — Iteration N

## Decision
<!-- CONTINUE | REVISE | VALIDATE | STOP -->

## Assessment of Previous Results
<!-- What was good, what was wrong, what's missing -->

## Next Actions (ordered)
1. **[ACTION_TYPE]** Description
   - Files to create/modify: ...
   - Expected output: ...
   - Success criterion: ...

## Hypotheses Update
<!-- Any new or revised hypotheses -->

## Open Concerns
<!-- Physics issues, numerical concerns, things to watch -->

## Stopping Criteria for This Iteration
<!-- When should Opus stop and send results back for review? -->
```

### EXECUTION_SUMMARY.md -- Opus -> Codex Communication

Written by the Execution Engineer (Opus) after each implementation cycle. Read by the Science Director (Codex).

```markdown
# Execution Summary — Iteration N

## Tasks Completed
- [x] P2.1: Implemented chi sweep script
- [x] P2.2: Generated fidelity heatmap

## Tasks Failed
- P2.3: GRAPE optimizer — ImportError on qutip.control
  - Attempted: pip install --user qutip, checked API
  - Resolution: Using scipy.optimize.minimize with L-BFGS-B instead

## Key Results
- Best fidelity: 0.987 (at chi = -2.91 MHz)
- Convergence: stable after 350 iterations
- [See figure: figures/fidelity_vs_chi.png]

## Result Digest (for Codex review)
<!-- 5-10 bullet points summarizing what was found -->

## Anomalies / Concerns
<!-- Anything unexpected that Codex should evaluate -->

## Updated File Manifest
<!-- Which files were created or changed -->

## Compute Notes
- Wall time: 23 min for parameter sweep
- Memory: ~2 GB peak
```

---

## Self-Debugging Protocol

When a task fails, Opus follows this escalation ladder:

```
Level 1: INSPECT (< 30 seconds)
  ├─ Read error message / traceback
  ├─ Check: syntax error? import error? file not found?
  └─ Classify: ENVIRONMENT | DEPENDENCY | SYNTAX | RUNTIME | PHYSICS | ASSUMPTION

Level 2: FIX (bounded: max 3 attempts per failure)
  ├─ ENVIRONMENT → check paths, permissions, Python version
  ├─ DEPENDENCY  → pip install --user, check installed version
  ├─ SYNTAX      → fix the code, re-run
  ├─ RUNTIME     → check array shapes, parameter ranges, NaN/Inf
  \-- PHYSICS \-> check Hamiltonian, units, parameter magnitudes
  \-- ASSUMPTION \-> log as potential issue for Codex review

Level 3: LOG & ESCALATE (after 3 failed attempts)
  ├─ Document: what was tried, what happened, full traceback
  ├─ Add to BLOCKERS.md with category tag
  ├─ Add to study_state.json failed_tasks
  ├─ Continue with next non-blocked task
  └─ Flag for Codex review in EXECUTION_SUMMARY.md

Level 4: STOP (only if ALL remaining tasks are blocked)
  ├─ Write comprehensive blocker report
  ├─ Save all partial results
  └─ Update study_state.json status = "BLOCKED"
```

---

## Research-Quality Stopping Criteria

The loop does NOT stop merely because code ran. Before marking COMPLETE, ALL must pass:

```
□ Scientific question answered?
  └─ Results directly address the study objective
  └─ Key figures exist that demonstrate the answer

□ Physics consistency?
  └─ Limiting cases match known results
  └─ Conservation laws satisfied
  └─ Parameter magnitudes physically reasonable

□ Diagnostics and controls?
  └─ Convergence verified (Hilbert space, time steps, iterations)
  └─ At least one sanity check documented
  └─ Literature comparison if applicable

□ Robustness?
  └─ Results stable to small parameter perturbations
  └─ Not obviously stuck in local minimum (if optimization)

□ Documentation complete?
  └─ report.tex written with mandatory appendix
  └─ IMPROVEMENTS.md finalized
  └─ study_state.json status = COMPLETE
  └─ All figures saved in .png and .pdf

□ Open questions documented?
  └─ IMPROVEMENTS.md has Open Questions section
  └─ report.tex has Limitations and Future Work section
```

---

## Token Efficiency Rules

| Principle | Implementation |
|-----------|---------------|
| Codex sees summaries, not raw data | EXECUTION_SUMMARY.md is capped at ~500 lines |
| Opus gets structured directives, not prose | FOLLOWUP_PROMPT.md uses action lists |
| State persists in files, not context | study_state.json is the single source of truth |
| Large files are never passed whole | Scripts summarized to function signatures + key results |
| Each iteration is self-contained | No reliance on chat history from prior iterations |
| Figures described, not embedded | Text description + file path, not base64 |
