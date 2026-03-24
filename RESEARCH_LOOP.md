# Continuous Research Loop — Architecture

> **System goal:** Click a study goal → autonomous two-model loop runs → produces validated, publication-quality results or a clear next-step plan with diagnostics.

---

## Two-Model Division of Labor

### Model 1: Science Director (Codex / GPT)

The **scientific brain**. Called only for high-value reasoning.

| Responsibility | When Called |
|---------------|------------|
| Understand cQED physics and cqed_sim deeply | Start of each review cycle |
| Propose research hypotheses | Planning and review phases |
| Design numerical experiments | Planning phase |
| Review execution outputs for physics correctness | After each implementation cycle |
| Judge result quality — is it actually good enough? | Review phase |
| Identify flaws, missing controls, stronger follow-ups | Review phase |
| Decide: continue, revise, branch, or stop | End of each review cycle |

**Input:** Compact structured state (study_state.json + key figure summaries + latest results digest)
**Output:** Structured directive (JSON/Markdown) for the Execution Agent

### Model 2: Execution Engineer (Opus 4.6)

The **research engineer + technical writer**. Handles all implementation.

| Responsibility | When Called |
|---------------|------------|
| Execute code changes, run simulations | Implementation phase |
| Generate figures and save data | Implementation phase |
| Debug failures (environment, syntax, runtime) | Self-debugging phase |
| Update documentation, reports, logs | After each task |
| Maintain reproducibility and folder structure | Continuously |
| Summarize results for Science Director review | End of implementation cycle |

**Input:** Structured directive from Science Director + full file access
**Output:** Updated study state + result summaries + figures

---

## Research Loop Protocol

```
┌──────────────────────────────────────────────────────────┐
│                    USER TRIGGERS STUDY                    │
│           (study goal + optional constraints)             │
└────────────────────────┬─────────────────────────────────┘
                         ▼
┌──────────────────────────────────────────────────────────┐
│              PHASE 0: BOOTSTRAP (Opus)                    │
│  • Create study folder structure                          │
│  • Initialize README, IMPROVEMENTS.md                     │
│  • Create task_runs/ state files                          │
│  • Write initial study_state.json                         │
│  • Lookup cqed_sim API for relevant modules               │
└────────────────────────┬─────────────────────────────────┘
                         ▼
┌──────────────────────────────────────────────────────────┐
│          PHASE 1: SCIENCE PLAN (Codex)                    │
│  • Read study goal + cqed_sim capabilities                │
│  • Classify problem (OPT/REP/DES/ANA)                    │
│  • Propose hypotheses and experiment design               │
│  • Define success criteria (quantitative)                 │
│  • Produce SCIENCE_DIRECTIVE.md                           │
│  • Estimate compute budget                                │
└────────────────────────┬─────────────────────────────────┘
                         ▼
┌──────────────────────────────────────────────────────────┐
│         PHASE 2: IMPLEMENT & EXECUTE (Opus)               │
│  • Read SCIENCE_DIRECTIVE.md                              │
│  • Write simulation scripts                               │
│  • Run simulations, save data                             │
│  • Generate figures                                       │
│  • Handle failures (self-debug, retry, log)               │
│  • Update study_state.json                                │
│  • Write EXECUTION_SUMMARY.md for Codex review            │
└────────────────────────┬─────────────────────────────────┘
                         ▼
┌──────────────────────────────────────────────────────────┐
│          PHASE 3: SCIENCE REVIEW (Codex)                  │
│  • Read EXECUTION_SUMMARY.md + key results                │
│  • Check physics correctness                              │
│  • Check whether results answer the real question          │
│  • Identify flaws, missing controls, better tests          │
│  • Decision:                                              │
│    ├─ CONTINUE → refine parameters, extend scope          │
│    ├─ REVISE   → change approach, new hypothesis          │
│    ├─ VALIDATE → results look good, run validation        │
│    └─ STOP     → blocked, needs human input               │
│  • Produce updated SCIENCE_DIRECTIVE.md                   │
└────────────────────────┬─────────────────────────────────┘
                         ▼
              ┌─── Decision ───┐
              │                │
        CONTINUE/REVISE    VALIDATE
              │                │
              ▼                ▼
         Back to          ┌────────────────────────┐
         Phase 2          │ PHASE 4: VALIDATE (Opus)│
                          │ • Sanity checks         │
                          │ • Convergence tests     │
                          │ • Literature comparison │
                          │ • Update study_state    │
                          └───────────┬────────────┘
                                      ▼
                          ┌────────────────────────┐
                          │ PHASE 5: REPORT (Opus)  │
                          │ • Write report.tex      │
                          │ • Compile PDF           │
                          │ • Finalize IMPROVEMENTS │
                          │ • Mark COMPLETE         │
                          └────────────────────────┘
```

---

## State Management

### study_state.json — Machine-Readable Study State

```json
{
  "study_name": "transmon_chi_optimization",
  "study_path": "studies/transmon_chi_optimization",
  "status": "IMPLEMENTING",
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
  "file_manifest": {
    "scripts": ["scripts/optimize_chi.py", "scripts/sweep_parameters.py"],
    "data": ["data/sweep_results.npz", "data/optimal_params.json"],
    "figures": ["figures/fidelity_vs_chi.png"],
    "report": "report/report.tex"
  }
}
```

### SCIENCE_DIRECTIVE.md — Codex → Opus Communication

Written by the Science Director after each review. Structured for machine parsing.

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

### EXECUTION_SUMMARY.md — Opus → Codex Communication

Written by the Execution Engineer after each implementation cycle. Compact.

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
  ├─ PHYSICS     → check Hamiltonian, units, parameter magnitudes
  └─ ASSUMPTION  → log as potential science issue for Codex review

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
| Opus gets structured directives, not prose | SCIENCE_DIRECTIVE.md uses action lists |
| State persists in files, not context | study_state.json is the single source of truth |
| Large files are never passed whole | Scripts summarized to function signatures + key results |
| Each iteration is self-contained | No reliance on chat history from prior iterations |
| Figures described, not embedded | Text description + file path, not base64 |
