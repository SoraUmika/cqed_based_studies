# Agentic-Assisted Research Workflow Analysis and Reproducible Setup

**cqed_sim + VS Code + Multi-Agent Models**

> A comprehensive technical report enabling external researchers to understand,
> reproduce, and adapt the autonomous agent-assisted research workflow developed
> for circuit quantum electrodynamics (cQED) simulation studies.

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Part 1 — cqed_sim: Source-of-Truth Simulator Analysis](#part-1--cqed_sim-source-of-truth-simulator-analysis)
   - [1.1 Architecture](#11-architecture)
   - [1.2 Physics & Conventions](#12-physics--conventions)
   - [1.3 Simulation Engine](#13-simulation-engine)
   - [1.4 Correctness & Validation](#14-correctness--validation)
   - [1.5 Agent Integration](#15-agent-integration)
3. [Part 2 — cqed_based_study: Agentic Research Layer](#part-2--cqed_based_study-agentic-research-layer)
   - [2.1 Study Structure](#21-study-structure)
   - [2.2 Execution Pipeline](#22-execution-pipeline)
   - [2.3 Agent Roles](#23-agent-roles)
   - [2.4 Outputs](#24-outputs)
   - [2.5 Failure Handling](#25-failure-handling)
4. [Part 3 — Agent Workflow Design](#part-3--agent-workflow-design)
   - [3.1 Multi-Agent System Design](#31-multi-agent-system-design)
   - [3.2 Iterative Loop](#32-iterative-loop)
   - [3.3 Prompting Strategy](#33-prompting-strategy)
   - [3.4 AGENTS.md System](#34-agentsmd-system)
5. [Part 4 — VS Code Setup (Reproducible Guide)](#part-4--vs-code-setup-reproducible-guide)
   - [4.1 Required Tools](#41-required-tools)
   - [4.2 Model Usage Strategy](#42-model-usage-strategy)
   - [4.3 Recommended Workflow in VS Code](#43-recommended-workflow-in-vs-code)
   - [4.4 Automated Loop (Advanced)](#44-automated-loop-advanced)
   - [4.5 Terminal / Execution Integration](#45-terminal--execution-integration)
6. [Part 5 — Generalization to Other Fields](#part-5--generalization-to-other-fields)
7. [Part 6 — Recommendations & Improvements](#part-6--recommendations--improvements)
8. [Appendix A — File Tree Reference](#appendix-a--file-tree-reference)
9. [Appendix B — Complete Agent Definition Files](#appendix-b--complete-agent-definition-files)
10. [Appendix C — Glossary](#appendix-c--glossary)

---

## 1. Executive Summary

This repository implements a **fully autonomous, agent-assisted research workflow** for computational physics studies in circuit quantum electrodynamics (cQED). The system combines:

- **`cqed_sim`** — a hardware-faithful, pulse-level cQED simulator built on QuTiP
- **`cqed_based_study`** — a structured research workspace with agent definitions, skills, tools, and study templates
- **Multi-model AI agents** — coordinated through VS Code's GitHub Copilot infrastructure using Claude Opus 4.6, OpenAI Codex, and Copilot Chat

The workflow operates as a **continuous research loop**: a user specifies a research question, and the system autonomously plans experiments, writes simulation code, runs simulations, validates results, and produces publication-quality LaTeX reports — all with structured state management that survives interruptions and enables resumable execution.

**Six completed studies** demonstrate the workflow's capabilities:
1. Dispersive readout pulse optimization (4 progressive sub-studies)
2. Gray-box adaptive control for cQED systems
3. Thermal noise cavity sensing
4. SQR gate design
5. Hybrid qubit-cavity control
6. Literature-informed selective pulse primitives

All studies have been re-run and validated end-to-end, producing 69/69 passing validation checks across the test suite.

---

## Part 1 — cqed_sim: Source-of-Truth Simulator Analysis

### 1.1 Architecture

#### Design Philosophy

`cqed_sim` is a **hardware-faithful, time-domain, pulse-level circuit-QED simulator**. It models the full dispersive Hamiltonian with explicit drive waveforms, hardware distortion effects, and Lindblad open-system dynamics. The design prioritizes:

1. **Physical fidelity** — models include higher-order dispersive terms (falling-factorial expansion), multi-mode coupling, and hardware signal-chain effects
2. **Composability** — frozen dataclass models, cached operators, and builder patterns allow mixing and matching physics components
3. **Reproducibility** — deterministic solver backends, unit-coherent conventions, and JSON-serializable gate sequences
4. **Extensibility** — plugin backend system (QuTiP, NumPy-dense, JAX-dense), modular subpackages

#### Module Structure

The framework is organized into clearly separated layers:

```
cqed_sim/
├── core/           ← Physics models, Hilbert space, ideal states/gates
│   ├── models.py          # DispersiveTransmonCavityModel, UniversalCQEDModel
│   ├── frames.py          # FrameSpec, rotating-frame definitions
│   ├── ideal_gates.py     # Perfect unitary gate library
│   ├── state_prep.py      # Initial state construction
│   └── energy_spectrum.py # Dressed energy-level solvers
│
├── pulses/         ← Pulse construction and calibration
│   ├── pulse.py           # Pulse dataclass (frozen)
│   ├── envelope.py        # GaussianEnvelope, FlatTopEnvelope, etc.
│   ├── builders.py        # build_displacement_pulse, build_rotation_pulse, ...
│   ├── calibration.py     # Analytic pulse calibration formulas
│   └── hardware.py        # HardwareConfig (IQ distortion, ZOH, DAC)
│
├── sequence/       ← Waveform compilation pipeline
│   └── compiler.py        # SequenceCompiler: pulse list → sampled timeline
│
├── sim/            ← Simulation engine
│   ├── engine.py          # simulate_sequence(), SimulationSession
│   ├── hamiltonian.py     # Hamiltonian assembly from model + drives
│   ├── noise.py           # NoiseSpec, Lindblad collapse operators
│   ├── result.py          # SimulationResult with state extractors
│   └── diagnostics.py     # Runtime checks and convergence monitors
│
├── backends/       ← Pluggable solver backends
│   ├── base.py            # BaseBackend ABC
│   ├── numpy_backend.py   # Dense matrix expm propagation
│   └── jax_backend.py     # JAX-accelerated dense backend
│
├── measurement/    ← Readout modeling
│   ├── qubit.py           # Ideal, confusion-matrix, and IQ measurement
│   ├── readout_chain.py   # Resonator + Purcell filter + amplifier
│   └── continuous.py      # SME-based continuous monitoring
│
├── floquet/        ← Periodic-drive Floquet analysis
├── analysis/       ← Parameter translation (bare → dispersive)
├── calibration/    ← SQR gate calibration, multitone validation
├── calibration_targets/  ← Surrogate experiments (spectroscopy, Rabi, Ramsey, ...)
├── gates/          ← Ideal unitary gate library (100+ gates)
├── operators/      ← Cached Pauli, ladder, Fock projector operators
├── observables/    ← Bloch, Fock-resolved, phase, Wigner diagnostics
├── plotting/       ← Publication-quality visualization
├── tomo/           ← State/process tomography, leakage matrix calibration
├── unitary_synthesis/  ← Gate-sequence optimal control
├── optimal_control/    ← Direct GRAPE with hardware-aware forward model
├── rl_control/     ← Reinforcement learning environment
├── system_id/      ← Calibration-informed priors
├── quantum_algorithms/ ← Holographic algorithm utilities
└── io/             ← JSON I/O for gate sequences
```

#### Separation of Concerns

| Layer | Responsibility | Key Classes |
|-------|---------------|-------------|
| **Physics Model** | Hamiltonian parameters, coupling specs, rotating frames | `UniversalCQEDModel`, `DispersiveTransmonCavityModel`, `FrameSpec` |
| **Pulse Construction** | Waveform generation, envelope shaping, calibration | `Pulse`, `PulseBuilder`, `HardwareConfig` |
| **Numerical Solver** | Time evolution, state propagation, noise channels | `simulate_sequence()`, `SimulationSession`, `NoiseSpec` |
| **Analysis / API** | Results extraction, plotting, tomography, optimization | `SimulationResult`, `UnitarySynthesizer`, `GrapeSolver` |

### 1.2 Physics & Conventions

#### Hamiltonian Construction

The dispersive Hamiltonian for a 2-mode (qubit + cavity) system in the rotating frame:

$$H_0/\hbar = \delta_c \hat{n}_c + \delta_q \hat{n}_q + \frac{\alpha}{2} \hat{b}^{\dagger 2}\hat{b}^2 + \frac{K}{2}\hat{n}_c(\hat{n}_c - 1) + \chi \hat{n}_c \hat{n}_q + \chi_2 \hat{n}_c(\hat{n}_c - 1)\hat{n}_q + \cdots$$

where:
- $\delta_c, \delta_q$ = detunings from the rotating frame
- $\alpha$ = transmon anharmonicity (negative, typically −255 MHz)
- $K$ = cavity self-Kerr (typically −28 kHz)
- $\chi$ = dispersive shift (−2.84 MHz in the reference parameter set)
- $\chi_2$ = second-order dispersive shift

**Higher-order terms** use the **falling-factorial** form: $\chi_{\text{higher}}[i]$ multiplies $\hat{n}(\hat{n}-1)\cdots(\hat{n}-i)$, not $\hat{n}^{i+1}$.

#### Sign Conventions

| Quantity | Convention | Physical Meaning |
|----------|-----------|------------------|
| $\chi$ | Positive $\chi$ → qubit frequency **increases** with photon number | $\omega_q(n) = \omega_q + \chi \cdot n$ |
| $\alpha$ | Negative for transmon | $\omega_{ef} = \omega_{ge} + \alpha$ |
| $K$ | Negative (cavity self-Kerr) | Photon-number-dependent frequency pull |
| Drive carrier | $\omega_{\text{carrier}} = -\omega_{\text{transition}}$ | Uses $e^{+i\omega t}$ convention |

#### Drive Hamiltonian

$$H_{\text{drive}} = \epsilon(t) \hat{O}^+ + \epsilon^*(t) \hat{O}^-$$

where $\hat{O}^{\pm}$ are raising/lowering operators resolved from string targets (`"qubit"`, `"storage"`, `"sideband"`) or structured specs (`TransmonTransitionDriveSpec`, `SidebandDriveSpec`).

#### Rotating Frames

`FrameSpec` defines the rotating frame for each mode. The default frame rotates at the mode frequencies, placing the dispersive shifts as the dominant energy scales. All pulse carriers are defined relative to the frame.

#### Alignment with Physics Conventions Document

The package includes `physics_and_conventions/conventions.py` with enum types (`UnitType`, `DetuningSign`, `TensorOrdering`) and enforcement decorators that ensure sign and ordering consistency across all modules. The physics conventions document (`physics_conventions_report.tex`) provides the full derivation chain from the Jaynes-Cummings model through the dispersive limit.

### 1.3 Simulation Engine

#### Solver Backend

**Primary path (QuTiP):**
- Pure states: `qutip.sesolve` (Schrödinger equation)
- Mixed states / open systems: `qutip.mesolve` (Lindblad master equation)
- Steady state: `qutip.steadystate` (for equilibrium calculations)

**Dense backend path (NumPy / JAX):**
- Piecewise-constant matrix exponential propagation (`scipy.linalg.expm` or `jax.scipy.linalg.expm`)
- For small systems and parity checks
- JAX backend enables GPU acceleration for GRAPE optimal control

#### Time Evolution Pipeline

```
┌──────────────┐     ┌──────────────┐     ┌──────────────────┐     ┌───────────────┐
│  Physics     │     │  Pulse       │     │  Sequence        │     │  Simulation   │
│  Model       │────▶│  Objects     │────▶│  Compiler        │────▶│  Engine       │
│              │     │              │     │                  │     │               │
│ • Parameters │     │ • Envelope   │     │ • Sample pulses  │     │ • Build H(t)  │
│ • Coupling   │     │ • Carrier    │     │ • Apply hardware │     │ • Solve ODE   │
│ • Frame      │     │ • Phase      │     │ • ZOH, lowpass,  │     │ • Extract     │
│              │     │ • DRAG       │     │   quantization   │     │   observables │
└──────────────┘     └──────────────┘     └──────────────────┘     └───────────────┘
                                                                          │
                                                                          ▼
                                                                   ┌───────────────┐
                                                                   │  Simulation   │
                                                                   │  Result       │
                                                                   │               │
                                                                   │ • States(t)   │
                                                                   │ • Populations │
                                                                   │ • Bloch coords│
                                                                   │ • Wigner      │
                                                                   │ • Fidelity    │
                                                                   └───────────────┘
```

**Detail of each stage:**

1. **Model construction**: Frozen dataclass with all Hamiltonian parameters. Three model classes wrap `UniversalCQEDModel`:

   | Model | Modes | Use Case |
   |-------|-------|----------|
   | `DispersiveTransmonCavityModel` | qubit + cavity | Single-mode readout or control |
   | `DispersiveReadoutTransmonStorageModel` | qubit + storage + readout | Full 3-mode experiments |
   | `UniversalCQEDModel` | Arbitrary N modes | General multi-mode simulation |

2. **Pulse construction**: `Pulse` dataclass with envelope (analytic or pre-sampled), carrier frequency, phase, amplitude, DRAG coefficient, and sample rate. Waveform: $\epsilon(t) = \text{amp} \cdot \text{env}(t_{\text{rel}}) \cdot e^{i(\text{carrier} \cdot t + \text{phase})}$

3. **Sequence compilation**: `SequenceCompiler.compile()` takes a list of `Pulse` objects and produces sampled waveform timelines, applying hardware effects: zero-order hold (ZOH), lowpass filtering, DAC quantization, timing quantization, IQ distortion, and crosstalk mixing.

4. **Simulation**: `simulate_sequence()` assembles the time-dependent Hamiltonian $H_0 + H_{\text{drive}}(t)$, adds Lindblad collapse operators from `NoiseSpec`, and calls the solver.

5. **Result extraction**: `SimulationResult` provides partial traces, Bloch coordinates, Fock populations, Wigner functions, fidelity calculations, and more.

#### High-Throughput Execution

- `SimulationSession` / `prepare_simulation()` pre-computes the Hamiltonian once
- `run_many()` / `simulate_batch()` provides parallel execution via `ProcessPoolExecutor`
- Session reuse pattern eliminates redundant operator construction across parameter sweeps

#### Noise Modeling

`NoiseSpec` supports:

| Channel | Parameters | Collapse Operator |
|---------|-----------|-------------------|
| Transmon $T_1$ | Aggregate or per-ladder-transition | $\sqrt{\gamma_1} \hat{b}$ |
| Transmon pure dephasing | $T_\phi$ (rate $\gamma_\phi = 1/2T_\phi$) | $\sqrt{\gamma_\phi} \hat{\sigma}_z$ |
| Cavity loss | $\kappa$, $n_{\text{th}}$ | $\sqrt{\kappa(n_{\text{th}}+1)} \hat{a}$, $\sqrt{\kappa n_{\text{th}}} \hat{a}^\dagger$ |
| Cavity dephasing | $T_\phi$ (rate $\gamma_\phi = 1/T_\phi$) | $\sqrt{\gamma_\phi} \hat{n}$ |

Supports `split_collapse_operators()` for separating monitored vs. unmonitored channels in stochastic master equation (SME) simulations.

### 1.4 Correctness & Validation

#### Test Infrastructure

The framework includes **58+ test files** organized by physics topic:

| Test Category | Files | What Is Tested |
|---------------|-------|----------------|
| Free evolution & sanity | `test_01` | No-drive identity, energy conservation |
| Cavity drive & Kerr | `test_02` | Displacement fidelity, Kerr rotation |
| Dispersive & Ramsey | `test_03` | χ-dependent phase accumulation |
| χ convention | `test_10` | Sign convention consistency |
| Model invariants | `test_11` | Parameter round-trip, frozen dataclass |
| Pulse semantics | `test_12` | Carrier/phase/envelope composition |
| Dissipation | `test_14` | T1 decay rates, steady-state photon numbers |
| Convergence regression | `test_07` | Step-size convergence, truncation stability |
| Hardware effects | `test_08` | Timeline compilation, ZOH, quantization |
| Gate library | `test_17`, `test_42` | Unitary gate correctness for 100+ gates |
| Three-mode model | `test_27`, `test_28` | Multi-mode Hamiltonian assembly |
| Universal model | `test_34` | Arbitrary coupling specifications |
| GRAPE optimal control | `test_40`, `test_41`, `test_51`, `test_52` | Convergence, gradient accuracy |
| Floquet analysis | `test_58` | Quasienergies, multiphoton resonances |
| API completeness | `test_37` | All public symbols exported correctly |

Additional specialized test subdirectories: `tests/analysis/`, `tests/calibration_targets/`, `tests/conventions/`, `tests/experiment/`, `tests/golden/`, `tests/quantum_algorithms/`, `tests/rl_control/`, `tests/sim/`, `tests/unitary_synthesis/`

#### Verification Methods

1. **Analytic limiting cases**: Zero-drive → identity evolution; weak-coupling → Jaynes-Cummings; harmonic limit → exact coherent state
2. **Cross-backend parity**: QuTiP solver results compared against dense NumPy/JAX backends
3. **Golden-file regression**: Key simulation outputs stored as reference data; CI checks for drift
4. **Convergence sweeps**: Hilbert space truncation, time-step refinement, and optimizer iteration counts
5. **Convention reconciliation**: Dedicated `tests/conventions/` ensuring all sign definitions are self-consistent

### 1.5 Agent Integration

Agents interact with `cqed_sim` through a structured skill system:

#### How Agents Write Simulation Code

1. **API lookup first**: The `cqed-sim-lookup` skill requires agents to read the API Reference before writing any simulation code. This ensures agents use existing functionality rather than duplicating it.

2. **Coverage assessment**: Agents classify required features as:
   - **Fully supported** → use `cqed_sim` directly
   - **Partially supported** → extend using `cqed_sim` as foundation, document gap
   - **Not supported** → write standalone code, document gap, suggest upstreaming

3. **Code generation**: Agents write Python scripts that import `cqed_sim` classes and follow established patterns:
   ```python
   from cqed_sim import DispersiveTransmonCavityModel, simulate_sequence, NoiseSpec
   
   model = DispersiveTransmonCavityModel(
       omega_c=5.241e9 * 2 * np.pi,
       omega_q=6.150e9 * 2 * np.pi,
       chi=-2.84e6 * 2 * np.pi,
       ...
   )
   result = simulate_sequence(model, pulses, config)
   ```

#### Iterative Refinement Loop

```
Agent writes script
       │
       ▼
Agent runs script  ──── Error? ──── Agent reads traceback
       │                                    │
       │                               Classifies error
       │                               (ENVIRONMENT / DEPENDENCY /
       │                                SYNTAX / RUNTIME / PHYSICS)
       │                                    │
       │                               Applies fix (max 3 attempts)
       │                                    │
       ▼                                    ▼
Results produced ◄──────────────── Re-runs script
       │
       ▼
Agent validates results
  (sanity checks, convergence, literature comparison)
       │
       ▼
Results pass? ── No ──▶ Log in IMPROVEMENTS.md, iterate
       │
      Yes
       │
       ▼
Generate figures + report
```

---

## Part 2 — cqed_based_study: Agentic Research Layer

### 2.1 Study Structure

Every study follows a standardized directory layout:

```
studies/<study_name>/
├── README.md           ← Study definition (goals, methods, status)
├── IMPROVEMENTS.md     ← Living log of limitations, ideas, failed approaches
├── study_state.json    ← Machine-readable state (for agent coordination)
├── scripts/            ← Python simulation and analysis code
├── data/               ← Raw and processed numerical outputs (.npz, .json)
├── figures/            ← Plots in both .png (300 dpi) and .pdf (vector)
└── report/
    ├── report.tex      ← LaTeX report following AGENTS.md template
    ├── references.bib  ← BibTeX bibliography
    └── report.pdf      ← Compiled PDF
```

#### README.md — Study Contract

Every study README contains mandatory sections that serve as the contract between the human researcher and the agents:

| Section | Purpose |
|---------|---------|
| **Problem Class** | OPT / REP / DES / ANA — determines which validation checks and appendix content apply |
| **Motivation** | Why this study matters; links to papers for REP-class |
| **Goals** | Numbered, concrete, falsifiable objectives |
| **Methods** | Which `cqed_sim` modules will be used; documented gaps |
| **Expected Outcomes** | Quantitative success criteria |
| **Known Limitations** | Updated throughout; feeds into report |
| **Status** | ACTIVE / COMPLETE / BLOCKED |

#### IMPROVEMENTS.md — The Living Improvement Log

This file is the bridge between current and future work. Agents update it **in real time** during implementation, not just at the end. Structure:

- **Critical Gaps (P1)** — things that could make results qualitatively wrong
- **Recommended Improvements (P2)** — meaningful accuracy or scope improvements
- **Nice-to-Haves (P3)** — lower-priority enhancements
- **Open Questions** — unresolved physics observations
- **What Was Tried and Did Not Work** — documents dead ends to prevent repetition
- **Compute & Resource Notes** — wall-clock times, memory usage, bottlenecks

### 2.2 Execution Pipeline

#### Standard Workflow

```
Step 1: Initialize  →  Create study folder, README, IMPROVEMENTS.md
Step 2: Plan        →  Check cqed_sim API, identify gaps, state assumptions
Step 3: Implement   →  Write scripts, run simulations, save data, generate figures
Step 4: Validate    →  Sanity checks ✓  Convergence ✓  Literature comparison ✓
Step 5: Report      →  Write report.tex with mandatory appendices → compile PDF
```

#### Parameter Sweeps and Optimization

Studies call into `cqed_sim` through several patterns:

| Pattern | Use Case | Example |
|---------|----------|---------|
| **Single simulation** | Verify a specific parameter point | `simulate_sequence(model, pulses, config)` |
| **Parameter sweep** | Map out fidelity landscape | Loop over $\chi$, $\kappa$, amplitude grids |
| **Batch execution** | High-throughput parallel sweeps | `SimulationSession.run_many()` |
| **Optimization** | Find optimal control pulse | `GrapeSolver`, `UnitarySynthesizer`, `scipy.optimize` |
| **Calibration** | Targeted parameter extraction | `cqed_sim.calibration_targets` surrogates |

#### Example: Gray-Box Adaptive Control Study

This study demonstrates the full pipeline:

1. **Phase 4**: Compare nominal, gray-box, perfect-knowledge, and black-box control across chi mismatch levels (0–40%)
2. **Phase 5**: Stress-test robustness under noise, readout confusion, probe budget reduction, drift, and Hamiltonian omissions
3. **Validation**: 12/12 checks pass; multistart GRAPE with 3 seeds; convergence verified under expanded Hilbert space
4. **Output**: 9 figures, 7 data files, validated LaTeX report

### 2.3 Agent Roles

The system uses a **two-model architecture** with clearly separated responsibilities:

#### Science Director (Codex / GPT)

| Capability | Description |
|------------|-------------|
| **Role** | The scientific brain — reasons about physics, not implementation |
| **When called** | Planning phase (design experiments) and review phase (evaluate results) |
| **Input** | Compact structured state (study_state.json + figure summaries + results digest) |
| **Output** | Structured SCIENCE_DIRECTIVE.md with ordered action items |
| **Strengths** | Physics reasoning, hypothesis generation, experiment design, quality judgment |
| **Does NOT** | Write code, run simulations, debug errors |

#### Execution Engineer (Claude Opus 4.6)

| Capability | Description |
|------------|-------------|
| **Role** | Research engineer + technical writer — handles all implementation |
| **When called** | Bootstrap, implement, validate, and report phases |
| **Input** | SCIENCE_DIRECTIVE.md + full file access |
| **Output** | Code, data, figures, documentation, EXECUTION_SUMMARY.md |
| **Strengths** | Code generation, debugging, documentation, LaTeX report writing, structured reasoning |
| **Does NOT** | Make physics judgments or decide research direction |

#### Copilot Chat (Codex-medium)

| Capability | Description |
|------------|-------------|
| **Role** | Fast-iteration assistant for quick edits and experiments |
| **When called** | Interactive editing, quick debugging, ad-hoc queries |
| **Strengths** | Speed, low latency, good for routine code changes |
| **Limitation** | Less reliable for deep physics reasoning or complex multi-step tasks |

### 2.4 Outputs

Each completed study produces:

| Artifact | Format | Purpose |
|----------|--------|---------|
| **Data files** | `.npz`, `.json` | Reproducible numerical results |
| **Figures** | `.png` (300 dpi) + `.pdf` (vector) | Publication-quality plots with colorblind-friendly palettes |
| **LaTeX report** | `report.tex` + `report.pdf` | Peer-review-style document with abstract, methods, results, validation, discussion, limitations, appendix |
| **IMPROVEMENTS.md** | Markdown | Actionable handoff document for future agents/researchers |
| **study_state.json** | JSON | Machine-readable state for agent coordination |

#### Report Structure (Mandatory)

```
Abstract → Introduction → System & Methods → Results → Validation
→ Discussion → Conclusion → Limitations & Future Work → References → Appendices
```

The appendix is **required** and contains detailed data (pulse shapes, full parameter tables, sweep data, cost landscapes) that supports the main text's findings and interpretation.

### 2.5 Failure Handling

#### Self-Debugging Protocol (4-Level Escalation)

```
Level 1: INSPECT (< 30 seconds)
  ├─ Read error traceback
  ├─ Classify: ENVIRONMENT | DEPENDENCY | SYNTAX | RUNTIME | PHYSICS | ASSUMPTION
  └─ Apply targeted fix

Level 2: FIX (max 3 attempts per failure)
  ├─ ENVIRONMENT → check paths, permissions, Python version
  ├─ DEPENDENCY  → pip install --user, check version
  ├─ SYNTAX      → fix code, re-run
  ├─ RUNTIME     → check shapes, parameter ranges, NaN/Inf
  ├─ PHYSICS     → check Hamiltonian, units, magnitudes
  └─ ASSUMPTION  → log as science issue for director review

Level 3: LOG & ESCALATE (after 3 failed attempts)
  ├─ Document in BLOCKERS.md with full traceback
  ├─ Add to study_state.json failed_tasks
  ├─ Continue with next non-blocked task
  └─ Flag in EXECUTION_SUMMARY.md for review

Level 4: STOP (only if ALL remaining tasks blocked)
  ├─ Write comprehensive blocker report
  ├─ Save all partial results
  └─ Set status = "BLOCKED"
```

#### Recovery from Interruptions

The file-based state system enables seamless recovery:
- `study_state.json` tracks completed, failed, pending, and blocked tasks
- `TASK_CHECKLIST.md` serves as the execution source of truth
- `PROGRESS_LOG.md` provides append-only checkpoints
- The `autonomous-resume` agent reconstructs state and continues from exactly where work stopped

---

## Part 3 — Agent Workflow Design

### 3.1 Multi-Agent System Design

#### Role Differentiation

| Role | Model | Responsibility | Invocation Method |
|------|-------|---------------|-------------------|
| **Science Director** | OpenAI Codex / GPT | Physics reasoning, experiment design, result review, quality judgment | `@science-director study=... phase=plan\|review` |
| **Execution Engineer** | Claude Opus 4.6 | Code writing, simulation execution, debugging, documentation, reports | `@execution-engineer study=... phase=bootstrap\|implement\|validate\|report` |
| **Research Loop Orchestrator** | Combined (switches hats) | End-to-end autonomous coordination | `@research-loop study=... goal='...'` |
| **Autonomous Planner** | General | Convert task documents into execution plans | `/Autonomous Plan task=... run=...` |
| **Autonomous Implementer** | General | Execute checklist items with state tracking | `/Autonomous Implement task=... run=...` |
| **Autonomous Resume** | General | Recover interrupted tasks from file state | `/Autonomous Resume task=... run=...` |

#### Why Two Models?

The two-model split is driven by complementary strengths:

| Dimension | Science Director (Codex/GPT) | Execution Engineer (Opus) |
|-----------|------------------------------|--------------------------|
| **Physics reasoning** | Deep domain knowledge, hypothesis generation | Follows directives, does not second-guess physics |
| **Code generation** | Limited to pseudocode / API references | Full implementation with debugging |
| **Tool access** | Read-only (search, read files) | Full access (read, write, execute, edit) |
| **Context efficiency** | Sees summaries, not raw data | Sees full files, writes updates |
| **Decision authority** | Makes research decisions (continue/revise/validate/stop) | Executes decisions, escalates issues |

This separation prevents a single model from both designing bad experiments and then "validating" them with confirmation bias.

### 3.2 Iterative Loop

#### Full Research Loop Protocol

```
┌──────────────────────────────────────────┐
│           USER TRIGGERS STUDY            │
│    (goal + optional constraints)         │
└──────────────┬───────────────────────────┘
               ▼
┌──────────────────────────────────────────┐
│     PHASE 0: BOOTSTRAP (Opus)            │
│  • Create study folder structure         │
│  • Initialize README, IMPROVEMENTS.md    │
│  • Create study_state.json               │
│  • Lookup cqed_sim API                   │
└──────────────┬───────────────────────────┘
               ▼
┌──────────────────────────────────────────┐
│     PHASE 1: SCIENCE PLAN (Codex)        │
│  • Classify problem (OPT/REP/DES/ANA)   │
│  • Propose hypotheses                    │
│  • Design experiments                    │
│  • Define success criteria               │
│  • Produce SCIENCE_DIRECTIVE.md          │
└──────────────┬───────────────────────────┘
               ▼
┌──────────────────────────────────────────┐
│   PHASE 2: IMPLEMENT & EXECUTE (Opus)    │
│  • Read SCIENCE_DIRECTIVE.md             │
│  • Write simulation scripts              │
│  • Run simulations, save data            │
│  • Generate figures                      │
│  • Self-debug failures                   │
│  • Write EXECUTION_SUMMARY.md            │
└──────────────┬───────────────────────────┘
               ▼
┌──────────────────────────────────────────┐
│    PHASE 3: SCIENCE REVIEW (Codex)       │
│  • Evaluate physics correctness          │
│  • Check result quality                  │
│  • Identify gaps and improvements        │
│  • DECISION:                             │
│    ├─ CONTINUE → refine, extend          │
│    ├─ REVISE → new approach/hypothesis   │
│    ├─ VALIDATE → ready for validation    │
│    └─ STOP → needs human input           │
└──────────────┬───────────────────────────┘
               ▼
         ┌── Decision ──┐
         │              │
    CONTINUE/         VALIDATE
    REVISE               │
         │               ▼
    Back to        ┌──────────────────┐
    Phase 2        │ PHASE 4: VALIDATE│
                   │ • Sanity checks  │
                   │ • Convergence    │
                   │ • Literature     │
                   └────────┬─────────┘
                            ▼
                   ┌──────────────────┐
                   │ PHASE 5: REPORT  │
                   │ • Write LaTeX    │
                   │ • Compile PDF    │
                   │ • Finalize docs  │
                   │ • Mark COMPLETE  │
                   └──────────────────┘
```

#### Handoff Mechanism

Communication between agents is entirely **file-based** — no chat history is shared:

| Direction | File | Contents |
|-----------|------|----------|
| Science Director → Execution Engineer | `SCIENCE_DIRECTIVE.md` | Experiments to run, parameters, success criteria |
| Execution Engineer → Science Director | `EXECUTION_SUMMARY.md` | Results digest, anomalies, compute notes |
| Both → State | `study_state.json` | Machine-readable single source of truth |
| Both → History | `PROGRESS_LOG.md` | Append-only record of what happened |

#### Error Detection and Correction

| Error Type | Detection | Correction |
|------------|-----------|------------|
| **Code errors** | Traceback from script execution | Execution Engineer self-debugs (3 attempts) |
| **Physics errors** | Non-physical values, wrong units, failed sanity checks | Science Director flags in review, issues REVISE directive |
| **Missing controls** | Incomplete experiment coverage | Science Director identifies in review, adds to next directive |
| **Local minima** | Flat cost landscape, multistart comparison | Multiple initial seeds, alternative optimizers |
| **Unconverged results** | Fidelity change > threshold on parameter doubling | Re-run with increased resolution/truncation |

#### Convergence Determination

The loop terminates only when **all** research-quality stopping criteria pass:

- [ ] Scientific question answered with evidence
- [ ] Physics consistency verified (limiting cases, conservation laws)
- [ ] Diagnostics and controls complete (convergence, sanity checks)
- [ ] Results robust to small parameter perturbations
- [ ] Documentation complete (report.tex, IMPROVEMENTS.md, all figures)
- [ ] Open questions documented

### 3.3 Prompting Strategy

#### Structured Context, Not Prose

Agents receive structured state files rather than conversational context:

```
SCIENCE_DIRECTIVE.md         ← Ordered action items with parameters
  ┌─ ## Study Objective
  ├─ ## Problem Classification
  ├─ ## Hypotheses
  ├─ ## Experiment Design    ← Per-experiment: purpose, method, params, expected outcome
  ├─ ## Execution Plan       ← Numbered [IMPLEMENT]/[RUN]/[ANALYZE]/[DOCUMENT] tasks
  ├─ ## Assumptions
  └─ ## Stopping Criteria
```

#### Token Efficiency Rules

| Principle | Implementation |
|-----------|---------------|
| Codex sees summaries, not raw data | EXECUTION_SUMMARY.md capped at ~500 lines |
| Opus gets directives, not prose | Action lists with specific file paths and parameters |
| State persists in files, not context | study_state.json is the single source of truth |
| Large files never passed whole | Function signatures + key results only |
| Each iteration is self-contained | No reliance on prior chat history |
| Figures described, not embedded | Text description + file path, not base64 |

#### Per-Agent Context Strategy

| Agent | Gets | Does Not Get |
|-------|------|-------------|
| Science Director | study_state.json, EXECUTION_SUMMARY.md, figure paths, result digest | Full scripts, raw data arrays, compilation logs |
| Execution Engineer | SCIENCE_DIRECTIVE.md, full file access, API reference | Prior chat history, other studies' data |
| Research Loop | Everything (switches context based on current hat) | N/A — manages its own context |

### 3.4 AGENTS.md System

#### How Rules Are Defined

`AGENTS.md` is the master rule document that all agents read. It defines:

1. **5 Non-Negotiable Rules** — always use cqed_sim, never skip steps, install packages as needed, no venvs, no notebooks
2. **Problem Classes** — OPT/REP/DES/ANA taxonomy with typical deliverables
3. **Workflow Steps** — Initialize → Plan → Implement → Validate → Report
4. **Report Template** — Complete LaTeX template with mandatory sections
5. **Figure Standards** — 300 dpi PNG + vector PDF, colorblind-friendly palettes
6. **Validation Checklist** — Three required checks (sanity, convergence, literature)
7. **Decision Trees** — Can I use standalone code? Should I install a package? Is the study complete?

#### How Agents Follow Constraints

Agents are configured through VS Code's `.github/` directory:

```
.github/
├── agents/              ← Agent definition files (.agent.md)
│   ├── science-director.agent.md
│   ├── execution-engineer.agent.md
│   ├── research-loop.agent.md
│   ├── autonomous-planner.agent.md
│   ├── autonomous-implementer.agent.md
│   └── autonomous-resume.agent.md
├── prompts/             ← Slash command prompt files (.prompt.md)
│   ├── autonomous-plan.prompt.md
│   ├── autonomous-implement.prompt.md
│   └── autonomous-resume.prompt.md
├── instructions/        ← Context-specific instructions (.instructions.md)
│   └── task-run-state.instructions.md
└── skills/              ← Reusable capability definitions
    ├── cqed-sim-lookup/SKILL.md
    ├── latex-report/SKILL.md
    ├── publication-figures/SKILL.md
    ├── study-init/SKILL.md
    └── validate-results/SKILL.md
```

Each agent file specifies:
- **Description**: When to use this agent
- **Tools**: Which tools the agent can access (read, search, edit, execute, todo)
- **Argument hint**: Expected invocation format
- **System prompt**: Detailed behavioral instructions

#### How Updates Propagate

When AGENTS.md is updated, all agents automatically receive the new rules because:
1. AGENTS.md is attached to the workspace context
2. Agent definition files reference it ("Read AGENTS.md Quick Reference section")
3. Skills reference it ("AGENTS.md Step 5", "AGENTS.md validation checklist")

---

## Part 4 — VS Code Setup (Reproducible Guide)

### 4.1 Required Tools

#### Essential Software

| Tool | Version | Purpose | Installation |
|------|---------|---------|-------------|
| **VS Code** | Latest stable | IDE and agent host | [code.visualstudio.com](https://code.visualstudio.com) |
| **GitHub Copilot** | Latest | Agent infrastructure (chat, agents, skills) | VS Code marketplace extension |
| **Python** | 3.12.x | Runtime for simulations | [python.org](https://python.org) — system install, **no venvs** |
| **Git** | Latest | Version control | [git-scm.com](https://git-scm.com) |

#### Python Dependencies

```bash
pip install --user numpy scipy matplotlib qutip pandas seaborn lmfit
```

| Package | Version Requirement | Purpose |
|---------|-------------------|---------|
| `numpy` | ≥ 1.24 | Array operations |
| `scipy` | ≥ 1.10 | Optimization, linear algebra |
| `qutip` | ≥ 5.0 | Quantum simulation backend |
| `matplotlib` | ≥ 3.8 | Plotting |
| `pandas` | ≥ 2.0 | Data handling |
| `seaborn` | (optional) | Statistical visualization |
| `lmfit` | (optional) | Curve fitting |

#### cqed_sim Installation

```bash
# Install as editable package from local source
pip install --user -e /path/to/cQED_simulation
```

> **Windows note**: QuTiP 5.x import may hang due to a `platform._wmi_query` stall during NumPy import. Workaround: patch `platform._wmi_query` to raise `OSError` before importing qutip. The repository includes a `runtime_compat.py` helper that applies this automatically.

#### Optional Tools

| Tool | Purpose |
|------|---------|
| **LaTeX** (MiKTeX or TeX Live) | Compile PDF reports |
| **JAX** | GPU-accelerated GRAPE optimal control |

### 4.2 Model Usage Strategy

This workflow uses **three different AI models**, each optimized for a specific part of the research loop:

#### Copilot Chat (Codex-medium) — The Fast Iterator

```
┌─────────────────────────────────────────────┐
│  COPILOT CHAT (Codex-medium)                │
│                                             │
│  Best for:                                  │
│  • Quick code edits and fixes               │
│  • Running experiments interactively        │
│  • Ad-hoc debugging                         │
│  • Code completion and IntelliSense         │
│                                             │
│  Tradeoffs:                                 │
│  ✓ Fast (< 5 second response)              │
│  ✓ Cheap (low token cost)                  │
│  ✓ Good for routine coding tasks           │
│  ✗ Less reliable for complex reasoning     │
│  ✗ May hallucinate physics formulas         │
│  ✗ Limited multi-step planning             │
│                                             │
│  Use in VS Code:                            │
│  • Ctrl+I for inline edit                   │
│  • Copilot Chat sidebar                     │
│  • Code completions (Tab)                   │
└─────────────────────────────────────────────┘
```

#### OpenAI Codex (High Effort) — The Physics Validator

```
┌─────────────────────────────────────────────┐
│  OPENAI CODEX (High Effort)                 │
│                                             │
│  Best for:                                  │
│  • Physics-critical implementation          │
│  • Numerical correctness verification       │
│  • Experiment design and hypothesis testing │
│  • Sign convention and unit checking        │
│                                             │
│  Tradeoffs:                                 │
│  ✓ Strong physics reasoning                │
│  ✓ Good at spotting numerical issues       │
│  ✓ Reliable for high-stakes correctness    │
│  ✗ Slower response time                    │
│  ✗ Higher cost per query                   │
│  ✗ Limited tool access in some modes       │
│                                             │
│  Use as @science-director:                  │
│  • Plan and review phases                  │
│  • Evaluate physics correctness            │
│  • Design numerical experiments            │
└─────────────────────────────────────────────┘
```

#### Claude Opus 4.6 — The Implementation & Documentation Expert

```
┌─────────────────────────────────────────────┐
│  CLAUDE OPUS 4.6                            │
│                                             │
│  Best for:                                  │
│  • Writing full simulation scripts          │
│  • Multi-step implementation tasks          │
│  • LaTeX report generation                  │
│  • Structured documentation                 │
│  • Debugging complex failures               │
│  • File organization and state management   │
│                                             │
│  Tradeoffs:                                 │
│  ✓ Excellent long-context reasoning        │
│  ✓ Best readability and structure          │
│  ✓ Reliable multi-step execution           │
│  ✓ Full tool access (read/write/execute)   │
│  ✗ Slower than Copilot for trivial edits   │
│  ✗ Higher cost per session                 │
│                                             │
│  Use as @execution-engineer:               │
│  • Bootstrap, implement, validate, report  │
│  • All hands-on research engineering       │
└─────────────────────────────────────────────┘
```

#### Model Selection Decision Tree

```
Is this a physics reasoning or experiment design task?
├── YES → Use Codex (high effort) as @science-director
└── NO → Is this a multi-step implementation or report writing task?
         ├── YES → Use Opus 4.6 as @execution-engineer
         └── NO → Is this a quick edit, completion, or ad-hoc question?
                  ├── YES → Use Copilot Chat (Codex-medium)
                  └── Both? → Use @research-loop (switches between modes)
```

### 4.3 Recommended Workflow in VS Code

#### Step-by-Step: Running Your First Agent-Assisted Study

**1. Open the repository in VS Code**

```
File → Open Folder → select the cqed_based_study directory
```

**2. Verify the workspace is recognized**

Check that VS Code shows the agent definitions in the Copilot Chat sidebar. You should see available agents like `@science-director`, `@execution-engineer`, `@research-loop`.

**3. Option A: Full Autonomous Loop (Recommended)**

In the Copilot Chat panel:

```
@research-loop study=studies/my_new_study goal='Optimize dispersive shift for 99.5% readout fidelity'
```

The research loop agent will:
1. Create the study folder structure
2. Plan experiments (Science Director hat)
3. Write and run simulation code (Execution Engineer hat)
4. Review results (Science Director hat)
5. Iterate until validated
6. Write the final report

**4. Option B: Step-by-Step Control**

For more control, drive each phase manually:

```bash
# Step 1: Initialize the study
# In VS Code: Terminal → Run Task → "Research: New Study"
# Or in chat:
@execution-engineer study=studies/my_study run=task_runs/my_study phase=bootstrap
```

```bash
# Step 2: Science Director plans
@science-director study=studies/my_study run=task_runs/my_study phase=plan
```

```bash
# Step 3: Execution Engineer implements
@execution-engineer study=studies/my_study run=task_runs/my_study phase=implement
```

```bash
# Step 4: Science Director reviews
@science-director study=studies/my_study run=task_runs/my_study phase=review
```

```bash
# Step 5: Repeat steps 3-4 until VALIDATE decision, then:
@execution-engineer study=studies/my_study run=task_runs/my_study phase=validate
@execution-engineer study=studies/my_study run=task_runs/my_study phase=report
```

**5. Check Study Status**

```bash
# Via VS Code task:
Terminal → Run Task → "Research: Study Status"

# Or via PowerShell:
.\tools\research_loop.ps1 -Action status -StudyName "my_study"
```

**6. Resume an Interrupted Study**

```
@research-loop study=studies/my_study resume
```

Or:
```bash
.\tools\research_loop.ps1 -Action resume -StudyName "my_study"
```

### 4.4 Automated Loop (Advanced)

#### PowerShell Research Loop Script

The `tools/research_loop.ps1` script provides command-line orchestration:

```powershell
# Initialize a new study
.\tools\research_loop.ps1 -Action init -StudyName "chi_optimization" -StudyGoal "Optimize chi for readout"

# Check status
.\tools\research_loop.ps1 -Action status -StudyName "chi_optimization"

# Resume interrupted study
.\tools\research_loop.ps1 -Action resume -StudyName "chi_optimization"

# Run specific phase
.\tools\research_loop.ps1 -Action execute -StudyName "chi_optimization"
.\tools\research_loop.ps1 -Action validate -StudyName "chi_optimization"
.\tools\research_loop.ps1 -Action report -StudyName "chi_optimization"
```

#### VS Code Tasks

Available from **Terminal → Run Task**:

| Task | What It Does |
|------|-------------|
| **Research: New Study** | Initialize study folder + state files |
| **Research: Study Status** | Show current loop state and next action |
| **Research: Resume Study** | Detect phase and continue |
| **Research: Run Loop Action** | Pick any phase to run |
| **Copilot: Init Task Run** | Bootstrap state files for resumable task |
| **Copilot: Show Task Run Status** | Display task run state |

#### Retry on Failure

The self-debugging protocol handles most failures automatically. For persistent blockers:

1. Check `BLOCKERS.md` in the task run directory
2. Address the blocker manually or provide guidance
3. Resume: `@research-loop study=studies/<name> resume`

#### Long-Running Study Execution

For studies that take multiple sessions:

1. The file-based state system (`study_state.json`, `TASK_CHECKLIST.md`, `PROGRESS_LOG.md`) survives VS Code restarts
2. The `autonomous-resume` agent reconstructs context from files, not chat history
3. Each iteration is self-contained — no dependence on prior conversation

### 4.5 Terminal / Execution Integration

#### Running Experiments from VS Code Terminal

All simulation scripts are designed to run from the terminal:

```powershell
# Navigate to study scripts
cd studies/gray_box_adaptive_control/scripts

# Run a simulation phase
python study_phase4.py

# Run validation
python validate_results.py
```

#### Monitoring Progress

- **study_state.json** — machine-readable current state
- **PROGRESS_LOG.md** — human-readable chronological record
- **TASK_CHECKLIST.md** — checkbox-based completion tracking
- **Console output** — real-time simulation progress

#### Resuming Interrupted Tasks

```powershell
# Check what's pending
Get-Content task_runs/my_study/TASK_CHECKLIST.md | Select-String "^\- \[ \]"

# Check blockers
Get-Content task_runs/my_study/BLOCKERS.md

# Resume via agent
# @autonomous-resume task=... run=task_runs/my_study
```

---

## Part 5 — Generalization to Other Fields

The agent-assisted research workflow is **domain-agnostic in design** — the cQED-specific components can be replaced while preserving the full orchestration infrastructure.

### What Stays the Same

These components transfer directly to any computational physics domain:

| Component | Description | Files |
|-----------|-------------|-------|
| **Agent infrastructure** | VS Code agents, prompts, skills, instructions | `.github/agents/`, `.github/prompts/`, `.github/skills/` |
| **Study folder structure** | Standardized `studies/<name>/` layout | AGENTS.md template |
| **State management** | study_state.json, TASK_CHECKLIST.md, PROGRESS_LOG.md, BLOCKERS.md | `task_runs/<name>/` |
| **Two-model loop** | Science Director (physics reasoning) + Execution Engineer (implementation) | RESEARCH_LOOP.md |
| **Self-debugging protocol** | 4-level escalation: inspect → fix → log → stop | Agent definitions |
| **Report template** | LaTeX with mandatory sections (abstract through appendix) | `.github/skills/latex-report/` |
| **Validation framework** | 3-check gate: sanity, convergence, literature | `.github/skills/validate-results/` |
| **IMPROVEMENTS.md pattern** | Living improvement log with P1/P2/P3 priority tags | AGENTS.md specification |
| **Figure standards** | PNG + PDF, colorblind-friendly, labeled axes with units | `.github/skills/publication-figures/` |
| **Automation scripts** | PowerShell orchestration, VS Code tasks | `tools/`, `.vscode/tasks.json` |

### What Changes Per Domain

| Component | cQED Version | Adaptation Required |
|-----------|-------------|---------------------|
| **Physics simulator** | `cqed_sim` (QuTiP-based) | Replace with domain-specific solver |
| **API Reference** | `API_REFERENCE.md` for cqed_sim | Write equivalent for new simulator |
| **Skill: sim-lookup** | `cqed-sim-lookup/SKILL.md` | Rewrite with new simulator's API |
| **Physics conventions** | Dispersive Hamiltonian, χ, Kerr, etc. | Domain-specific Hamiltonian and conventions |
| **Validation checks** | cQED-specific sanity checks (unitarity, dispersive limit) | Domain-specific limiting cases |
| **Default parameters** | Transmon/cavity frequencies, coupling strengths | Domain-specific parameter table |
| **Problem classes** | OPT/REP/DES/ANA (cQED-oriented) | May need domain-specific categories |
| **Agent physics knowledge** | cQED expertise in Science Director prompt | Update system prompt with domain expertise |

### Domain-Specific Adaptation Guide

#### Plasma Physics

| cQED Component | Plasma Physics Replacement |
|---------------|---------------------------|
| `cqed_sim` | PlasmaPy, GENE, or custom MHD solver |
| Hamiltonian construction | MHD equations, Vlasov-Poisson system |
| Noise modeling (T1, T2, κ) | Collisional transport, resistivity, radiation losses |
| Dispersive shift validation | Alfvén wave dispersion, MHD stability criteria |
| Fidelity metrics | Energy conservation, momentum conservation, growth rates |

#### Nuclear Theory

| cQED Component | Nuclear Theory Replacement |
|---------------|---------------------------|
| `cqed_sim` | LAMMPS, Geant4, nuclear shell model codes |
| QuTiP solver | Many-body Schrödinger equation solvers, density functional theory |
| Pulse sequences | External field protocols, bombarding energy profiles |
| Validation | Nuclear mass tables, cross-section databases, Bethe-Weizsäcker |

#### Condensed Matter

| cQED Component | Condensed Matter Replacement |
|---------------|-------------------------------|
| `cqed_sim` | Kwant, PySCF, VASP interfaces, tight-binding solvers |
| Dispersive model | Band structure Hamiltonians, Hubbard model |
| Readout simulation | Transport calculation, optical spectra |
| Validation | Band gap databases, DFT benchmarks |

#### AMO Physics

| cQED Component | AMO Replacement |
|---------------|-----------------|
| `cqed_sim` | QuTiP directly (no cQED wrapper needed), ARC (alkali Rydberg calculator) |
| Transmon/cavity model | Atom + laser field, Rydberg blockade |
| Noise (T1, T2) | Spontaneous emission, Doppler broadening, collisional dephasing |
| Validation | Atomic spectral databases (NIST), known transition strengths |

### Minimal Adaptation Checklist

To adapt this workflow to a new domain:

- [ ] **Replace the simulator**: Install your domain's solver, write an `API_REFERENCE.md`
- [ ] **Update the sim-lookup skill**: Rewrite `SKILL.md` with new API documentation
- [ ] **Update AGENTS.md**: Replace cQED-specific rules, default parameters, and problem classes
- [ ] **Update Science Director prompt**: Add domain expertise to the agent definition
- [ ] **Update validation skill**: Replace sanity checks with domain-appropriate tests
- [ ] **Update figure skill**: Adjust default axes labels and color conventions (if needed)
- [ ] **Keep everything else**: The study structure, state management, two-model loop, report template, and automation scripts work as-is

---

## Part 6 — Recommendations & Improvements

### 6.1 Strengths

| Strength | Evidence |
|----------|---------|
| **Reproducibility** | All 6 studies re-run successfully; 69/69 validation checks pass |
| **Structured output** | Standardized study folders, mandatory appendices, living improvement logs |
| **Failure recovery** | File-based state survives interruptions; resume agent reconstructs context |
| **Physics correctness** | Two-model separation prevents self-confirming validation |
| **Documentation quality** | LaTeX reports with mandatory limitations and future work sections |
| **Knowledge preservation** | IMPROVEMENTS.md and failed-approach logs prevent repeating mistakes |

### 6.2 Weaknesses

#### Agent Reliability

| Issue | Impact | Severity |
|-------|--------|----------|
| **Hallucinated physics** | Agents may generate plausible but incorrect formulas | Medium — caught by validation checks, but not always |
| **Code-first, think-later** | Execution Engineer may write code before fully understanding the physics | Medium — mitigated by Science Director review |
| **Context window limits** | Large studies may exceed model context, causing lost state | Low — file-based state compensates |
| **Non-deterministic outputs** | Same prompt may produce different code across runs | Low — validation catches issues, but complicates exact reproducibility |

#### Validation Gaps

| Gap | Description | Risk |
|-----|-------------|------|
| **No experimental data** | All validation is simulation-to-theory, not simulation-to-experiment | Medium — results are internally consistent but may not match real devices |
| **Limited cross-validation** | Only QuTiP + dense backend comparison; no third-party simulator | Low — QuTiP is well-validated |
| **Automated validation coverage** | Not all studies have full validation scripts | Medium — dispersive readout and SQR gate lack programmatic validation |

#### Simulation-to-Experiment Mismatch

| Factor | Current Treatment | Gap |
|--------|-------------------|-----|
| **Decoherence** | Aggregate T1, T2, κ | No hardware-calibrated noise model |
| **Fabrication variation** | Fixed parameter set | No Monte Carlo over parameter uncertainty |
| **Crosstalk** | HardwareConfig supports it | Not exercised in most studies |
| **Measurement back-action** | Strong-readout model available | Not universally applied |

### 6.3 Recommended Improvements

#### P1 — Critical

| Improvement | Description | Difficulty |
|-------------|-------------|-----------|
| **Experimental validation loop** | Connect to lab data for at least one study | HIGH |
| **Automated end-to-end CI** | Run all study validation scripts in CI pipeline | MEDIUM |
| **Physics unit tests for agent output** | Automatically check generated code for unit consistency | MEDIUM |

#### P2 — Meaningful

| Improvement | Description | Difficulty |
|-------------|-------------|-----------|
| **Multi-agent orchestration API** | Replace file-based handoff with structured API calls | MEDIUM |
| **Parameter uncertainty propagation** | Monte Carlo sweeps over fabrication variation | MEDIUM |
| **Science Director memory** | Persistent cross-study knowledge base for the physics agent | LOW |
| **Automated figure quality check** | Verify axes labels, units, colorblind-friendliness | LOW |
| **LaTeX compilation in CI** | Auto-compile reports on study completion | LOW |

#### P3 — Nice-to-Have

| Improvement | Description | Difficulty |
|-------------|-------------|-----------|
| **Web dashboard** | Real-time study status visualization | HIGH |
| **Cross-study meta-analysis** | Agent that synthesizes findings across related studies | HIGH |
| **Natural language query interface** | "What's the optimal chi for 99% readout fidelity?" | MEDIUM |
| **Automatic theme detection** | Cluster related studies into themes automatically | MEDIUM |

### 6.4 Orchestration Improvements

The current workflow relies on VS Code Copilot Chat for agent invocation, which has limitations:

1. **No true background execution** — agents run in the chat context, not as background services
2. **Manual phase transitions** — the user must invoke each phase (or use the research-loop agent)
3. **No parallel studies** — one study at a time per VS Code window

**Proposed improvement**: A lightweight process manager that:
- Runs phase transitions automatically based on study_state.json
- Supports multiple concurrent studies
- Provides a CLI for monitoring and intervention
- Integrates with GitHub Actions for CI-triggered study runs

### 6.5 Reproducibility Improvements

| Current State | Improvement |
|--------------|-------------|
| System Python (no venv) | Pin exact package versions in `requirements.txt` |
| Manual cqed_sim install | Add `pyproject.toml` with cqed_sim as a dependency |
| Windows-specific workarounds | Cross-platform CI testing (Linux + macOS) |
| No data versioning | DVC or Git LFS for large data files |

---

## Appendix A — File Tree Reference

### Complete Repository Structure

```
cqed_based_study/
├── AGENTS.md                      ← Master agent rules document
├── README.md                      ← Workspace overview
├── RERUN_SUMMARY.md               ← Full re-run verification report
├── RESEARCH_LOOP.md               ← Two-model loop architecture
│
├── .github/
│   ├── agents/                    ← Agent definition files
│   │   ├── science-director.agent.md
│   │   ├── execution-engineer.agent.md
│   │   ├── research-loop.agent.md
│   │   ├── autonomous-planner.agent.md
│   │   ├── autonomous-implementer.agent.md
│   │   └── autonomous-resume.agent.md
│   ├── prompts/                   ← Slash command definitions
│   │   ├── autonomous-plan.prompt.md
│   │   ├── autonomous-implement.prompt.md
│   │   └── autonomous-resume.prompt.md
│   ├── instructions/              ← Context-specific rules
│   │   └── task-run-state.instructions.md
│   └── skills/                    ← Reusable capability packs
│       ├── cqed-sim-lookup/SKILL.md
│       ├── latex-report/SKILL.md
│       ├── publication-figures/SKILL.md
│       ├── study-init/SKILL.md
│       └── validate-results/SKILL.md
│
├── .vscode/
│   ├── settings.json              ← Workspace settings
│   └── tasks.json                 ← VS Code task definitions
│
├── tools/
│   ├── research_loop.ps1          ← Research loop CLI orchestrator
│   ├── copilot_task_run.ps1       ← Task run bootstrapper
│   └── test_params.ps1            ← Parameter testing utility
│
├── studies/
│   ├── dispersive_readout_optimization/
│   ├── gray_box_adaptive_control/
│   ├── hybrid_qubit_cavity_control/
│   ├── literature_informed_selective_primitives/
│   ├── sqr_gate_design/
│   ├── thermal_noise_cavity_sensing/
│   └── themes/                    ← Cross-study synthesis documents
│
└── task_runs/
    └── thermal_noise_cavity_sensing/   ← State files for active run
        ├── BLOCKERS.md
        ├── PROGRESS_LOG.md
        └── TASK_CHECKLIST.md
```

### cqed_sim Framework (External)

```
cQED_simulation/                   ← Source-of-truth simulator
├── API_REFERENCE.md               ← Canonical API documentation
├── cqed_sim/                      ← Python package
│   ├── core/                      ← Models, frames, gates, states
│   ├── pulses/                    ← Pulse construction + calibration
│   ├── sequence/                  ← Waveform compilation
│   ├── sim/                       ← Simulation engine
│   ├── backends/                  ← QuTiP / NumPy / JAX solvers
│   ├── measurement/               ← Readout modeling
│   ├── floquet/                   ← Periodic-drive analysis
│   ├── analysis/                  ← Parameter translation
│   ├── calibration/               ← SQR gate calibration
│   ├── calibration_targets/       ← Surrogate experiments
│   ├── gates/                     ← 100+ ideal unitary gates
│   ├── operators/                 ← Cached quantum operators
│   ├── observables/               ← Diagnostics and metrics
│   ├── plotting/                  ← Visualization
│   ├── tomo/                      ← Tomography
│   ├── unitary_synthesis/         ← Gate-sequence optimal control
│   ├── optimal_control/           ← Direct GRAPE
│   ├── rl_control/                ← RL environment
│   └── io/                        ← JSON I/O
├── physics_and_conventions/       ← Sign conventions documentation
└── tests/                         ← 58+ test files
```

---

## Appendix B — Complete Agent Definition Files

### Agent: Science Director

- **Model**: OpenAI Codex / GPT (high-effort)
- **Tools**: read, search, todo (read-only)
- **Invocation**: `@science-director study=studies/<name> run=task_runs/<slug> phase=plan|review`
- **Two phases**:
  - **PLAN**: Reads study goal + API capabilities → produces SCIENCE_DIRECTIVE.md with hypotheses, experiment design, execution plan, assumptions, stopping criteria
  - **REVIEW**: Reads EXECUTION_SUMMARY.md → evaluates physics correctness → decides CONTINUE/REVISE/VALIDATE/STOP → produces updated SCIENCE_DIRECTIVE.md

### Agent: Execution Engineer

- **Model**: Claude Opus 4.6
- **Tools**: read, search, edit, execute, todo (full access)
- **Invocation**: `@execution-engineer study=studies/<name> run=task_runs/<slug> phase=bootstrap|implement|validate|report`
- **Four phases**:
  - **BOOTSTRAP**: Create study folder structure, initialize state files
  - **IMPLEMENT**: Read directive → write scripts → run simulations → generate figures → write EXECUTION_SUMMARY.md
  - **VALIDATE**: Run 3-check validation gate (sanity, convergence, literature)
  - **REPORT**: Write report.tex → compile PDF → finalize IMPROVEMENTS.md → mark COMPLETE

### Agent: Research Loop

- **Model**: Combined (switches hats between Science Director and Execution Engineer)
- **Tools**: read, search, edit, execute, todo (full access)
- **Invocation**: `@research-loop study=studies/<name> goal='Research question'` or `study=studies/<name> resume`
- **Manages the full loop autonomously**, switching between physics reasoning and implementation

### Agent: Autonomous Planner

- **Tools**: read, search, edit, todo
- **Invocation**: `/Autonomous Plan task=<path> run=<dir>`
- **Creates**: EXECUTION_PLAN.md, TASK_CHECKLIST.md, PROGRESS_LOG.md, BLOCKERS.md

### Agent: Autonomous Implementer

- **Tools**: read, search, edit, execute, todo
- **Invocation**: `/Autonomous Implement task=<path> run=<dir>`
- **Executes**: Next unchecked item(s) from TASK_CHECKLIST.md

### Agent: Autonomous Resume

- **Tools**: read, search, edit, execute, todo
- **Invocation**: `/Autonomous Resume task=<path> run=<dir>`
- **Recovers**: State from files, continues from exact stopping point

---

## Appendix C — Glossary

| Term | Definition |
|------|-----------|
| **Agent** | An AI model instance configured with specific tools, instructions, and system prompts for a defined role |
| **cqed_sim** | The source-of-truth circuit QED simulation framework built on QuTiP |
| **cQED** | Circuit quantum electrodynamics — the study of superconducting qubits coupled to microwave resonators |
| **Dispersive regime** | Operating regime where qubit-cavity coupling is much weaker than their detuning ($g \ll \Delta$) |
| **GRAPE** | Gradient Ascent Pulse Engineering — an optimal control algorithm for quantum gate synthesis |
| **IMPROVEMENTS.md** | Living document logging limitations, failed approaches, and actionable improvements |
| **Opus 4.6** | Claude Opus 4.6 — the large language model used as the Execution Engineer |
| **Problem class** | Classification of study type: OPT (optimization), REP (reproduction), DES (design), ANA (analysis) |
| **QuTiP** | Quantum Toolbox in Python — the numerical backend for quantum dynamics simulation |
| **Research Loop** | The continuous Plan → Implement → Review → (Iterate or Validate) → Report cycle |
| **Science Director** | The physics-reasoning agent role (OpenAI Codex/GPT) |
| **Execution Engineer** | The implementation agent role (Claude Opus 4.6) |
| **Skill** | A reusable capability definition that provides domain-specific instructions to agents |
| **SQR** | Selective Quantum Rotation — a class of qubit-cavity conditional operations |
| **study_state.json** | Machine-readable JSON file serving as the single source of truth for study progress |
| **Transmon** | A type of superconducting qubit with reduced charge noise sensitivity |
| **χ (chi)** | Dispersive shift — the photon-number-dependent qubit frequency shift |
| **κ (kappa)** | Cavity (resonator) decay rate / linewidth |

---

*Report generated: 2026-03-23*
*Repository: cqed_based_study (Shyam Shankar Quantum Circuits Group)*
*Simulator: cqed_sim v0.1.0*
