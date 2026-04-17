# Realistic Universal Hybrid Control in the Dispersive cQED Regime

## Problem Class
ANA | DES

## Motivation
The ideal universal-control claim for a hybrid transmon-cavity system is usually framed in terms of three primitives:

1. arbitrary selective qubit rotations,
2. cavity displacements,
3. unconditional qubit rotations.

This study asks whether that claim survives the realistic dispersive Hamiltonian once the leading dispersive shift `chi`, the higher-order correction `chi'`, cavity self-Kerr, and finite-duration phase accumulation are enforced. The objective is not to add one more black-box optimizer pass, but to identify what the control problem itself allows or forbids and to determine whether a physically meaningful non-GRAPE route to universal control actually exists.

## Problem Class Rationale
- `ANA`: the study is primarily a first-principles and evidence-synthesis analysis of what is and is not physically realizable.
- `DES`: the study is also protocol-facing because it identifies which structured pulse families remain credible building blocks for future experimental design.

## Goals
1. Derive the timing and phase-accumulation constraints imposed by the realistic dispersive Hamiltonian.
2. Consolidate the validated repository evidence on unconditional displacement, spectator qubit rotations, practical selective pulses, strict ideal SQR, relaxed CPSQR, arbitrary conditional control, and sequence-level replay.
3. Decide whether the strict ideal primitive gate set survives once realistic Hamiltonian dynamics are enforced.
4. State explicitly which primitives survive only approximately, which survive only after gauge relaxation, and which require replacement by a stronger native primitive.
5. Produce a final report, machine-readable summary, figures, and a reproducibility notebook that future agents can use without re-running every upstream study.

## Study Composition
This study is a synthesis layer over the following completed studies:

| Inherited component | Role in this synthesis |
|---|---|
| `waveform_level_gate_realization_dispersive_cqed` | Unconditional displacement, spectator qubit rotations, echoed displacement negative result |
| `literature_informed_selective_primitives` | Best noisy practical SQR and geometric SNAP references |
| `native_rich_multitone_sqr_cpsqr_feasibility` | Strict full-joint ideal SQR versus relaxed CPSQR |
| `strong_validation_arbitrary_fock_conditional_rotations` | Arbitrary blockwise SU(2) versus left-Z-gauge-relaxed conditional control |
| `hybrid_unitary_native_entangling_evolution` | Sequence-level replay readiness and runtime limitations |

## Methods
- `cqed_sim` API audit via the local `API_REFERENCE.md` and installed package inspection.
- Structured data aggregation from the inherited study artifacts and summaries.
- First-principles analytic timing and phase-budget calculation for:
  - photon-number-dependent qubit detuning through `chi` and `chi'`,
  - cavity branch mismatch during displacement,
  - Kerr-induced phase shear.
- Unified machine-readable output generation through:
  - [`build_synthesis_dataset.py`](scripts/build_synthesis_dataset.py)
  - [`build_reproducibility_notebook.py`](scripts/build_reproducibility_notebook.py)
- Publication figure generation:
  - [`timescale_hierarchy`](figures/timescale_hierarchy.pdf)
  - [`phase_budget`](figures/phase_budget.pdf)
- Report writing and compilation to PDF in [`report/report.tex`](report/report.tex)

## Analytic Preliminary
The common effective Hamiltonian is

\[
\frac{H}{\hbar} =
\omega_c a^\dagger a
+ \omega_q b^\dagger b
+ \frac{\alpha}{2} b^{\dagger 2} b^2
+ \chi\, a^\dagger a\, b^\dagger b
+ \chi' a^\dagger a(a^\dagger a - 1)b^\dagger b
+ \frac{K}{2} a^\dagger a(a^\dagger a - 1)
+ H_d^{(q)}(t)
+ H_d^{(c)}(t).
\]

Two competing control conditions follow immediately.

For qubit control on cavity manifold `n`,

\[
\delta_q(n)/2\pi \approx n\chi/2\pi + n(n-1)\chi'/2\pi.
\]

Selective number resolution therefore prefers `|chi| T / 2pi ~ 1` or larger, while spectator-indifferent qubit control prefers `|delta_q(n)| T / 2pi << 1` across the populated manifold window.

For cavity displacement, branch-independent action likewise requires `|chi| T / 2pi << 1`, while Kerr-induced shear scales as `|K| \bar{n} T / 2pi`.

The study's main analytic hypothesis is therefore:
- the same dispersive shift that enables number selectivity also obstructs literal unconditionality,
- so the realistic primitive library must become phase aware rather than remain identical to the ideal abstraction.

## cqed_sim Gap Analysis
| Functionality | Needed? | Available in cqed_sim? | Plan |
|---|---|---|---|
| Dispersive Hamiltonian conventions for `chi`, `chi'`, and Kerr | Yes | Yes | Use the local API reference and inherited validated studies |
| Machine-readable primitive-level results | Yes | Partial | Load from completed studies where available; keep synthesis-local summary JSON |
| Waveform bridge for `SNAP` and `FreeEvolveCondPhase` | Needed for future full universal stack, not for this synthesis | No public bridge | Record as a continuing blocker for future pulse-backed universality work |
| `ConditionalDisplacement` primitive support | Relevant to broader universal-stack follow-up | Not available as a public synthesis gate in the inherited evidence | Document as future-work gap rather than a blocker for this synthesis |
| New heavy pulse simulation | No | N/A | Reuse validated inherited artifacts instead of duplicating expensive studies |

## Assumptions
- The inherited repository studies are the source of record for their own underlying simulations and validations.
- The common nominal device point is:
  - `omega_q / 2pi = 6.150 GHz`
  - `omega_c / 2pi = 5.241 GHz`
  - `alpha / 2pi = -255 MHz`
  - `chi / 2pi = -2.84 MHz`
  - `chi' / 2pi = -21 kHz`
  - `K / 2pi = -28 kHz`
- This study evaluates universality at the level of primitive viability and sequence-level inherited evidence, not by claiming a new global optimum over all conceivable control families.
- GRAPE is treated as a benchmark or assistive lower bound, not as the conceptual answer.

## Compute & Resource Strategy
- Upfront estimate: this synthesis should be lightweight because it aggregates existing artifacts rather than rerunning the full repository.
- Chosen strategy: reuse validated JSON summaries and only add new analytic post-processing plus figure/report generation.
- No GPU backend was needed.
- No new packages were installed.
- Realized cost:
  - synthesis data build: a few seconds,
  - notebook generation: under one second,
  - notebook execution via `python -m jupyter nbconvert --execute`: about five seconds,
  - report compilation: a few `pdflatex` / `bibtex` passes on CPU.

## Expected Outcomes
- A precise claim boundary for the ideal primitive gate set.
- A machine-readable summary of which primitives survive literally, approximately, or only after gauge relaxation.
- A practical recommendation for the strongest current phase-aware constructive route.
- A scientifically honest statement of why a fully validated non-GRAPE universal-control stack is still open.

## Known Limitations
- This is a synthesis study, not a new end-to-end re-optimization over every control family.
- Most inherited strict-control studies are closed-system and optimize in a two-level transmon loop.
- The sequence-level runtime conclusion is strongest for low-dimensional logical windows.
- The report deliberately avoids treating unconstrained GRAPE as the conceptual answer, so it does not claim a tight ultimate upper bound on constrained control performance.

## Validation
- [x] Sanity checks
  - The derived timescale hierarchy matches the inherited numerical operating points: near-unconditional primitives live at `|chi| T / 2pi << 1`, while successful selective primitives live at `|chi| T / 2pi ~ 1-5`.
  - The integrated verdict is consistent with the inherited primitive results:
    - short two-tone displacement succeeds as an approximate unconditional primitive,
    - vacuum-calibrated qubit rotations fail as truly unconditional operations once photons occupy the cavity,
    - relaxed selective control is much stronger than strict full-joint selective control.
- [x] Convergence
  - The new synthesis uses only saved validated artifacts plus deterministic analytic post-processing.
  - The inherited studies already carried their own truncation and time-step checks; this study does not introduce a new heavy solver loop that would require a separate Hilbert-space convergence sweep.
- [x] Literature comparison (if applicable)
  - The timing conflict is consistent with the standard dispersive number-resolution logic from the transmon, SNAP, and bosonic-control literature cited in the report.
  - The phase-aware interpretation is consistent with the inherited repository evidence that conditional-phase-relaxed control is much easier than strict full-joint control.
- [x] Reproducibility notebook
  - `scripts/reproducibility_notebook.ipynb` was generated and executed successfully through `python -m jupyter nbconvert --to notebook --execute --inplace`.

## Suggested Upstreaming
- Add first-class waveform support for `SNAP`, `FreeEvolveCondPhase`, and future conditional-displacement-like primitives so a fully pulse-backed universal stack can be tested fairly.
- Add package-level helpers for primitive-level claim-boundary synthesis from validated study artifacts.
- Add a public phase-aware composite-sequence objective that targets strict joint control and relaxed conditional-phase control side by side.

## Status
COMPLETE
