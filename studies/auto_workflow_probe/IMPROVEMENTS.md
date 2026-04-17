# Improvement Log: auto_workflow_probe

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | MEDIUM]** Rewrite the report around substantive sections and integrate the independent cross-check outputs. The current report still contains placeholder scaffold text, so the new evidence will not matter until the report phase replaces the template language.

## Recommended Improvements (P2)
- **[P2 | LOW]** Repeat the same probe with nonzero $\chi'$ or Kerr terms after the minimal run succeeds, to show explicitly how the ideal $t_{\pi} = \pi/|\chi|$ result deforms outside the strictly first-order low-occupation limit.
- **[P2 | LOW]** If a later review still wants a dynamics-level check beyond the static Hamiltonian exponential, compile an explicit idle segment through the sequence simulator and compare it against the current cross-check artifact.

## Nice-to-Haves (P3)
- **[P3 | LOW]** Add a tiny reproducibility-notebook cell that reloads the saved artifact and recomputes the analytic $\pi$ time as a smoke test for the notebook path.

## Open Questions
- None yet.

## What Was Tried and Did Not Work
- The first revision attempt evaluated the static-Hamiltonian cross-check with `frame=None`, which retained the bare cavity carrier and produced a wrapped disagreement of about `3.05 rad` against the dispersive phase law. Re-running the same cross-check in the rotating frame defined by the bare qubit and cavity frequencies fixed the mismatch and isolated the intended dispersive phase.

## Compute & Resource Notes
- Planned compute budget: keep the implementation phase below 2 minutes on CPU.
- Expected runtime for the actual probe script is well under 10 seconds because only small dense unitaries and a single simple figure are required.
- Implementation update: `scripts/free_dispersive_pi_probe.py` now performs the helper-based phase scan with `n_cav=2` and `n_cav=3`, writes one CSV plus one JSON artifact, and saves one dual-format figure.
- Measured implementation runtime was 1.12 s on CPU for the full script, including CSV, JSON, PNG, and PDF writes.
- Validation update: the saved artifact confirms a maximum wrapped analytic-versus-helper mismatch of 8.882e-16 rad and a zero wrapped difference between `n_cav=2` and `n_cav=3`; no additional time-step convergence sweep was required because the probe uses a direct dispersive helper evaluation.
- Revision implementation update: the same script now also evaluates explicit static-Hamiltonian evolution on the three-level cavity cutoff, writes a second CSV plus JSON artifact, and saves a second dual-format comparison figure.

## Resolved
- **[P1 | LOW]** The reviewer's helper-only evidence concern was addressed during the revision implement phase by adding a static-Hamiltonian evolution cross-check, recorded in `artifacts/free_dispersive_hamiltonian_cross_check.json` and the matching comparison figure.
