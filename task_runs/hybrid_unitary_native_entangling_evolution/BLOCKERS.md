# Blockers

## Active
- No active blockers. The runtime-validation milestone is complete; remaining work concerns surrogate quality, duration reduction, and broader robustness studies rather than a missing execution path.

## Resolved
- Built a replay-backed validation path for the top native-heavy candidates by replacing the exact qubit Hadamard with replayable `QubitRotation` pairs, replaying `FreeEvolveCondPhase` as explicit idle evolution, and substituting GRAPE-derived replayable surrogates for the non-replayable `exact_hc` / `A_local` local blocks.
- Established the exact U_target and logical basis from the earlier hybrid study before creating the new study, avoiding a duplicated or inconsistent target definition.
- Confirmed that Bloch and Wigner diagnostics are available in the installed cqed_sim API, so later diagnostic phases can stay within the framework.
- Re-verified replay support claims in the installed cqed_sim build. The check is complete; the result is a concrete replay-path blocker rather than an unknown API question.
- Completed the requested symbolic depth diagnostics, so the study is no longer blocked on missing Bloch-versus-depth or Wigner-versus-depth analysis.