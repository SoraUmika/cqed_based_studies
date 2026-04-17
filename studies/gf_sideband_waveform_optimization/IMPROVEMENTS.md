# Improvement Log: Waveform Optimization for the `|f,n-1> <-> |g,n>` Sideband Interaction in cQED

> Written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | HIGH]** No microscopic pump model in the native workflow: The final ranking still uses the effective sideband-control abstraction available in the simulator. That is adequate for waveform ranking inside the present control layer, but it does not derive the sideband rate, AC Stark shift, or parasitic channels from the physical pump line. A future study should replace or augment the effective sideband operator with a pump-aware Hamiltonian.
- **[P1 | HIGH]** The sideband-reset example does not export device-matched transmon coherence values: The current transmon-noise study uses a matched local tomography workflow as a sensitivity anchor (`T1 = 9.813 us`, `T2,Ramsey = 6.325 us`). That is much better than ignoring transmon decoherence, but it is still not a device-verified parameter set for the exact sideband example. A future experiment-facing pass should attach measured transmon `T1/T2` to the same device tuple used for the sideband study.
- **[P1 | HIGH]** Constant simultaneous two-tone driving is not a generic multi-photon transfer primitive: The new extension shows that the clean `|g,0,1> -> |g,1,0>` bright-state picture is specific to the single-photon manifold. For `n > 1`, the constant overlapping drive opens a longer conversion chain and the nominal target population stalls around `0.61-0.69` even in the closed system. Any future protocol aimed at arbitrary storage Fock states must explicitly break or control that larger ladder.

## Recommended Improvements (P2)
- **[P2 | MEDIUM]** Add unitary-level optimization on top of the new gate-oriented diagnostic: The study now exports a gate-oriented ranking and shows explicitly that no analytic family in the current set is a phase-clean coherent SWAP across `n = 1,2,3`. The next step is no longer to add the diagnostic, but to optimize directly against it.
- **[P2 | MEDIUM]** Upstream the study-local envelope optimizer and noise-scenario utilities: The current study now contains reusable local helpers for arbitrary family replay, transmon-noise scenario generation, and projected-subspace diagnostics. Turning those into public simulator utilities would reduce duplicated study-local code in future waveform projects.
- **[P2 | HIGH]** Replace the constant two-tone square pulse with shaped counter-intuitive overlap: The detuned single-photon Raman-like case does suppress the intermediate `|f>` occupation to about `0.12`, but it takes `348 ns` and falls to about `0.56` peak target probability once readout decay is included. A future extension should test STIRAP-like timing, independent leg envelopes, and open-system optimal control to keep the protocol both dark-state-like and faster than the readout lifetime.
- **[P2 | HIGH]** Add photon-number-selective structure to the two-tone transfer if `n > 1` matters: The present simultaneous protocol treats `n = 2,3` only as diagnostics and shows strong chain leakage. A follow-up study should test number-selective detuning, sequential transfer, or shaped pulses that isolate one ladder segment at a time.

## Nice-to-Haves (P3)
- **[P3 | LOW]** Extend the practical reranking to optimal-control ansatz families: The present reranking is strong for the analytic waveform families requested here. A natural next extension is to test whether a short optimal-control pulse can recover selectivity under the same open-system conditions.

## Open Questions
- Does the best selective family remain Gaussian once a microscopic pump Hamiltonian and stronger Stark dressing are included?
- Can any short optimal-control pulse recover a selective readout-sideband primitive under the present linewidth and transmon-coherence scales?
- How much better would the storage selective regime become if the exact device-matched transmon `T1/T2` were known to be significantly longer than the matched local reference?
- Can a counter-intuitive two-tone sequence retain the low intermediate-state occupancy of the detuned Raman-like case while restoring a transfer time short compared with the readout lifetime?
- Is there a simple photon-number-dependent detuning rule that isolates the `n = 1` ladder without opening the longer `n > 1` conversion chain?

## What Was Tried and Did Not Work
- **Hard-gating the selective regime on projected swap fidelity**: Treating coherent SWAP quality as a hard eligibility requirement eliminated almost every selective candidate. The better interpretation is that the present analytic family set provides strong state-transfer pulses, while coherent-gate quality must be treated as a separate optimization target.
- **Assuming the sideband-reset example implicitly fixed transmon decoherence**: It does not. The example exports storage and readout noise but no transmon `T1/T2`. The study therefore had to add a separate, explicitly labeled transmon-reference sensitivity pass.
- **Using large common detuning to make the two-tone transfer dark-state-like with a constant square pulse**: In the single-photon case, increasing the common detuning does reduce the peak intermediate `|f>` population from about `0.50` to `0.12`, but it stretches the transfer from `29.5 ns` to `348 ns`. Under the measured readout linewidth, that slower protocol loses most of its practical advantage and falls to roughly `0.56` peak target probability in the noisy replay.
- **Applying the same constant simultaneous two-tone pulse logic to `n = 2,3`**: The reduced three-state picture breaks immediately because the drive couples into a longer ladder. The best resonant closed-system target probability then stays well below unity, so the protocol cannot be treated as a straightforward higher-photon generalization.

## Compute & Resource Notes
- Planned single-process execution because simulator import cost is large on this Windows machine.
- Main closed-system sweep plus finalist enrichment: `103.5 s`
- Extended run with transmon-noise scenarios and transmon-reference reranking: `534.8 s`
- Simultaneous two-tone extension scan plus open-system replays: `21.1 s`
- Validation pass after extension: `3.6 s`
- Notebook execution before the extension: `~6.6 s`
- Notebook execution after the two-tone extension: `~6.8 s`

## Resolved
- **[Resolved | 2026-04-10] Selective-versus-fast definitions made explicit**: The study uses machine-readable thresholds for selective, fast/unselective, and gate-oriented control.
- **[Resolved | 2026-04-10] Both storage and readout mode sidebands covered in one workflow**: The final study compares the two mode types side by side under the same device model, metric definitions, and figure/report pipeline.
- **[Resolved | 2026-04-10] Device-parameter provenance exported explicitly**: The study now writes a parameter table and device manifest that state exactly which local cQED parameters were used.
- **[Resolved | 2026-04-10] Transmon decoherence incorporated into the open-system study**: The workflow now replays the finalist pulses under mode-only noise and under matched local transmon-reference scenarios, and it reranks the shortlist under the transmon-reference model.
- **[Resolved | 2026-04-10] Gate-oriented objective added as a saved study artifact**: The study now exports gate-oriented winner tables and shows that the analytic family set remains unsuitable for phase-clean coherent SWAP control.
- **[Resolved | 2026-04-10] Selective readout-sideband impracticality promoted from suspicion to quantified result**: The report now states clearly that selective readout control fails under both the mode-only and transmon-inclusive noisy replays, and that no readout family remains threshold-valid under the matched transmon-reference reranking.
- **[Resolved | 2026-04-10] Simultaneous two-tone storage-to-readout transfer quantified for the local device**: The study now includes a dedicated extension showing that the single-photon resonant bright-state protocol is fast and strong, the detuned Raman-like variant is too slow once readout decay is included, and the constant-drive picture does not generalize cleanly to `n > 1`.
