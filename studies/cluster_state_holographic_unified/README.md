# Holographic Cluster-State Synthesis: Unified Study

## Problem Class
OPT, ANA

## Motivation

This unified study now contains two linked control problems built on the same structured gate families.

1. The baseline branch optimized the full logical transfer unitary `SWAP · CZ · (H ⊗ I)` on `{|g,0>, |g,1>, |e,0>, |e,1>}` within the exclusive families `D + R + SQR` and `D + R + CPSQR`.
2. The current follow-up replaces that full joint-unitary target with the restricted ground-sector transfer objective
   `U_joint(|g> ⊗ |psi>) ≈ |g> ⊗ H_c |psi>`
   on `{|g,0>, |g,1>}`, where `H_c` is the logical cavity Hadamard.

That objective change is not cosmetic. The old full-target winners achieve near-perfect fidelity on the joint logical map, but when they are rescored on the restricted ground-sector task they fail because they push most amplitude through the excited ancilla sector. The follow-up therefore introduces explicit ground-sector fidelity, ancilla-excitation, and support-leakage diagnostics, then reruns a targeted local restricted-target refinement for one representative per family and block depth.

## Goals

1. Implement a restricted ground-sector fidelity objective and leakage diagnostics in the shared study helpers.
2. Re-score the saved structured-family representatives under both the old full-target metric and the new ground-sector metric.
3. Warm-start a local restricted-target refinement at `N_cav = 12` for one representative per family and block count.
4. Determine, for each family, the best absolute fidelity, the minimum depth above `0.99`, and the most physically plausible candidate after ancilla, support, and truncation checks.
5. Regenerate the follow-up artifacts, figures, report conclusions, and reproducibility material.

## Methods

- Shared model, target, and gate builders: `scripts/common.py`
- Baseline exclusive-family search: `scripts/run_design_space_study.py`
- Ground-sector follow-up rescoring and local refinement: `scripts/run_ground_sector_followup.py`
- Primitive-level GRAPE diagnostics: `scripts/run_primitive_grape_diagnostics.py`
- Reproducibility notebook: `scripts/reproducibility_notebook.ipynb`

Retained baseline outputs:

- `data/corrected_scope_summary.json`
- `data/corrected_design_space_best.csv`
- `data/corrected_target_99_finalists.csv`
- `data/min_depth_summary.json`
- `artifacts/corrected_best_sqr.json`
- `artifacts/corrected_best_cpsqr.json`
- `data/grape_frontier_extension.json`
- `data/primitive_grape_diagnostics.json`

New ground-sector follow-up outputs:

- `data/ground_sector_followup_summary.json`
- `data/ground_sector_followup_candidates.csv`
- `artifacts/ground_sector_best_sqr.json`
- `artifacts/ground_sector_best_cpsqr.json`
- `figures/ground_sector_block_summary.pdf`
- `figures/ground_sector_objective_scatter.pdf`

Ground-sector follow-up workflow:

- Load the best saved full-target representative for each family and block count.
- Re-evaluate each saved representative at `N_cav = 10, 12, 14` under:
  - the original full-target fidelity,
  - the restricted ground-sector fidelity,
  - ancilla excitation,
  - cavity-support leakage.
- Warm-start a local restricted-target refinement at `N_cav = 12` for each block representative.
- Re-check the refined solutions at `N_cav = 10, 12, 14` and rank them by `N_cav = 12` fidelity, ancilla leakage, outside-target leakage, and truncation stability.

Ground-sector target convention:

- Ground-sector basis: `{|g,0>, |g,1>}`
- Cavity-only target: logical cavity Hadamard `H_c`
- Objective: `U_joint(|g> ⊗ |psi>) ≈ |g> ⊗ H_c |psi>`

## Expected Outcomes

- A direct numerical comparison between the old full-target winners and the new ground-sector winners.
- A family-by-family statement of minimum depth above `0.99` under the restricted objective.
- A physically grounded ranking that separates best absolute fidelity from the most plausible implementation.
- Updated local figures, artifacts, and summary tables for the follow-up branch.

## Known Limitations

- The ground-sector follow-up uses a targeted warm-start strategy rather than a brand-new exhaustive search. One representative per family and block count was refined locally at `N_cav = 12`; a fresh unrestricted search may still uncover different low-depth structures.
- The new rankings are still closed-system ideal-unitary results. No open-system process-fidelity study has yet been run for the ground-sector winners.
- Only one local seed and refinement budget (`maxiter = 16`) was used in the restricted-target pass. The winners are numerically clean, but local-minimum risk is not fully excluded.
- The follow-up validates truncation stability through `N_cav = 10, 12, 14`, ancilla excitation, and support leakage, but it does not yet include a dedicated noise-aware or pulse-level replay analysis for the ground-sector task.
- The separate GRAPE frontier remains a full-target reference benchmark, not a like-for-like ground-sector optimization baseline.

## Validation

- [x] Sanity checks — the old full-target winners were rescored under the restricted objective and shown to have poor ground-sector transfer because they strongly excite the ancilla.
- [x] Convergence — every restricted-target refined representative was checked at `N_cav = 10, 12, 14`; the retained winners have convergence deltas far below `1e-2`.
- [x] Literature / prior-study comparison — the original full-target structured study and the validated GRAPE frontier are retained as explicit baselines.

## Key Outcomes

| Category | Result |
|---|---|
| Objective shift | The old full-target winners are not good ground-sector transfers: their `N_cav = 12` restricted fidelities collapse to about `0.50`, with ancilla excitation close to `1`. |
| Best absolute `D + R + SQR` | 5 blocks, `RDS`, `n_active = 4`, levels `{0,1,2,3}`, `F_ground,12 = 0.9999260`, outside-target leakage worst `1.77e-4`. |
| Minimum-depth `D + R + SQR` above `0.99` | 3 blocks, `RDS`, `n_active = 4`, levels `{0,1,2,3}`, `F_ground,12 = 0.9997453`. |
| Most plausible `D + R + SQR` | Same as the best absolute SQR winner; ancilla excitation worst `4.13e-5`, support leakage worst `1.45e-4`, truncation delta `7.85e-10`. |
| Best absolute `D + R + CPSQR` | 4 blocks, `RDCP`, `n_active = 4`, levels `{0,1,2,3}`, `F_ground,12 = 0.9999284`, outside-target leakage worst `1.76e-4`. |
| Minimum-depth `D + R + CPSQR` above `0.99` | 3 blocks, `DRCP`, `n_active = 3`, levels `{0,1,2}`, `F_ground,12 = 0.9978771`. |
| Most plausible `D + R + CPSQR` | Same as the best absolute CPSQR winner; ancilla excitation worst `9.04e-5`, support leakage worst `8.56e-5`, truncation delta `1.19e-11`. |
| Two-block frontier | Both families approach but do not cross the `0.99` target at two blocks: SQR `0.9875652`, CPSQR `0.9886810`. |
| Cross-family best absolute | CPSQR 4-block `RDCP` is the strongest restricted-target solution in the current pass, by about `2.4e-6` over the SQR 5-block winner. |
| Cross-family minimum depth | Both families reach the `0.99` threshold at 3 blocks under the restricted target. |
| Full-target compatibility | The restricted-target winners generally have poor full-target fidelities (`~0.17` to `0.41`), confirming that the reduced objective is a genuinely different control problem. |
| Separate GRAPE reference | The saved GRAPE frontier remains a full-target benchmark only; its best validated point stays at replay `0.9760` and open-system process fidelity `0.9187` at `800 ns`. |

## Next-Stage Priority

1. Run a fresh restricted-target search at 2- and 3-block depth, rather than only warm-starting the best full-target representatives, to test whether either family can push the 2-block frontier above `0.99`.
2. Perform open-system validation for the 3-block and best-absolute ground-sector winners before treating the current ranking as hardware-relevant.
3. Expand ordering and level-subset refinements around the new shallow winners:
   - `D + R + SQR`: 3-block `RDS`
   - `D + R + CPSQR`: 3-block `DRCP` and 4-block `RDCP`
4. Compare the 3-block SQR and CPSQR winners under a common duration-aware or decoherence-aware metric to decide which family is the better practical compromise.
5. Keep the full-target and ground-sector branches separate in future reporting, because the winner sets are not interchangeable.

Current saved depth evidence from `data/ground_sector_followup_summary.json`:

- `D + R + SQR` first clears the `0.99` threshold at 3 blocks after restricted-target refinement.
- `D + R + CPSQR` also first clears the `0.99` threshold at 3 blocks after restricted-target refinement.
- Two-block representatives in both families remain sub-threshold but close, which makes that depth the next high-value search target.

## Reproduction

Full baseline design-space rerun:

```bash
cd studies/cluster_state_holographic_unified/scripts
python run_design_space_study.py
```

Ground-sector follow-up rerun:

```bash
cd studies/cluster_state_holographic_unified/scripts
python run_ground_sector_followup.py
```

Fast ground-sector rescoring without local refinement:

```bash
cd studies/cluster_state_holographic_unified/scripts
python run_ground_sector_followup.py --skip-refine
```

JAX-enabled GRAPE frontier rerun:

```bash
cd studies/cluster_state_holographic_unified/scripts
python run_grape_frontier_extension.py --durations-ns 500 600 800 1000 --engine auto
```

Primitive-level GRAPE diagnostic rerun:

```bash
cd studies/cluster_state_holographic_unified/scripts
python run_primitive_grape_diagnostics.py --engine auto --maxiter 200 --seeds 17 42 73
```

Notebook path:

```bash
cd studies/cluster_state_holographic_unified/scripts
jupyter notebook reproducibility_notebook.ipynb
```

## Status

COMPLETE
