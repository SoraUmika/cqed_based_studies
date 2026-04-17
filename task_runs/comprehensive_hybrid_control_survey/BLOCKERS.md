# Blockers: Comprehensive Hybrid Control Survey

## Active Blockers

- **TargetStateMapping broken**: Returns objective=2.0 for all inputs. Blocks state-to-state GRAPE benchmarking. Workaround: use `make_target()` unitary + `UnitaryObjective`. Upstream fix required.
- **ConditionalDisplacement not a valid gateset**: Prevents ECD-style synthesis benchmarks. Upstream addition required.

## Resolved Blockers
- [x] **scipy version conflict**: scipy 1.17.1 broke cqed_sim import via scipy.signal. Fixed by downgrading to 1.14.1.
- [x] **Dimension mismatch (45 vs 30)**: Subspace.qubit_cavity_block assumes qubit_dim=2, incompatible with n_tr=3 models. Fixed by using n_tr=2 for synthesis tasks.
- [x] **GrapeResult API**: No `objective` attribute — used `objective_value` and `nominal_final_unitary` instead.
- [x] **Qobj to numpy conversion**: `np.array(Qobj)` returns 0-dim scalar. Fixed with `.full().flatten()`.
