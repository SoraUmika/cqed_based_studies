# Native / Rich Multitone Feasibility for Ideal SQR and CPSQR

## Problem Class
OPT | REP | DES | ANA

## Motivation
This study is the final patched-package follow-up to the multitone SQR line of work. Earlier studies answered neighboring questions about arbitrary block-diagonal targets, residual-Z cancellation, direct-vs-echoed ideal-SQR constructions, and multi-input validation, but they did not yet give one unified final answer to the narrow question asked here:

1. Can a native direct multitone waveform realize an ideal x-axis SQR?
2. If not, can a relaxed conditional-phase selective qubit rotation (CPSQR) be realized?
3. Does a composite or echoed multitone construction materially change that answer?
4. Does any apparent success survive when the full joint qubit-cavity action is enforced rather than only a reduced qubit-only criterion?

The study uses the current patched `cqed_sim` only and re-audits all patch-sensitive conventions before trusting earlier helper code or saved metrics.

## Goals
1. Audit the patched `cqed_sim` conventions relevant to SQR/CPSQR and record any stale study-local helper mismatches that must not be reused.
2. Define fixed strict-SQR and relaxed-CPSQR target families and keep them unchanged throughout the study.
3. Implement a unified diagnostic stack that separates:
   - reduced qubit-only success,
   - full-state success,
   - full joint-unitary success.
4. Screen direct native and richer multitone waveform families on a representative structured grid.
5. Test composite and echoed constructions under matched fairness conventions.
6. Compare `chi_only` and `chi_plus_chiprime`.
7. Estimate the minimum gate duration at which strict SQR or relaxed CPSQR reaches convincing thresholds, if such thresholds are reached at all.
8. Produce a final report, notebook, figures, and machine-readable artifacts with a scientifically honest conclusion.

## Methods
- Simulation framework: patched `cqed_sim` only.
- Convention audit sources:
  - `cqed_sim.core.ideal_gates`
  - `cqed_sim.core.conventions`
  - `cqed_sim.core.frequencies`
  - `cqed_sim.calibration.conditioned_multitone`
  - `cqed_sim.calibration.sqr`
  - `cqed_sim.calibration.targeted_subspace_multitone`
  - `cqed_sim.pulses.calibration`
  - `cqed_sim.sequence`
  - `physics_and_conventions/physics_conventions_report.tex`
  - selected upstream regression tests for tensor ordering, Gaussian IQ convention, and additive SQR amplitude correction
- Analytic preliminary:
  - derive the direct multitone first-order manifold-resolved drive matrix and identify residual spectator-induced `Z` structure;
  - derive the ideal echoed cancellation argument for instantaneous `X_pi` refocusing pulses and the finite-pulse obstruction;
  - define CPSQR as strict ideal x-axis SQR modulo a per-manifold cavity block phase and a post-rotation qubit `Z` phase:
    `U_CPSQR = sum_n exp(i gamma_n) |n><n| tensor R_z(delta_n) R_x(theta_n)`.
- Numerical workflow:
  - family-screening stage on the harder patched `chi_plus_chiprime` model for `N_active = 2, 3`, `|chi| T / 2pi = 1, 3, 5`, and structured targets `smooth_x`, `staggered_x`;
  - selected-family comparison stage on `chi_only` and `chi_plus_chiprime`, including `N_active = 4` and a small random `random_x` ensemble;
  - duration-refinement stage for the most promising direct and echoed families.
- Local study code is used only for:
  - corrected run-config construction that avoids stale `fock_fqs_hz` wrappers,
  - CPSQR fitting and metric extraction,
  - richer sampled-envelope family optimization,
  - report/notebook generation.

## Patched Convention Audit
- `cqed_sim.core.ideal_gates.sqr_op` and the tensor-product regression tests confirm that the runtime tensor order is `qubit tensor cavity`.
- The full-space flat computational basis is qubit-major, while the logical addressed-subspace ordering used throughout this study is blockwise:
  `( |g,0>, |e,0>, |g,1>, |e,1>, ... )`.
- `phi_n` in the conditioned-multitone helpers is the intended in-plane rotation-axis parameter used by the XY-rotation helper, not the final Bloch-sphere azimuth extracted after evolution.
- The patched carrier convention is generated from the negative manifold transition frequency in the chosen rotating frame.
- The patched `cqed_sim.pulses.calibration.build_sqr_tone_specs(...)` path expects `fock_fqs_hz` to be absolute manifold transition frequencies when an override is supplied. Earlier study-local wrappers that passed frame-shifted frequencies are stale and must not be reused.
- Reduced effective-qubit diagnostics in this study are extracted from the full propagator by projecting each addressed manifold to a qubit channel and then fitting either the strict `R_x(theta_n)` block or the relaxed CPSQR block.
- Full-joint diagnostics use the restricted addressed-subspace propagator directly and therefore retain inter-manifold phase relations that the reduced metrics intentionally discard.

## Analytic Preliminary
### Direct Multitone Feasibility
Within each addressed Fock manifold `n`, the driven effective qubit Hamiltonian takes the form
`H_n^(eff)(t) = (1/2) [ Omega_x,n(t) X + Omega_y,n(t) Y + Delta_n^(ac)(t) Z ]`
after moving into the rotating frame and collecting the manifold-resolved resonant terms. The desired ideal-SQR limit requires the time-integrated transverse term to reproduce `theta_n` about a common `x` axis while the residual coherent `Z` accumulation
`Phi_n = int dt Delta_n^(ac)(t)`
is simultaneously suppressed for every addressed `n`. A direct multitone family therefore has to fit one transverse target angle per manifold while also canceling or equalizing manifold-dependent Stark/spectator phases. This already suggests why strict SQR is hard but CPSQR can remain viable: if the dominant mismatch is `Phi_n`, then allowing a manifold-dependent post-rotation `Z` phase relaxes the hardest constraint.

### Echo Argument
For the idealized schedule
`half-SQR -> X_pi -> half-SQR -> X_pi`,
instantaneous manifold-independent `X_pi` refocusing pulses leave the desired `X` rotation invariant while toggling `Z -> -Z`, so the first-order residual `Z` accumulation cancels when the two driven halves are symmetric. The obstruction in the real multitone setting is that the inserted `X_pi` pulse is itself implemented on a dispersive manifold ladder: it is not perfectly manifold independent, and its own coherent error need not commute with the target. Once the refocusing pulse becomes manifold dependent, the clean toggling-frame cancellation argument no longer guarantees strict ideal-SQR success.

### CPSQR Interpretation
The relaxed family fixed for this study is
`U_CPSQR = sum_n exp(i gamma_n) |n><n| tensor R_z(delta_n) R_x(theta_n)`.
This keeps the target rotation angle `theta_n` fixed on each addressed manifold but allows the dominant residual error channel to be a manifold-dependent cavity-block phase `gamma_n` together with a post-rotation qubit `Z` phase `delta_n`. Strict ideal SQR corresponds to the special case `gamma_n = delta_n = 0` for all addressed manifolds.

### Input-Consistency Logic
Matching only `|g,n>` can hide coherent axis errors because population transfer alone is insensitive to many relative phases. A two-state check still misses general Bloch-sphere shear and conditional-phase errors. The quartet
`{ |g,n>, |e,n>, |+x,n>, |+y,n> }`
is therefore used as the practical reduced-state analogue of qubit process tomography on each addressed manifold.

## Assumptions and Convergence Criteria
- Main Hamiltonian models: `chi_only` and `chi_plus_chiprime`; cavity Kerr and `n_tr > 2` are reserved for follow-up sensitivity checks only if the main closed-system answer warrants them.
- Addressed cavity window: contiguous `n = 0, ..., N_active - 1` with `N_active in {2, 3, 4}` across the staged grid.
- Time discretization: fixed simulation time step from the shared study constants, with duration sweeps reported in the dimensionless form `|chi| T / 2pi`.
- Strict success thresholds emphasized in the duration study:
  - reduced strict process fidelity `>= 0.99` and `>= 0.995`;
  - reduced CPSQR process fidelity `>= 0.99`;
  - full addressed-subspace strict joint process fidelity `>= 0.99`;
  - full addressed-subspace CPSQR joint process fidelity `>= 0.99`.
- A convincing strict-SQR claim requires agreement across reduced quartet validation, full-state quartet validation, and the full addressed-subspace joint operator.

## Suggested Upstreaming
- Add an explicit public helper or validation guard for the `fock_fqs_hz` absolute-frequency override so older study-local frame-shifted wrappers cannot silently produce incorrect carriers.
- Expose a package-level reduced-channel extraction helper for manifold-resolved process diagnostics to reduce duplicated study-local metric code.

## Expected Outcomes
- Direct native multitone should remain the clean baseline for strict ideal-SQR testing, but previous studies suggest that coherent manifold-dependent `Z` structure will limit strict success even when reduced qubit action looks promising.
- Richer direct envelopes may improve reduced metrics without fully preserving the joint cavity block structure.
- Echo may help primarily by refocusing first-order `Z` accumulation, but any benefit may weaken once finite manifold-dependent `X_pi` pulses are included.
- CPSQR is expected to be more achievable than strict ideal SQR if the dominant remaining error is conditional `Z` structure rather than generic transverse mismatch.

## Known Limitations
- The main workflow is closed-system and uses the two-level qubit model (`n_tr = 2`) for the principal grid.
- The family screen is intentionally staged rather than brute-force over every family on every case.
- Duration thresholds are finite-grid estimates refined only for the most promising families.
- Any optional cavity-Kerr or higher-transmon-level follow-up is secondary to the main closed-system answer and may remain future work.

## Final Outcome
1. **Can direct native/rich multitone realize ideal SQR?**
   - Yes, but only conditionally. The best direct family (`reduced_unitary_direct`) achieved full joint ideal-SQR success on the easier addressed windows:
     - `chi_only`, `N_active = 2`: strict joint fidelity `0.9988` on `smooth_x` at `|chi|T/2pi = 5`.
     - `chi_plus_chiprime`, `N_active = 2`: essentially the same best strict joint fidelity (`0.9988`) on `smooth_x` at `|chi|T/2pi = 5`.
     - `chi_only`, `N_active = 3`: strict joint fidelity remained above `0.99` on the refined smooth-target duration grid and on the structured comparison grid.
   - The same claim does **not** generalize across the full tested space. With `chi_plus_chiprime` and `N_active = 3`, the best direct strict joint fidelity on the refined smooth-target duration grid was only about `0.9846`, and with `chi_plus_chiprime` and `N_active = 4` the best direct strict joint fidelities on the comparison grid stayed below `0.95`.
2. **Can CPSQR be realized?**
   - Yes, convincingly. The best echoed family (`echoed_independent`) reached essentially unit CPSQR joint fidelity on many structured cases, including `chi_plus_chiprime`, `N_active = 2`, `staggered_x`, `|chi|T/2pi = 5`, where the CPSQR joint fidelity was `0.9999997+` while the strict joint fidelity on the same pulse remained only about `0.868`.
3. **Does echo help?**
   - Yes for the relaxed CPSQR target, but not as a general strict-SQR fix. Echoed families repeatedly converted poor or moderate strict joint fidelity into near-unit CPSQR fidelity, which is strong evidence that the dominant residual mismatch is conditional phase structure rather than an arbitrary transverse control failure.
4. **Reduced-only versus full-joint success**
   - The reduced manifold-resolved qubit picture is consistently more optimistic than the full joint-unitary picture.
   - Echoed families often look almost perfect on reduced strict/CPSQR quartet metrics while still failing the strict full-joint metric by a wide margin.
   - The strict-SQR question therefore has different answers in the two problem versions:
     - ignoring cavity subspace: many more cases look successful;
     - enforcing full joint qubit-cavity action: strict ideal SQR is only supported on a limited subset of the tested grid.
5. **Minimum-duration findings**
   - Best direct strict-SQR family (`reduced_unitary_direct`):
     - `chi_only`, `N_active = 2`: first reached strict joint fidelity `>= 0.99` at `|chi|T/2pi = 3.0`.
     - `chi_plus_chiprime`, `N_active = 2`: also first reached strict joint fidelity `>= 0.99` at `|chi|T/2pi = 3.0`.
     - `chi_only`, `N_active = 3`: first reached strict joint fidelity `>= 0.99` at `|chi|T/2pi = 3.0`.
     - `chi_plus_chiprime`, `N_active = 3`: never reached strict joint fidelity `>= 0.99` on the refined smooth-target duration grid.
   - Best echoed CPSQR family (`echoed_independent`):
     - `chi_only`, `N_active = 2`: first reached CPSQR joint fidelity `>= 0.99` at `|chi|T/2pi = 1.5`.
     - `chi_plus_chiprime`, `N_active = 2`: also first reached CPSQR joint fidelity `>= 0.99` at `|chi|T/2pi = 1.5`.
     - `chi_only`, `N_active = 3`: first reached CPSQR joint fidelity `>= 0.99` at `|chi|T/2pi = 2.0`.
     - `chi_plus_chiprime`, `N_active = 3`: also first reached CPSQR joint fidelity `>= 0.99` at `|chi|T/2pi = 2.0`.
   - No echoed duration sweep reached strict joint fidelity `>= 0.99`; the echoed benefit is therefore a CPSQR benefit, not a strict-SQR duration win.
6. **Overall best current answer**
   - The current best strict ideal-SQR family is `reduced_unitary_direct`.
   - The current best CPSQR family is `echoed_independent`.
   - The strongest scientifically honest conclusion is:
     - strict ideal SQR is achievable only on a restricted subset of the tested direct multitone problem;
     - CPSQR is broadly easier and is the more faithful description of what echoed multitone waveforms can realize;
     - once the full joint cavity action is enforced, broad positive claims about ideal SQR do not survive the entire tested grid.

## Status
COMPLETE
