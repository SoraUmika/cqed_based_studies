# Science Directive — Iteration 1

## Study Objective
Determine, with the smallest possible analytic and numerical workload, whether free dispersive idle evolution in a qubit-storage system can produce a $\pi$ conditional phase difference between the $n=0$ and $n=1$ cavity manifolds over a finite time, and prepare the minimal next-step implementation needed to validate the staged research loop.

## Problem Classification
ANA. This is a system-analysis study: the task is to characterize a fixed physical mechanism and its finite-time phase consequence, not to optimize a control or reproduce a large published dataset.

## Physics Context
The minimal model is a dispersive qubit-storage Hamiltonian in the rotating frame,
$$
H/\hbar = \chi \hat{n} |e\rangle\langle e|,
$$
where $\hat{n}$ is the cavity photon-number operator and $\chi$ is the storage-qubit dispersive shift. In the excited-state branch, cavity manifold $n$ acquires phase at rate $\chi n$, so the central question reduces to whether the $n=1$ and $n=0$ manifolds accumulate a relative phase of $\pi$ after a finite idle time.

## Analytic Preliminary
From the first-principles dispersive Hamiltonian above, the manifold-resolved excited-branch phase is
$$
\phi_n(t) = -\chi n t.
$$
Therefore,
$$
\Delta \phi_{1,0}(t) = \phi_1(t) - \phi_0(t) = -\chi t \pmod{2\pi}.
$$
The first finite idle time at which the relative phase reaches $\pi$ is
$$
t_{\pi} = \frac{\pi}{|\chi|}, \qquad \chi \neq 0.
$$
So the idealized answer is yes: free dispersive evolution alone is sufficient whenever the dispersive shift is nonzero.

Controlled approximations and validity conditions:
- First-order dispersive Hamiltonian only.
- Closed-system idle evolution with no drives and no decoherence.
- Focus on $n=0$ and $n=1$, so the self-Kerr term $K n(n-1)/2$ vanishes identically for the target manifolds.
- The dispersive approximation remains valid only when $|g/\Delta| \ll 1$ and higher-order corrections remain perturbative over the chosen occupations and idle times.

## Hypotheses
1. If $\chi \neq 0$, free dispersive evolution alone produces a $\pi$ conditional phase difference at $t_{\pi} = \pi/|\chi|$.
2. A tiny cqed_sim confirmation using `DispersiveTransmonCavityModel` and `dispersive_phase` will reproduce the analytic phase within wrapped-phase numerical precision.
3. For the representative shift $\chi = -2\pi \times 2.84\,\mathrm{MHz}$, the first $\pi$ crossing will occur near $176\,\mathrm{ns}$, keeping the numerical workload trivial.

## Experiment Design
### Experiment 1: Analytic idle-phase derivation
- **Purpose:** Establish the finite-time $\pi$ condition before any numerics.
- **Method:** Derive the manifold-resolved phase accumulation under the minimal dispersive Hamiltonian.
- **Parameters:** Symbolic $\chi$, then representative $\chi = -2\pi \times 2.84\,\mathrm{MHz}$.
- **Expected outcome:** $t_{\pi} = \pi/|\chi|$ for nonzero $\chi$.
- **Success criterion:** The derivation makes the sign convention explicit and yields a unique first positive $t_{\pi}$.

### Experiment 2: cqed_sim unitary confirmation
- **Purpose:** Verify the analytic result with the lightest framework-backed computation.
- **Method:** Use `cqed_sim.gates.coupled.dispersive_phase` with qubit-first ordering and cavity dimension at least 2; use `DispersiveTransmonCavityModel` only as needed for sign-convention or energy cross-checks.
- **Parameters:** `cavity_dim = 2` or `3`, `qubit_dim = 2`, representative $\chi = -2\pi \times 2.84\,\mathrm{MHz}$, and a short idle-time scan around $t_{\pi}$.
- **Expected outcome:** The relative phase between the $n=1$ and $n=0$ manifolds crosses $\pi$ near $176\,\mathrm{ns}$.
- **Success criterion:** The wrapped analytic and numerical phases agree to within floating-point precision near $t_{\pi}$.

### Experiment 3: Minimal output generation
- **Purpose:** Produce the minimum artifacts needed for the loop.
- **Method:** Save one JSON artifact with parameters, analytic $t_{\pi}$, and sampled numerical phases; generate one simple phase-versus-idle-time figure with the $\pi$ crossing marked.
- **Parameters:** Roughly 25-51 time samples over $[0, 2 t_{\pi}]$.
- **Expected outcome:** One artifact and one figure, both lightweight.
- **Success criterion:** Outputs are saved under the study directories and total runtime stays far below two minutes.

## Execution Plan
1. **[IMPLEMENT]** Write a single lightweight script that evaluates the analytic $t_{\pi}$ and the cqed_sim dispersive unitary.
   - Files to create: `studies/auto_workflow_probe/scripts/free_dispersive_pi_probe.py`
   - Expected output: one JSON artifact and one PNG/PDF figure
2. **[RUN]** Execute the tiny idle-time scan at the representative dispersive shift.
   - Script: `studies/auto_workflow_probe/scripts/free_dispersive_pi_probe.py`
   - Expected runtime: under 10 s on CPU
3. **[ANALYZE]** Confirm the $\pi$ crossing and annotate the saved outputs.
   - Generate figures: `phase_difference_vs_idle_time.png` and `phase_difference_vs_idle_time.pdf`
4. **[VALIDATE]** Compare analytic and numerical phases and repeat with one minimal Hilbert-space variation.
5. **[DOCUMENT]** Update `EXECUTION_SUMMARY.md`, `TASK_CHECKLIST.md`, `PROGRESS_LOG.md`, and `study_state.json` after implementation, without entering report generation yet.

## Assumptions and Approximations
- First-order dispersive Hamiltonian only; no $\chi'$ or higher-order Kerr corrections in the minimal run.
- Closed-system free evolution with no decoherence channels.
- Internally consistent rad/s and seconds units.
- Tiny Hilbert-space truncation is sufficient because only the $n=0$ and $n=1$ manifolds are probed.

## Known Risks
- A sign-convention mismatch between the analytic notes and cqed_sim helper conventions (`n_e` versus `z`) could create an apparent minus sign unless it is recorded explicitly during implementation.
- The skill file's configured API-reference path was stale on this machine; future invocations should continue using the discovered local Box path unless the instruction is updated.
- Because this probe is intentionally lightweight, a later reviewer may request a full idle-sequence cross-check if the helper-based confirmation is judged too abstract.

## Stopping Criteria for This Iteration
- `SCIENCE_DIRECTIVE.md` is written.
- The README contains the mandatory analytic-first planning sections.
- The study state is updated to a planned status with explicit pending implementation tasks.
- No implementation, validation, or reporting work is performed in this invocation.

## Compute Budget Estimate
- Analytic work: negligible.
- Numerical confirmation: well under 1 s for the core unitary evaluation and under 10 s end-to-end including file writes.
- Total budget for the next phase: comfortably below the user's two-minute CPU limit.