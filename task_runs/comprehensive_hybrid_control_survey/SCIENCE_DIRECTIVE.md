# Science Directive: Comprehensive Hybrid Control Survey — Iteration 1

## Physics Context
In dispersive cQED, the transmon (ancilla) and cavity (storage) interact via the
dispersive coupling chi. The effective Hamiltonian (rotating frame) is:

H_0/hbar = (alpha/2) b^dag^2 b^2 + chi * n_c * n_q + (K/2) * n_c * (n_c - 1)

Control is exerted through drives on the qubit and cavity channels.

The key insight: different control paradigms exploit different terms in this Hamiltonian
and are therefore optimal in different physical regimes.

## Control Paradigms to Benchmark

### 1. Native Constructive Gate Set (SQR + Displacement + Rotation)
**Hamiltonian resource:** chi-dependent photon-number-selective qubit transitions.
Each Fock manifold |n> has a distinct qubit transition frequency omega_ge(n) = omega_ge(0) + n*chi.
When |chi| >> 1/T_gate, individual manifolds can be addressed selectively.

**Primitive gates:**
- Displacement D(alpha): cavity phase-space translation
- Qubit rotation R(theta, phi): unconditional qubit rotation
- SQR(thetas, phis): number-selective qubit rotations (multitone pulse addressing all manifolds simultaneously)
- SNAP(phases): number-selective phase gate (realized via SQR with specific angles + displacements)

**Expected regime:** Strong dispersive (|chi/2pi| > 1 MHz)

### 2. GRAPE Optimal Control
**Hamiltonian resource:** Full time-dependent control of qubit + cavity drives simultaneously.
No assumption about dispersive selectivity — directly optimize the piecewise-constant
control Hamiltonian.

**Expected regime:** Any, but most impactful where constructive decomposition fails
or requires too many gates.

### 3. ECD-Style Control (Echoed Conditional Displacement)
**Hamiltonian resource:** Conditional displacement (qubit-state-dependent cavity translation)
interspersed with qubit rotations. The key operation is:
CD(beta) |g>|psi_c> -> |g> D(+beta/2)|psi_c>, |e>|psi_c> -> |e> D(-beta/2)|psi_c>

**Expected regime:** Weak-to-intermediate dispersive, where large cavity displacements
are acceptable and individual Fock manifold selectivity is poor.

### 4. Hybrid: Constructive Decomposition + GRAPE Refinement
Decompose target into native gates, then use GRAPE to optimize the slowest/hardest
gates or the full sequence.

**Expected regime:** Large logical subspaces where pure constructive depth is high.

## Benchmark Tasks (Concrete)

### Task A: Ancilla X_pi gate with cavity spectator
- Target: X_pi on qubit, cavity in |0> (trivial) and |2> (non-trivial due to chi)
- Methods: Gaussian pulse, DRAG pulse, GRAPE
- Metrics: fidelity, leakage to |f>, gate time
- Sweep: chi from 0.1 MHz to 5 MHz

### Task B: Fock state preparation |0> -> |n>
- Target: prepare |g,n> from |g,0> for n = 1, 2, 3
- Methods: SQR+displacement decomposition, GRAPE, ECD-style
- Metrics: state fidelity, total duration, max cavity population, leakage

### Task C: Conditional phase gate
- Target: diag(1, 1, ..., 1, e^{i*phi}, ...) on qubit-cavity subspace — phase on |e> sector
- Methods: Free dispersive evolution, SQR, GRAPE
- Metrics: process fidelity in logical subspace, leakage, gate time

### Task D: SNAP gate (cavity-only phase map)
- Target: apply specific phases to first N Fock levels
- Methods: SQR-based SNAP decomposition, GRAPE
- Metrics: fidelity, total duration, number of tones

### Task E: Scaling test — Fock |0> -> |5>
- Target: prepare |g,5> from |g,0>
- Compare constructive depth vs GRAPE for increasing target Fock number
- Metrics: fidelity, gate count, total time

## Parameter Regimes

### Baseline system parameters (from AGENTS.md defaults)
- omega_q/2pi = 6.150 GHz
- omega_c/2pi = 5.241 GHz
- alpha/2pi = -255 MHz
- chi/2pi = -2.84 MHz (strong dispersive baseline)
- K/2pi = -28 kHz
- n_cav = 15 (truncation), n_tr = 3

### Chi sweep values
- Weak: chi/2pi = -0.1, -0.3 MHz
- Intermediate: chi/2pi = -1.0, -2.0 MHz
- Strong: chi/2pi = -2.84, -5.0 MHz

### Coherence (for open-system benchmarks)
- T1 = 10 us (qubit)
- T_phi = 15 us (qubit pure dephasing)
- kappa = 1/(500 us) (cavity decay, i.e. cavity T_cav = 500 us)

### Control budgets
- Max qubit drive: 50 MHz (amplitude in rad/s units: 2pi * 50e6)
- Max cavity drive: 20 MHz
- Pulse duration range: 20 ns to 5000 ns
- dt = 2 ns (compilation time step)

## Ordered Action Items

1. Create `scripts/benchmark_common.py` — shared model setup, metrics functions, plotting utilities
2. Create `scripts/benchmark_ancilla_control.py` — Task A (Gaussian, DRAG, GRAPE for qubit gates)
3. Create `scripts/benchmark_fock_prep.py` — Task B (SQR+D, GRAPE, ECD for Fock state prep)
4. Create `scripts/benchmark_conditional_phase.py` — Task C (dispersive, SQR, GRAPE)
5. Create `scripts/benchmark_snap_gate.py` — Task D (SQR-SNAP, GRAPE)
6. Create `scripts/benchmark_scaling.py` — Task E (scaling test)
7. Create `scripts/chi_sweep.py` — chi parameter sweep across methods
8. Create `scripts/open_system_comparison.py` — open-system benchmarks
9. Create `scripts/validation_convergence.py` — truncation convergence and sanity checks
10. Create `scripts/generate_figures.py` — combined figure generation
11. Collect all results, write the survey + numerical report

## Success Criteria
- At least 4 control methods benchmarked on at least 3 tasks each
- Chi sweep covers at minimum 3 distinct regimes
- Open-system simulation for at least the top 2 methods
- Convergence check for truncation dimension
- Clear recommendation matrix produced
