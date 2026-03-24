# Thermal Noise Cavity Sensing in a cQED System

## Problem Class

ANA | DES | REP

## Motivation

A post cavity in a cQED system can act as a **narrowband thermal noise sensor**:
when a target dissipative device (e.g., an attenuator) is coupled to the cavity, the
cavity acquires a nonzero photon occupation. By probing the cavity — directly or via
a dispersively coupled ancilla qubit — one can in principle infer the effective bath
temperature / occupation of the target device.

This study provides a first-principles theoretical and numerical investigation of
this sensing scheme, answering: under what conditions is the target device's thermal
occupation identifiable, and what are the fundamental limitations?

## Goals

1. Build and validate a finite-temperature cavity model (multiple thermal baths,
   Lindblad master equation) against known analytic results.
2. Map out sensing feasibility: how large must the target coupling fraction be for
   the cavity occupation to shift measurably above the background floor?
3. Simulate ancilla-based measurement protocols (number-selective spectroscopy,
   Ramsey) and assess bias from qubit imperfections.
4. Solve the inverse problem: given simulated observables, what bath parameters
   can actually be recovered, and are there fundamental degeneracies?

## Methods

### cqed_sim usage and documented gaps

**cqed_sim is used for:**
- Experimental parameter reference: `DispersiveTransmonCavityModel` constants (ω_c,
  ω_q, χ, T1, T2) match the existing `sqr_pulse_waveform_design` study.
- Physical frame/convention alignment with the existing codebase.

**Documented gaps requiring QuTiP directly:**

| Gap | Reason |
|-----|--------|
| Cavity-only model (no transmon required) | `DispersiveTransmonCavityModel` always includes a transmon; cavity-only physics studied with `qt.destroy(N)` directly |
| Multiple thermal baths per mode | `NoiseSpec` supports a single `(kappa, nth)` per mode; multi-bath decomposition uses custom `c_ops` list |
| Steady-state solving | cqed_sim is time-domain only; uses `qt.steadystate()` for exact steady states |
| Dispersive ancilla Hamiltonian without pulses | Phase 3 uses a minimal `qt.tensor` dispersive Hamiltonian outside the pulse infrastructure |

All gaps are scientifically motivated. The multi-bath model is the central physics of
the study. Suggested upstreaming: a `MultiThermalBath` spec and `qt.steadystate`
wrapper in `cqed_sim.sim`.

### Key physics

**Master equation** (rotating frame, zero drive):
```
ρ̇ = Σ_j κ_j(n_j+1) D[a]ρ + Σ_j κ_j n_j D[a†]ρ
```

**Steady-state mean photon number:**
```
n̄_ss = (Σ_j κ_j n_j) / (Σ_j κ_j)
```

**Transient occupation:**
```
n̄(t) = n̄_ss + (n̄(0) − n̄_ss) exp(−κ_tot t)
```

**Thermal photon distribution:**
```
P_n = n̄^n / (1 + n̄)^{n+1}
```

**Dispersive ancilla Hamiltonian (rotating frame):**
```
H/ℏ = χ a†a |e⟩⟨e|    (χ < 0 in this system)
```

## Expected Outcomes

- Numerical steady-state and transient results agree with analytic formulas to < 10⁻⁶.
- Thermal P_n distribution matches analytic formula to < 10⁻⁵ (KL divergence).
- Sensing feasibility: κ_target/κ_tot ≥ 0.1 needed for n_target ≥ 0.1 to produce
  a ≥ 0.01 photon shift above background.
- Fundamental degeneracy: n_target and κ_target cannot be separately inferred from
  cavity photon number alone; κ_tot from transient rate breaks partial degeneracy.
- Ancilla spectroscopy recovers P_n accurately for χ T2_q >> 1 (satisfied here:
  χ T2_q ≈ 2π × 56.8).

## Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| ω_c / (2π) | 5.241 GHz | Post cavity frequency |
| ω_q / (2π) | 6.197 GHz | Ancilla qubit frequency |
| \|χ\| / (2π) | 2.84 MHz | Dispersive coupling |
| T₁ (qubit) | 20 μs | Qubit relaxation |
| T₂ (qubit) | 20 μs | Qubit coherence |
| κ_tot / (2π) | 100 kHz | Total cavity decay rate |
| κ_target | 0.0–1.0 × κ_tot | Swept |
| n_target | 0–10 | Target bath occupation (swept) |
| n_bg | 0.01 | Background bath (cold) |
| n_int | 0.0 | Internal loss (zero-T) |
| N (Fock) | 30 | Default cavity truncation |

## Status

COMPLETE — all 51 validation tests pass; report compiled to `report/report.pdf` (11 pages)

## Suggested Upstreaming

- `cqed_sim.sim.noise.MultiThermalBath` dataclass accepting a list of `(kappa, nth)` pairs
- `cqed_sim.sim.steadystate(model, noise)` wrapping `qt.steadystate()`
- `cqed_sim.measurement.spectroscopy_signal(pn, chi, gamma_q, omega_probe)` helper
