# Echo Comparison Summary: multitone_sqr_arbitrary_fock_conditional_rotations

## What changed relative to the previous report

- The previous report tested only the single Gaussian multitone SQR ansatz.
- This extension adds an explicit echoed sequence with the time-ordered schedule: half SQR -> X_pi -> half SQR -> X_pi.
- The inserted X_pi pulses are finite Gaussian pulses of duration 40.0 ns about the x axis.
- Two fairness conventions are included: fixed total gate duration and fixed active SQR duration.

## Headline results

- Best single pulse: chi_plus_chiprime_na2_chiT1p0_familyC with fidelity 0.873625.
- Best echoed case: echoed_fixed_total on chi_only_na3_chiT5p0_familyD_seed317160 with fidelity 0.606132.
- Strongest random-target fidelity gain: echoed_fixed_total on chi_only_na3_chiT5p0_familyD_seed317160 with delta fidelity +0.428582.
- Strongest residual-Z reduction: echoed_fixed_total on chi_only_na4_chiT3p0_familyD_seed318160 with delta residual Z -0.635356 rad.

## Branch means

- single_pulse: mean fidelity 0.350448, best fidelity 0.873625, mean worst-block fidelity 0.459944, mean residual Z 0.952334 rad, mean transverse 1.460624 rad
- echoed_fixed_active: mean fidelity 0.224832, best fidelity 0.417515, mean worst-block fidelity 0.413314, mean residual Z 0.955247 rad, mean transverse 1.574112 rad
- echoed_fixed_total: mean fidelity 0.219042, best fidelity 0.606132, mean worst-block fidelity 0.421859, mean residual Z 0.971618 rad, mean transverse 1.547813 rad

## Echo minus single medians

- Family C, echoed_fixed_active: delta fidelity median -0.588711, delta worst-block median -0.189823, delta residual Z median +0.098722 rad, delta transverse median +0.301154 rad, improved-fidelity count 0/18, reduced-residual-Z count 2/18
- Family C, echoed_fixed_total: delta fidelity median -0.560351, delta worst-block median -0.158777, delta residual Z median +0.141102 rad, delta transverse median +0.247512 rad, improved-fidelity count 0/18, reduced-residual-Z count 0/18
- Family D, echoed_fixed_active: delta fidelity median +0.002700, delta worst-block median -0.002352, delta residual Z median -0.040052 rad, delta transverse median +0.111162 rad, improved-fidelity count 38/72, reduced-residual-Z count 38/72
- Family D, echoed_fixed_total: delta fidelity median -0.002212, delta worst-block median +0.000008, delta residual Z median +0.026359 rad, delta transverse median +0.010483 rad, improved-fidelity count 32/72, reduced-residual-Z count 31/72
