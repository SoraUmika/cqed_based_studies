# Summary: ideal_sqr_direct_vs_echoed_multitone

## Audit Findings
- The earlier baseline arbitrary-rotation study optimized arbitrary SU(2) block targets rather than an ideal x-axis SQR gate.
- The earlier residual-Z follow-up used only one mid-sequence X_pi pulse, so it did not test the requested half-SQR -> pi -> half-SQR -> pi construction.
- A direct ideal SQR requires Theta_n = theta_n and Phi_n = 0 for every active manifold.
- A symmetric echoed sequence cancels first-order Z-type phase if the two halves accumulate the same Phi_n.
- Once the inserted X_pi pulses become manifold dependent, the clean toggling-frame cancellation is spoiled and the echoed protocol can fail even when the ideal first-order algebra is favorable.

## Best Overall
- direct_multitone on chi_only_smooth_x_na2_chiT5p0 with average gate fidelity 0.724539, mean residual-Z 0.026186 rad, mean transverse 1.225238 rad.

## Construction Means
- direct_multitone: mean fidelity 0.595142, best fidelity 0.724539, mean residual-Z 0.058152 rad, mean transverse 1.425383 rad.
- echoed_asymmetric: mean fidelity 0.248217, best fidelity 0.311341, mean residual-Z 0.318238 rad, mean transverse 1.300471 rad.
- echoed_independent: mean fidelity 0.248217, best fidelity 0.311341, mean residual-Z 0.318240 rad, mean transverse 1.300469 rad.
- echoed_symmetric: mean fidelity 0.248211, best fidelity 0.311322, mean residual-Z 0.318233 rad, mean transverse 1.300564 rad.

## Largest Echo Gain
- echoed_asymmetric on chi_plus_chiprime_staggered_x_na2_chiT3p0: delta fidelity -0.200881, delta residual-Z +0.399155 rad.

## Smallest Residual-Z Penalty
- echoed_asymmetric on chi_plus_chiprime_smooth_x_na3_chiT5p0: delta residual-Z +0.107005 rad, delta fidelity -0.403272.
