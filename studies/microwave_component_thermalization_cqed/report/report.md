# Thermalization of Microwave Components in cQED

## Problem statement
This study asks how a hot microwave component in a superconducting cQED stack maps onto measurable qubit and cavity observables, which observable is the best thermometer, when thermal loading becomes dangerous for coherence or readout, and what the next attenuator, cable, Teflon, and multimode experiments should measure first.

The model hierarchy is intentionally skeptical. We separate:

- an exact dispersive model, which captures thermal cavity occupation and thermal-photon dephasing but cannot heat the qubit population by itself;
- a weakly dressed qubit-readout model, which adds a small residual exchange channel so that qubit excitation can become a thermometer;
- a reduced multimode model, where a hot auxiliary mode stands in for a cable mode or thermally imperfect resonant structure coupled to the local readout mode;
- a transient step model, which captures only the intrinsic quantum response after the effective bath occupation changes, not the slower macroscopic thermalization of the hardware component itself.

## Model definitions
For the single-mode dispersive study we used a three-level transmon and one readout mode with

\[
H_{\mathrm{disp}}/\hbar = \omega_q b^\dagger b + \frac{\alpha}{2} b^\dagger b^\dagger b b + \omega_r a^\dagger a + \chi\, a^\dagger a\, b^\dagger b.
\]

The bath occupation is

\[
\bar n_{\mathrm{th}}(T,\omega_r)=\frac{1}{e^{\hbar \omega_r/(k_B T)}-1},
\]

with Lindblad cavity loss

\[
\kappa (\bar n_{\mathrm{th}}+1)\,\mathcal D[a]\rho + \kappa \bar n_{\mathrm{th}}\,\mathcal D[a^\dagger]\rho.
\]

For the weakly dressed thermometer model we replaced the strict dispersive coupling with a small exchange interaction

\[
H_{\mathrm{dress}}/\hbar = H_q + \omega_r a^\dagger a + g(a^\dagger b + a b^\dagger),
\]

so that thermal cavity photons can admix into the qubit sector and generate a nonzero excited-state population. This is the minimal extension needed to make qubit-population thermometry physically meaningful.

For the multimode study we used one local readout mode and one hot auxiliary mode with

\[
H_{\mathrm{multi}}/\hbar = H_q + \omega_r a_r^\dagger a_r + \omega_s a_s^\dagger a_s + \chi_r a_r^\dagger a_r\, b^\dagger b + \chi_s a_s^\dagger a_s\, b^\dagger b + J(a_s^\dagger a_r + a_s a_r^\dagger),
\]

plus a dressed variant where the qubit couples weakly to the readout mode through a residual exchange channel. The auxiliary mode is the only hot bath in that model; the qubit environment and local readout bath remain cold.

The baseline parameters were representative rather than exact slide values: \(\omega_q/2\pi=\SI{5.2}{GHz}\), \(\omega_r/2\pi=\SI{7.0}{GHz}\), \(\alpha/2\pi=\SI{-220}{MHz}\), \(\chi/2\pi=\SI{-0.25}{MHz}\), \(g/2\pi=\SI{30}{MHz}\), and \(\kappa_r^{-1}=\SI{150}{ns}\). Missing project-specific values were treated as sweep axes rather than hidden assumptions.

## Main results

### A. The cavity is the cleanest thermometer; the qubit is a nonlinear secondary thermometer
The exact dispersive steady-state cavity occupation tracked the Bose-Einstein bath occupation almost one-to-one across the full temperature sweep. Representative values were:

- \(\SI{50}{mK}\): \(\bar n_{\mathrm{th}}\approx 1.2\times 10^{-3}\), simulated cavity occupation \(n_r\approx 1.2\times 10^{-3}\)
- \(\SI{100}{mK}\): \(\bar n_{\mathrm{th}}\approx 0.036\), simulated cavity occupation \(n_r\approx 0.036\)
- \(\SI{200}{mK}\): \(\bar n_{\mathrm{th}}\approx 0.229\), simulated cavity occupation \(n_r\approx 0.229\)
- \(\SI{500}{mK}\): \(\bar n_{\mathrm{th}}\approx 1.044\), simulated cavity occupation \(n_r\approx 1.044\)
- \(\SI{2}{K}\): \(\bar n_{\mathrm{th}}\approx 5.47\), simulated cavity occupation \(n_r\approx 5.27\)

The weakly dressed qubit excited-state population rose with temperature but saturated strongly:

- \(\SI{50}{mK}\): \(P_e\approx 1.2\times 10^{-3}\)
- \(\SI{100}{mK}\): \(P_e\approx 3.3\times 10^{-2}\)
- \(\SI{200}{mK}\): \(P_e\approx 0.153\)
- \(\SI{500}{mK}\): \(P_e\approx 0.290\)
- \(\SI{2}{K}\): \(P_e\approx 0.332\)

The key takeaway is that the cavity occupation is the most faithful thermometer because it stays monotone and nearly linear in the underlying bath occupation, whereas the qubit population becomes a compressed, saturating proxy once the component is hot. In the present parameter set, even the simple temperature-resolution proxy favored cavity occupation over qubit population across the full sweep. The qubit is still useful because it is experimentally accessible and because any nonzero excited-state population immediately certifies that the system is not in the strictly dispersive idealization.

Two practical secondary observables also emerged:

- the thermal spectroscopy broadening proxy exceeded \(\SI{0.1}{MHz}\) near \(\SI{0.2}{K}\);
- the simple thermal readout-noise penalty dropped below \(0.9\) near \(\SI{0.15}{K}\), meaning the hot-mode noise floor had already increased enough to matter for readout fidelity.

### B. Thermal coherence becomes unacceptable around 0.1 K for the chosen readout parameters
The Ramsey-like dispersive simulations showed a clear crossover from benign low-occupation dephasing to rapid coherence loss:

- \(\SI{30}{mK}\): \(\Gamma_\phi \approx 7.6\times 10^3\,\mathrm{s^{-1}}\), \(T_2 \approx \SI{131}{us}\)
- \(\SI{50}{mK}\): \(\Gamma_\phi \approx 8.0\times 10^3\,\mathrm{s^{-1}}\), \(T_2 \approx \SI{100}{us}\)
- \(\SI{70}{mK}\): \(\Gamma_\phi \approx 1.1\times 10^4\,\mathrm{s^{-1}}\), \(T_2 \approx \SI{41}{us}\)
- \(\SI{100}{mK}\): \(\Gamma_\phi \approx 2.1\times 10^4\,\mathrm{s^{-1}}\), \(T_2 \approx \SI{12.5}{us}\)
- \(\SI{150}{mK}\): \(T_2 \approx \SI{4.1}{us}\)
- \(\SI{200}{mK}\): \(T_2 \approx \SI{2.2}{us}\)

The perturbative analytic estimate

\[
\Gamma_{\phi,\mathrm{ana}} \approx \frac{4\chi^2}{\kappa}\bar n_{\mathrm{th}}(\bar n_{\mathrm{th}}+1)
\]

tracked the numerical low-temperature trend within a factor of about three over \(\SI{50}{mK}\) to \(\SI{100}{mK}\), which is acceptable for a first-pass regime classifier but not precise enough to replace the full simulation once \(\bar n_{\mathrm{th}}\) is no longer very small.

Separating mechanisms mattered. In the exact dispersive model the qubit coherence loss is pure dephasing only. In the weakly dressed model the thermally induced upward rate became comparable to or larger than the dephasing contribution by \(\SI{100}{mK}\), and clearly dominated the total coherence budget by \(\SI{150}{mK}\). That means a hot component can look like "just extra dephasing" only in the idealized dispersive limit. Once even modest dressing or leakage is present, actual qubit heating becomes part of the same failure mechanism.

For the chosen device parameters, a conservative boundary for "safe" operation is below about \(\SI{70}{mK}\), while \(\SI{100}{mK}\) is already borderline and \(\SI{150}{mK}\) is clearly dangerous.

### C. A hot auxiliary mode is dangerous as soon as it is appreciably exchange-coupled
The multimode maps were intentionally harsh: one cold local readout mode, one hot auxiliary mode, and a storage-readout exchange \(J\) standing in for cable-mode or package-mode leakage. Under those assumptions the hot auxiliary mode strongly contaminated the local readout once \(J\) became nonzero.

At the representative hot auxiliary temperature of \(\SI{350}{mK}\):

- with \(J=0\), the readout stayed at its cold baseline occupation \(n_r\approx 1.3\times 10^{-3}\);
- with \(J/2\pi=\SI{1.5}{MHz}\), the induced readout occupation jumped to roughly \(0.30\) to \(0.33\) across the detuning range;
- with \(J/2\pi=\SI{12}{MHz}\), the dressed qubit excited-state population reached as high as \(0.167\) at the most dangerous detuning.

Under a stringent safe criterion of

- \(n_r < 0.05\),
- dressed \(P_e < 0.01\),
- dephasing proxy \(< 1/(\SI{20}{us})\),

only the \(J=0\) line remained safe in this parameter set. The design implication is blunt: a hot resonant structure that remains appreciably exchange-coupled to the local readout mode is not a mild correction, it is a direct failure channel.

The linewidth sensitivity scan reinforced that interpretation. At fixed detuning and \(J/2\pi=\SI{6}{MHz}\), broadening the hot auxiliary mode from \(\SI{300}{ns}\) to \(\SI{60}{ns}\) drove the induced readout occupation from \(0.220\) to \(0.428\) and the dephasing proxy from \(4.7\times 10^5\) to \(1.05\times 10^6\,\mathrm{s^{-1}}\). Broader hot modes were worse, not better, because they coupled thermal weight into the readout more efficiently.

### D. The internal quantum response is sub-microsecond, so slow VTS traces must be external
After an instantaneous step in bath occupation, the simulated intrinsic response times were:

- step to \(\SI{0.20}{K}\): cavity \(\tau \approx \SI{0.15}{us}\), qubit \(\tau \approx \SI{0.21}{us}\)
- step to \(\SI{0.50}{K}\): cavity \(\tau \approx \SI{0.15}{us}\), qubit \(\tau \approx \SI{0.089}{us}\)

This is one of the most useful experimental conclusions in the study. If the actual VTS or component-heating experiment shows milliseconds, seconds, or longer transients, those slow time constants are not being generated by the internal cQED Lindblad response. They must come from macroscopic thermalization of the component, its housing, or its interfaces. That means the quantum model should be coupled to a separate thermal transport model rather than stretched to explain slow hardware dynamics by itself.

## Validation
All targeted validation checks passed:

- zero-temperature limit: both cavity occupation and qubit excitation vanished numerically;
- Bose consistency: at \(\SI{0.20}{K}\) the simulated cavity occupation matched the target \(\bar n_{\mathrm{th}}\) to better than \(10^{-3}\) relative error;
- thermometry truncation: the \(\SI{2}{K}\) cavity occupation changed by only about \(2.0\%\) when the cavity cutoff increased from 30 to 36 photons;
- dephasing truncation: the fitted pure-dephasing rate at \(\SI{0.20}{K}\) changed by less than \(10^{-6}\) relative when the dynamic cutoff increased from 12 to 14 photons;
- analytic scaling: the low-temperature dephasing formula stayed within a factor of three over \(\SI{50}{mK}\) to \(\SI{100}{mK}\);
- multimode weak-coupling limit: at \(J=0\) the readout reverted to the cold-bath occupation;
- multimode truncation: representative readout-heating and dressed-excitation spot checks changed by less than \(3\%\) when the multimode cutoff increased from \(6\) to \(8\).

## Experimental interpretation
The simulations support the following measurement priorities.

1. Measure cavity-occupation or spectroscopy-broadening proxies first.
Their behavior is closest to the actual thermal occupation and remains informative even when the qubit thermometer saturates.

2. Use qubit excited-state population as a secondary alarm channel, not the primary calibrated thermometer.
It is excellent for detecting that a leakage channel exists, but it is less linear and more model-dependent than the cavity occupation.

3. Treat \(\SI{100}{mK}\) effective mode temperature as a serious warning threshold for the chosen parameter set.
By that point the simulated thermal-limited coherence has already fallen to about \(\SI{12.5}{us}\).

4. In multimode experiments, prioritize aggressive mode detuning and coupling suppression before trying to infer small thermal effects.
Once the hot auxiliary mode remains appreciably exchange-coupled, the thermal back-action is large enough to overwhelm "precision thermometry" logic.

5. Interpret slow VTS transients as thermal-transport data, not as intrinsic cQED response.
The quantum subsystem equilibrates in sub-microsecond times under the present model, so long transients should be fit with a separate thermal network whose output feeds the bath occupation in the cQED layer.

## Limitations and future work
This study does not model boundary resistance, distributed cable temperature profiles, dielectric relaxation, or non-Markovian thermal noise. The hot component appears only through an effective bath occupation or a hot auxiliary mode. Exact project slide parameters for the SNAILmon-cable system were not available in the workspace, so the multimode numbers should be read as regime-defining rather than device-final.

The cleanest next extension is modular:

1. Build a lumped thermal RC or finite-element model for the component and package.
2. Convert that model into a time-dependent effective bath occupation or mode temperature.
3. Feed the resulting \(n_{\mathrm{th}}(t)\) into the present cQED layer to predict qubit excitation, dephasing, and readout degradation.

That division of labor is realistic because the present study already showed that the quantum subsystem responds much faster than the hardware thermalization.
