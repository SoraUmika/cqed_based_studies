# Literature Design Memo: Hybrid Universal Control in a 2x2 Qubit+Cavity Logical Subspace

## Scope
This memo positions the main control paradigms relevant to the present study:

1. `SNAP + displacement + qubit rotations`
2. `SQR`-like selective hybrid control
3. `ECD`-like conditional-displacement control
4. `GRAPE` / direct waveform-level optimal control
5. Fast conditional-displacement variants such as `CNOD`
6. Native activated nonlinear bosonic interactions as an upper-benchmark hardware direction

The benchmark question is deliberately practical: not "what is universal in principle?", but "what is the best realizable route on the local device model already encoded in `cqed_sim` for a first 2x2 qubit+cavity logical block?"

## Core References
- Reinier W. Heeres et al., ["Cavity State Manipulation Using Photon-Number Selective Phase Gates"](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.137002), *Phys. Rev. Lett.* 115, 137002 (published September 22, 2015).
- Stefan Krastanov et al., ["Universal Control of an Oscillator with Dispersive Coupling to a Qubit"](https://arxiv.org/abs/1502.08015), *Phys. Rev. A* 92, 040303(R) (2015).
- Reinier W. Heeres et al., ["Implementing a universal gate set on a logical qubit encoded in an oscillator"](https://www.nature.com/articles/s41467-017-00045-1), *Nature Communications* 8, 94 (published July 21, 2017).
- Alec Eickbusch et al., ["Fast universal control of an oscillator with weak dispersive coupling to a qubit"](https://www.nature.com/articles/s41567-022-01776-9), *Nature Physics* 18, 1464-1469 (2022).
- A. Wallraff et al., ["Sideband Transitions and Two-Tone Spectroscopy of a Superconducting Qubit Strongly Coupled to an On-Chip Cavity"](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.050501), *Phys. Rev. Lett.* 99, 050501 (published July 30, 2007).
- Asaf A. Diringer et al., ["Conditional-not Displacement: Fast Multioscillator Control with a Single Qubit"](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.14.011055), *Phys. Rev. X* 14, 011055 (2024).
- Ni-Ni Huang et al., ["On-demand transposition across light-matter interaction regimes in bosonic cQED"](https://www.nature.com/articles/s41467-024-50201-7), *Nature Communications* 15 (2024).
- Jordan Huang et al., ["Fast Sideband Control of a Weakly Coupled Multimode Bosonic Memory"](https://arxiv.org/abs/2503.10623), arXiv:2503.10623 (2025).
- Navin Khaneja et al., ["Optimal control of coupled spin dynamics: design of NMR pulse sequences by gradient ascent algorithms"](https://doi.org/10.1016/j.jmr.2004.11.004), *J. Magn. Reson.* 172, 296-305 (2005).
- `cqed_sim` software reference: [SoraUmika/qubox_cQEDsim](https://github.com/SoraUmika/qubox_cQEDsim).

## Short Positioning Table
| Family | Universality status | Natural support for exact hybrid qubit+cavity control | Expected speed on dispersive hardware | Calibration burden | First-pass expectation |
| --- | --- | --- | --- | --- | --- |
| SNAP + D + qubit rotations | Strong theoretical and experimental basis for oscillator universality with ancilla | Good once combined with ancilla rotations and native dispersive phase | Usually slowest because of number selectivity | Moderate; primitives are interpretable and familiar | Strong baseline, especially for local cavity control |
| SQR-like selective hybrid control | Naturally hybrid, but still selective | Good entangling structure, especially number-conditioned qubit action | Slow-to-moderate because it is still selective | Moderate-to-high | Potential depth reduction, but duration may dominate |
| ECD-like conditional displacement | Strong oscillator-control route, especially in weak dispersive regime | Hybrid control is plausible when paired with fast qubit rotations and dispersive phase | Fast | Moderate if compiled to echoed pulses; low if treated ideally only | Serious candidate, but encoding-dependent |
| CNOD-style fast conditional displacement | Very expressive fast control family | Good for multimode/hybrid conditional action | Fast | Likely high until native tooling exists | High-upside future direction |
| GRAPE | No gate-library universality statement needed; direct pulse optimization | Excellent in principle | Potentially very fast | Highest | Best as a lower-bound benchmark, not automatically the best experimental choice |
| Native activated nonlinear interactions | Hardware-native expressivity can be extremely high | Potentially excellent | Fastest conceptual upper benchmark | Depends on hardware availability, often high integration cost | Useful benchmark, not current-device baseline |
| SWAP-/sideband-style exchange | Naturally hybrid and strongly tied to native JC/sideband dynamics | Good for state transfer and entangling primitives, less direct for exact logical phase structure | Fast | Moderate-to-high | Serious native candidate; likely strong for transfer, less automatic for exact logical CX |

## Family-by-Family Notes

## 1. SNAP + displacement + qubit rotations
The Heeres 2015 SNAP paper established the selective number-dependent arbitrary phase gate as an experimentally practical dispersive primitive and explicitly noted that combining SNAP with oscillator displacements yields arbitrary unitary control of oscillator-encoded qubits. Krastanov et al. then gave the general universality proof and constructive recipes for state preparation and arbitrary oscillator unitaries in the strong-dispersive regime.

For the present study, the main practical implications are:

- The route is theoretically mature and experimentally grounded.
- It is especially strong for local cavity control because `D` mixes neighboring Fock levels while `SNAP` repairs the number-dependent phases.
- It is not automatically the fastest route for hybrid entanglers because selective operations cost time.
- It ports cleanly to experiment because each primitive is interpretable and calibratable.

Conclusion for the benchmark: this is the right baseline library, and any winning alternative should beat it on either speed, leakage, or calibration simplicity without sacrificing too much exactness.

## 2. SQR-like selective hybrid control
The selective qubit-rotation idea is not the original bosonic-universality construction, but it is highly relevant for hybrid control because it acts directly on the ancilla conditioned on cavity number. That is exactly the structural feature needed for logical qubit+cavity entanglers.

Practical expectations:

- Better natural access to hybrid entanglers than pure `SNAP + D`.
- Still selective, so likely slower than echoed or activated conditional-displacement routes.
- May produce the "right entangling block" but with additional number-dependent logical phases that must be corrected separately.

Conclusion for the benchmark: plausible entangler candidate, but probably not the best all-around universal route unless paired with a cheap cavity-phase correction.

## 3. ECD-style conditional displacement
Eickbusch et al. showed that echoed conditional displacement can provide fast universal oscillator control even when the bare dispersive coupling is weak, precisely because the method avoids relying on narrow selective spectral resolution for every control action. This makes ECD especially attractive whenever `chi` is not large enough to make selective control cheap.

Practical expectations:

- Faster than SNAP-like selective constructions.
- Very promising when the logical encoding itself likes displacement-based action, especially cat-like or displaced-state encodings.
- More subtle in a strict low-Fock logical block, where a large conditional displacement may create apparent leakage unless it is paired with exactly the right phase-space geometry.

Conclusion for the benchmark: ECD is a serious candidate in principle, but the first-pass Fock `{|0>,|1>}` encoding is not obviously its most natural operating point.

## 4. GRAPE / direct optimal control
GRAPE remains the cleanest pulse-level benchmark for "what might be achievable in principle" under the modeled Hamiltonian and control channels. The 2017 oscillator-logical-gate paper by Heeres et al. is especially relevant here because it demonstrates that a holistic Hamiltonian-aware optimal-control strategy can implement a universal logical gate set on an oscillator-encoded qubit. In other words, direct optimal control is not just a theoretical lower bound; it has already informed experimentally useful bosonic control.

Practical expectations:

- Excellent benchmark for short logical entanglers.
- Can reveal new effective primitives or compressions.
- Usually carries the heaviest calibration and transfer burden, especially when the pulse depends delicately on model accuracy.

Conclusion for the benchmark: mandatory as a reference, but not automatically the preferred final implementation route.

## 5. Fast conditional-displacement variants / CNOD
The 2024 CNOD paper by Diringer et al. pushes the conditional-displacement idea toward fast multioscillator control with a single ancilla. Even if the exact primitive is not yet implemented locally, the conceptual lesson is important: unselective but structured hybrid primitives can outperform slow selective decompositions by moving the complexity into the geometric form of the drive instead of into spectral selectivity.

Practical expectations:

- High potential for future hybrid or multimode scaling.
- Attractive if the study later expands beyond a single cavity logical qubit.
- Requires first-class framework support before it can be compared fairly against the more mature libraries.

Conclusion for the benchmark: not yet a fair apples-to-apples current-device baseline, but a strong future-work direction and an argument against overcommitting to selective control for larger logical spaces.

## 6. SWAP-like / sideband / exchange-style native control
SWAP-like native control is the most direct alternative to selective dispersive gates when the device can transiently access red-sideband, blue-sideband, Jaynes-Cummings, or other exchange-dominated dynamics. The foundational circuit-QED sideband work by Wallraff et al. already framed sidebands as a route to controllable qubit-photon and qubit-qubit entanglement. More recent work pushes this much further: Huang et al. demonstrate regime switching between dispersive and resonant light-matter dynamics on demand, while the 2025 fast-sideband-control preprint makes the specific case that transmon-cavity SWAP gates can be accelerated well beyond the bare dispersive scale in weakly coupled multimode bosonic memories.

Practical expectations:

- Native exchange is naturally hybrid, not merely oscillator-only.
- The most natural unitary is often an iSWAP- or SWAP-like transfer, not a logical controlled-phase or controlled-X in the computational basis.
- Gate speed can be excellent because the interaction is activated directly rather than built from spectral selectivity.
- Leakage management is central: exchange primitives tend to connect neighboring excitation manifolds, so the physical interaction is attractive exactly where exact logical isolation becomes harder.

Conclusion for the benchmark: SWAP-/sideband-style primitives are serious native candidates and deserve explicit benchmarking, but they should be expected to excel first at fast transfer and entangling structure rather than at exact low-dimensional logical gauge matching.

## 7. Native activated nonlinear interactions
Activated nonlinear hardware routes matter mainly as an upper benchmark for what "best practical" could mean if the device exposes stronger native bosonic interactions. Even if that is not the present platform, these results show that some of the control burden carried by `SNAP`, `SQR`, or GRAPE in dispersive architectures can disappear when the hardware itself supplies a more direct nonlinearity.

Practical expectations:

- Fastest conceptual route.
- Lowest circuit depth if the interaction is already native.
- Usually not directly portable unless the hardware stack already supports the activated interaction.

Conclusion for the benchmark: useful for perspective, but not the primary recommendation for the current study because the local simulator/device assumptions are still dispersive-first.

## Cross-Cutting Answers to the Required Questions

## Which approaches are universal only for the cavity mode, and which naturally support hybrid control?
- `SNAP + D` is fundamentally an oscillator-universality construction with an ancilla. Hybrid qubit+cavity control is obtained by adding qubit rotations and exploiting the same dispersive coupling.
- `ECD` is also usually presented as a route to fast oscillator control with an ancilla, but it naturally creates hybrid correlations and therefore can support hybrid unitary synthesis when paired with ancilla-local operations.
- `SQR`-like primitives are hybrid by construction because they act on the qubit conditioned on cavity number.
- `GRAPE` side-steps this distinction and directly targets the full hybrid Hilbert space.

## Which primitives are selective and therefore slow?
- `SNAP`
- `SQR`
- Any number-resolved qubit-conditioned operation that depends on spectral selectivity over `chi`

## Which are unselective / echoed / activated and therefore potentially faster?
- `ECD`
- `CNOD`
- Native activated nonlinear interactions
- Some GRAPE pulses, when they discover a strong-drive geometric solution rather than a selective one

## Which routes are easiest to port from simulation to experiment?
- First: calibrated native dispersive waiting plus fast local rotations/displacements
- Second: `SNAP + D + qubit rotations`
- Third: `SQR` if it is already calibrated in the lab stack
- Last: raw GRAPE, unless a pulse family can be compressed into an interpretable primitive

## Which routes are likely most robust to weak-dispersive constraints, finite ancilla coherence, and cavity truncation?
- `ECD` and `CNOD` are attractive under weak-dispersive constraints because they reduce reliance on spectral selectivity.
- SWAP-/sideband-style native exchange is also attractive in weak-dispersive hardware because it converts the problem from high spectral selectivity to activated hybrid transfer, but it introduces a stronger leakage-management burden.
- `SNAP` and `SQR` become costly when selectivity is expensive.
- Raw GRAPE can sometimes outperform hand-built decompositions, but robustness is not guaranteed unless it is optimized explicitly for uncertainty.
- Low-Fock truncation can make displacement-based schemes look artificially leaky if the encoding is not chosen to match the primitive.

## Native vs Selective: What the literature suggests before numerics
- Selective routes (`SNAP`, `SQR`) are favored when exact number resolution and logical phase control are more important than speed.
- Native routes (`ECD`, JC/SWAP-like exchange, blue-sideband, activated cubic interactions) are favored when hardware realism rewards stronger interactions and shorter exposure to decoherence.
- Native exchange is especially compelling when the task is state transfer, encoding, or iSWAP-like entanglement.
- Selective control remains easier to reason about and calibrate for exact logical basis-change operations in a low-photon Fock code.
- The likely practical winner for a 2x2 Fock benchmark is therefore a hybrid answer: exact local control from selective dispersive primitives plus fast entanglers from native interactions.

## Design Implications for the Present Benchmark
The literature suggests a clear expectation before numerics:

1. A selective dispersive route should remain the safest baseline for exact control in the tiny Fock-encoded 2x2 logical block.
2. A native dispersive entangler may outperform more elaborate libraries for the specific entangling target.
3. SWAP-/sideband-style exchange should be benchmarked explicitly because it is the most natural native hybrid primitive after dispersive conditional phase, but it may suffer from leakage or gauge-mismatch when judged against an exact logical `CX`.
4. ECD-like control may be disadvantaged by the chosen encoding even if it becomes favorable later for cat-like logical bases or larger logical dimensions.
5. GRAPE should be treated as a benchmark and a source of intuition, not as the assumed winner.

That is exactly the comparison implemented in this study.
