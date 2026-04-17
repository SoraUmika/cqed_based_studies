# Auto Workflow Probe Study

## Problem Class
ANA

## Motivation
This probe study exercises the staged research loop on the smallest physically meaningful cQED question in the repository. The target question is whether free dispersive idle evolution alone can accumulate a conditional qubit phase difference of $\pi$ between the cavity $n=0$ and $n=1$ manifolds over a finite time. The scope stays deliberately minimal so the loop can be validated without heavy numerics: derive the condition analytically, confirm it with two tiny cqed_sim-backed calculations, and save a compact figure-artifact chain. This study is intentionally framed as a repository-internal model-validation note rather than as a calibrated hardware-realism claim.

## Goals
1. Derive the first-principles analytic condition for a $\pi$ relative qubit phase between the $n=0$ and $n=1$ cavity manifolds under free dispersive evolution.
2. Confirm the analytic result numerically with the lightest viable cqed_sim workflow on a tiny Hilbert space.
3. Add one independent framework-backed check based on static-Hamiltonian evolution rather than only the closed-form phase helper.
4. Keep the implementation lightweight enough to finish comfortably within the requested two-minute CPU budget.
5. Produce the minimum research-loop outputs needed for downstream validation and report revision.

## Methods
- Use the local cqed_sim API reference to anchor the study design before implementation.
- Model the system with `cqed_sim.core.model.DispersiveTransmonCavityModel` using a minimal qubit-storage dispersive Hamiltonian.
- Use `cqed_sim.gates.coupled.dispersive_phase` as the primary lightweight numerical confirmation path for the idle unitary.
- Use the model's static Hamiltonian in the rotating frame defined by the bare qubit and cavity frequencies for an independent idle-evolution cross-check on the same tiny Hilbert space.
- Use basis-energy and manifold-frequency helpers only as sign-convention and branch-shift audits.

## Analytic Preliminary
Start from the first-order dispersive Hamiltonian in the rotating frame,
$$
H/\hbar = \chi \hat{n} |e\rangle\langle e|,
$$
where $\hat{n}$ is the cavity photon-number operator and $\chi$ is the storage-qubit dispersive shift. In this convention, the excited-state branch in cavity manifold $n$ accumulates phase
$$
\phi_n(t) = -\chi n t.
$$
The relative phase between the $n=1$ and $n=0$ manifolds is therefore
$$
\Delta \phi_{1,0}(t) = \phi_1(t) - \phi_0(t) = -\chi t \pmod{2\pi}.
$$
Hence a conditional phase difference of $\pi$ occurs at the finite idle time
$$
t_{\pi} = \frac{\pi}{|\chi|}, \qquad \chi \neq 0.
$$
This is the minimal analytic answer: free dispersive evolution alone is sufficient whenever the dispersive shift is nonzero.

This probe therefore tests a model identity within the first-order dispersive Hamiltonian itself. It does not, by itself, establish that any particular calibrated device satisfies the dispersive-regime validity condition because no explicit coupling-detuning pair is instantiated in this study.

Controlled approximations used for this probe:
- First-order dispersive model only.
- Closed-system idle evolution with no drive and no decoherence.
- Low-occupation focus on $n=0,1$, for which the cavity self-Kerr contribution $K n(n-1)/2$ vanishes identically.
- Validity of the dispersive picture requires the usual regime assumptions, e.g. $|g/\Delta| \ll 1$ and low enough occupation that higher-order corrections remain perturbative.

## cqed_sim Gap Analysis
| Functionality | Needed? | Available in cqed_sim? | Plan |
|---|---|---|---|
| Minimal qubit-cavity dispersive model | Yes | Yes | Use `DispersiveTransmonCavityModel` for the implementation script. |
| Closed-form idle dispersive unitary | Yes | Yes | Use `dispersive_phase` for the fastest numerical confirmation. |
| Full pulse-sequence simulation of idle evolution | Optional | Yes | Skip for the minimal probe unless the lightweight helper check is insufficient. |
| Artifact and figure generation | Yes | Partially | Use local study script plus matplotlib for saved outputs. |

## Assumptions
- Representative parameters will use a storage-qubit dispersive shift near the repository default, e.g. $\chi = -2\pi \times 2.84\,\mathrm{MHz}$ expressed in rad/s.
- Internal calculations will use consistent rad/s and seconds units.
- A cavity truncation of at least two Fock states is sufficient for the target $n=0$ versus $n=1$ phase-difference question.
- Numerical agreement will be judged after phase wrapping, because the relevant observable is the relative phase modulo $2\pi$.
- The study makes no separate claim about a specific hardware value of $g/\Delta$; it is restricted to the implemented first-order dispersive model.

## Compute & Resource Strategy
The study is intentionally tiny. The main deliverable is analytic, and the numerical confirmation only requires evaluating a small dispersive unitary plus a static-Hamiltonian exponential on a 2x2 or 2x3 cavity space. No GPU, multiprocessing, or long-running sweeps are needed. Expected runtime for the revision implementation phase remains under 10 seconds end-to-end, with the full probe remaining far below the two-minute cap.

## Expected Outcomes
- Analytic result: $t_{\pi} = \pi/|\chi|$ for nonzero $\chi$.
- Representative numerical result: for $\chi = -2\pi \times 2.84\,\mathrm{MHz}$, the first $\pi$ crossing should occur near $176\,\mathrm{ns}$.
- One machine-readable artifact containing the chosen parameters, analytic $t_{\pi}$, and sampled numerical phase differences.
- One simple line plot of the relative phase versus idle time with the $\pi$ crossing marked.
- One additional figure-artifact pair showing that explicit static-Hamiltonian evolution agrees with the closed-form phase law over the sampled idle-time scan.

## Known Limitations
- This probe addresses only the ideal closed-system first-order dispersive model.
- It does not include higher-order dispersive terms, decoherence, readout-mode effects, or pulse-shape imperfections.
- It is a workflow-validation study rather than a comprehensive physics benchmark.
- It does not quantify a hardware dispersive-regime bound such as $g/\Delta$ and therefore should not be read as a calibrated device statement.

## Validation
- [x] Sanity checks - PASSED. The saved probe data confirm $\Delta\phi(0)=0$, the first $\pi$ crossing at $t_{\pi} = 176.056338\,\mathrm{ns}$, and machine-precision agreement between the analytic phase law, the closed-form framework helper, and the independent static-Hamiltonian evolution path, with a maximum wrapped mismatch of $8.882 \times 10^{-16}\,\mathrm{rad}$.
- [x] Convergence - PASSED. Increasing the cavity truncation from $n_\mathrm{cav}=2$ to $n_\mathrm{cav}=3$ changed the wrapped phase trace by $0.0\,\mathrm{rad}$ over all 51 samples, and the helper and Hamiltonian traces at the three-level cavity cutoff agree exactly after phase wrapping; no time-step sweep is required because both numerical paths use direct unitary evaluation rather than ODE integration.
- [x] Literature comparison (if applicable) - Not applicable for this internal workflow probe, which validates the first-principles dispersive result against a minimal cqed_sim confirmation rather than against an external benchmark.

## Status
ACTIVE
