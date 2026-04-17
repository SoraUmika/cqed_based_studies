# Study Prompt: Multitone SQR No-Go, Decoupled-Block Limit, and Echo Alternative in Dispersive cQED

Study Name: multitone_sqr_no_go_decoupled_echo_dispersive_cqed

## Context

You are working in an environment where **many closely related studies have already been attempted**, including analyses of SQR, CPSQR, multitone control, tomography-based validation, and various simulation studies. **Do not assume those prior studies are correct.** Treat them as potentially useful prior art, but **verify all important claims from first principles, simulation, code inspection, and consistency checks**.

The goal of this study is to rigorously examine the claim that, in a dispersive cQED system, an **ideal SQR gate cannot be achieved exactly by a simultaneous multitone ansatz** if we **do not allow artificial per-tone detunings** to cancel induced \(Z\)-precession, and if the only control knobs are:
- per-tone amplitude correction
- per-tone azimuthal-axis correction

You may assume:
- dispersive regime
- no leakage between neighboring cavity Fock states for the baseline analytical argument
- block-diagonal structure in Fock number under the dispersive Hamiltonian
- the multitone drive is applied on a **shared qubit control line**

Then, after establishing the no-go result carefully, you must study two alternatives:

1. A **stronger decoupled-block approximation** in which, for each addressed Fock block, one keeps only its resonant tone and drops all off-resonant tones exactly. Show that under this stronger approximation the no-go statement is false and ideal blockwise SQR becomes achievable. Also explain that this is effectively no longer the same simultaneous-shared-line multitone problem; operationally it corresponds to addressing blocks one at a time or otherwise assuming perfect decoupling.

2. A proposed **multitone echo pulse sequence** for SQR:
   \
   \[
   \text{half-SQR} \rightarrow \pi \rightarrow \text{half-SQR} \rightarrow \pi
   \]
   Study whether and when this can mitigate the unwanted \(Z\)-precession while preserving the desired conditional \(XY\) rotation.

This study must be done skeptically and carefully. If previous studies in the environment appear incomplete, inconsistent, or under-validated, say so explicitly and correct them.

---

## Primary Objectives

You must complete **all** of the following objectives.

### Objective A — Formal analytical no-go result for simultaneous multitone SQR

Consider the dispersive Hamiltonian in a Fock-resolved form:
\[
H_0 = \sum_n \frac{\omega_n}{2}\,\sigma_z \otimes |n\rangle\langle n|
\]
with
\[
\omega_n = \omega_q + \chi n + \chi' n(n-1) + \cdots
\]

Use a multitone shared-line qubit drive:
\[
H_d(t) = \sum_{m \in \mathcal S} \frac{\Omega_m}{2}
\left(
 e^{-i(\omega_m t - \phi_m)} \sigma_+ +
 e^{+i(\omega_m t - \phi_m)} \sigma_-
\right)
\]
with **no artificial detuning**, meaning each tone is placed exactly at the chosen block frequency and there is no independent per-block \(Z\)-compensation term.

The target ideal SQR should be defined as:
\[
U_{\mathrm{SQR}}^{\mathrm{ideal}} = \bigoplus_n R_{\hat n(\varphi_n)}(\theta_n)
\]
where each block is a pure \(XY\)-plane qubit rotation with no residual \(Z\)-phase.

You must:
1. derive the interaction-frame Hamiltonian block-by-block,
2. keep the exact influence of off-resonant tones rather than discarding them prematurely,
3. show that off-resonant tones generate effective \(Z\)-precession terms,
4. explain why amplitude and azimuth knobs are insufficient to cancel those terms generically,
5. present a clean two-block proof first,
6. then generalize to many addressed Fock blocks,
7. state precisely what “generically impossible” means,
8. clearly identify all assumptions used.

The proof should explicitly use either:
- Magnus expansion / average Hamiltonian theory, or
- an equivalent controlled perturbative argument

and should make it completely clear **where** the obstruction comes from.

---

### Objective B — Numerical verification and attempted falsification of the no-go claim

The study must not stop at the analytical argument. You must try to **numerically falsify** the no-go claim under realistic simulation assumptions.

Use the simulation stack available in the repository/environment and do not rely on handwaving. Build a numerical experiment that compares:

1. **Exact simultaneous multitone evolution** under the shared-line dispersive Hamiltonian
2. The target **ideal SQR**
3. The effective blockwise description inferred from the analytical argument

At minimum:
- truncate to a manageable Fock space
- choose parameters representative of the user’s working regime if available
- otherwise use clearly stated representative dispersive parameters
- simulate for multiple gate durations and rotation targets
- vary amplitude/phase knobs extensively
- explicitly test whether optimization over amplitudes and azimuths alone can recover the ideal SQR unitary or ideal logical action

You must report:
- best achieved fidelity to the ideal SQR
- whether any apparent “success” disappears when examined blockwise
- induced \(Z\)-phases per Fock block
- dependence on addressed subspace size
- sensitivity to gate duration

The numerical section should be framed as:
- **attempted falsification**
- if the no-go still holds numerically, explain why
- if an apparent counterexample is found, investigate whether it is due to approximation, metric choice, underconstrained optimization, leakage tolerance, or hidden effective detuning

Do not merely optimize a single scalar fidelity and declare victory. Inspect the blockwise action carefully.

---

### Objective C — Prove the stronger decoupled-block approximation *does* allow ideal SQR

Now impose the stronger approximation:

> In each Fock block, keep only the resonant tone for that block and drop all off-resonant tones exactly.

Under this approximation, show analytically that the no-go result no longer applies.

You must:
1. write the reduced blockwise Hamiltonian,
2. show that each block becomes an independently driven \(XY\)-rotation,
3. prove that ideal SQR is then realizable by choosing amplitude and azimuth appropriately,
4. explain clearly why this model is stronger than the physical simultaneous multitone shared-line model,
5. connect this approximation to the physical interpretation that one is effectively driving Fock levels in sequence or otherwise assuming exact decoupling.

Then verify this numerically:
- simulate the decoupled-block model directly,
- confirm that ideal SQR can be matched,
- compare against the full simultaneous multitone model,
- make the distinction crystal clear.

This section is important because it explains why some earlier studies may have “worked” if they implicitly or explicitly assumed blockwise decoupling.

---

### Objective D — Study the proposed multitone echo SQR sequence

Study the proposed alternative pulse sequence:

\[
\text{half-SQR} \rightarrow \pi \rightarrow \text{half-SQR} \rightarrow \pi
\]

Interpret this carefully in the dispersive cQED context. The purpose is to investigate whether echoing can cancel or strongly suppress the unwanted \(Z\)-precession induced during multitone SQR while preserving the target conditional \(XY\) rotation.

You must:
1. define the exact sequence convention used,
2. specify whether the \(\pi\) pulses are broadband, selective, or idealized instantaneous,
3. derive the toggling-frame / echo-frame picture,
4. analyze which \(Z\)-terms are canceled and which are not,
5. determine whether the desired SQR action survives with the correct net angle and axis,
6. compare the echoed sequence against plain simultaneous multitone SQR,
7. quantify the cost in gate duration and complexity,
8. identify failure modes and hidden assumptions.

If the echo sequence only works approximately, say so explicitly. If it works under a narrower set of assumptions than expected, spell that out.

---

## Required Deliverables

You must produce the following outputs.

### 1. A rigorous written report
The report should include:
- statement of the problem
- assumptions
- analytical derivations
- numerical methodology
- results
- interpretation
- limitations
- final conclusion

The writing should be explicit and skeptical. If prior work in the environment made stronger claims than justified, correct the record.

### 2. Reproducible code and/or notebooks
Provide simulation code that can be rerun. The code should:
- define the full simultaneous multitone Hamiltonian
- define the ideal target SQR
- define the stronger decoupled-block approximation model
- define the echoed multitone sequence
- evaluate suitable fidelity / blockwise metrics
- extract effective blockwise \(X_n, Y_n, Z_n\) behavior where appropriate

### 3. Clear plots
At minimum include plots of:
- achieved fidelity vs gate duration
- blockwise induced \(Z\)-phase vs duration
- comparison of full model vs decoupled approximation
- comparison of plain multitone vs echo multitone
- any optimization landscape or failure case plots that help make the conclusion convincing

### 4. Final verdict section
End with a concise section answering these questions directly:

1. Is exact ideal simultaneous multitone SQR possible under the physical shared-line dispersive model with only amplitude and azimuth knobs and no artificial detuning?
2. Under what assumptions is it impossible?
3. Under what stronger approximation does ideal SQR become achievable?
4. Does the echoed multitone sequence rescue the situation exactly, approximately, or only in special regimes?

---

## Strong Guidance on Rigor

### Do not trust earlier related studies by default
The environment likely contains previous attempts. You may use them for context, but:
- inspect their assumptions,
- check whether they silently dropped off-resonant terms,
- check whether their metric ignored blockwise \(Z\)-phase,
- check whether “success” only meant a loose state-transfer metric instead of an actual gate/unitary match,
- check whether they accidentally introduced effective detuning or hidden compensation.

If something is invalid or incomplete, state that plainly.

### Distinguish exact from approximate claims
You must be extremely careful to separate:
- exact realizability
- approximate realizability
- realizability under RWA
- realizability under a stronger decoupled-block approximation
- realizability only on a restricted input family
- realizability only up to a tolerated block-dependent \(Z\)-phase

Do not blur these together.

### Use multiple metrics
A single fidelity number is not sufficient. Use multiple diagnostics, such as:
- logical-subspace unitary fidelity
- blockwise axis-angle reconstruction
- blockwise effective \(Z\)-phase accumulation
- target-state fidelity on a set of probe states
- robustness to different initial cavity superpositions

If a metric can hide the failure mode, say so.

---

## Suggested Technical Structure

A good structure would be:

1. **Model setup**
   - dispersive Hamiltonian
   - multitone drive ansatz
   - ideal SQR definition

2. **Formal no-go derivation**
   - interaction-frame block decomposition
   - off-resonant tone terms
   - second-order / Magnus \(Z\)-term
   - two-block impossibility proof
   - many-block generalization

3. **Numerical falsification attempt**
   - simulation details
   - optimization over amplitudes and azimuths
   - results and failure analysis

4. **Stronger decoupled-block approximation**
   - analytical proof of realizability
   - direct numerical confirmation
   - why this is not the same physical problem

5. **Echoed multitone SQR**
   - sequence definition
   - toggling-frame analysis
   - numerical comparison
   - when it helps and when it does not

6. **Synthesis and final verdict**

---

## Parameter Guidance

If there are known user-relevant parameters in the environment, prefer them. Otherwise choose representative values and state them clearly. You may include dispersive nonlinearities such as:
- \(\chi\)
- \(\chi'\)
- self-Kerr if relevant

But the baseline no-go argument should already emerge in the simpler dispersive block-resolved setting.

If you include more realistic corrections later:
- clearly mark them as extensions,
- do not let them obscure the core result.

---

## What to Watch Out For

Common failure modes in this kind of study include:
- silently dropping off-resonant tones too early
- using a fidelity that is insensitive to block-dependent \(Z\)-phase
- validating only a narrow set of initial states
- conflating “good approximation” with “exact realizability”
- claiming the echo works without checking the preserved target rotation carefully
- using an optimization that only finds a locally good-looking but fundamentally non-ideal solution

Avoid all of these.

---

## Minimum Standard for Acceptance

Your study is not complete unless it does all of the following:
- presents a serious analytical argument,
- attempts to numerically disprove that argument,
- explains why the stronger decoupled-block model *does* permit ideal SQR,
- studies the proposed multitone echo sequence carefully,
- and produces a final answer that is cautious, explicit, and technically defensible.

If results remain inconclusive in any part, say so explicitly and explain exactly what remains unresolved.

---

## Final Instruction

Be skeptical. Be explicit. Do not overclaim. If prior work in this environment was suggestive but not actually valid, correct it. The purpose of this study is to determine what is **actually true** about simultaneous multitone SQR, the stronger decoupled-block approximation, and the multitone echo alternative in dispersive cQED.
