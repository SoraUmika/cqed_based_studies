# Review Request

Date: 2026-04-08
Study: `studies/storage_active_cooling_gf_sideband`
Run: `task_runs/storage_active_cooling_gf_sideband`
Status: Ready for review

## What Is Ready
- Technical report written and compiled:
  - `studies/storage_active_cooling_gf_sideband/report/report.tex`
  - `studies/storage_active_cooling_gf_sideband/report/report.pdf`
- Reproducibility notebook created and executed:
  - `studies/storage_active_cooling_gf_sideband/scripts/reproducibility_notebook.ipynb`
- Full artifact set exported under:
  - `studies/storage_active_cooling_gf_sideband/data/`
  - `studies/storage_active_cooling_gf_sideband/artifacts/`
  - `studies/storage_active_cooling_gf_sideband/figures/`

## Main Claims To Review
1. The ladder `|g,0_r,n_s> -> |f,0_r,n_s-1> -> |g,1_r,n_s-1> -> |g,0_r,n_s-1>` is viable in the effective local device model through `n_s <= 4`.
2. The correct Step A mechanism in the present framework is an effective storage red sideband, not a direct transmon `g-f` carrier.
3. The best practical Step A family is paper-motivated `bump` control for `n=1,3,4`, while Step B is best served by a `cosine_squared` pulse for all `n`.
4. The full open-system primitive achieves about `96-97%` single-cycle success and cools basis, coherent, and thermal-like inputs under repeated application.
5. The report keeps the critical modeling boundary explicit: the study validates effective-control feasibility, not microscopic pump generation.

## Requested Reviewer Focus
- Evidence-claim mapping for the main viability claim and the experiment-facing recommendation
- Whether the comparison to `arXiv:2503.10623v1` is specific enough and accurately distinguishes analogy from protocol identity
- Whether the `|q,n_r,n_s>` notation, unit conventions, and internal-vs-reported ordering are clear throughout
- Whether the Floquet/Stark discussion is appropriately cautious given the truncation-boundary warning
- Whether the remaining limitations around calibration transfer and spectral crowding are stated strongly enough
