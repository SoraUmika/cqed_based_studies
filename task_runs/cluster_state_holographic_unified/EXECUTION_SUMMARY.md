# Execution Summary

- Implemented the corrected-scope rerun in studies/cluster_state_holographic_unified/scripts/run_corrected_scope_study.py and added a pure D+R+CPSQR builder in common.py.
- Restricted the structured-family comparison to D + R + SQR and D + R + CPSQR only, with physical re-evaluation at N_cav = 10, 12, 14 and default final evidence at N_cav = 12.
- Applied the explicit active-subspace rule with threshold 1e-3 and 99.9% captured population reporting for each shortlisted candidate.
- Corrected outcome:
	- Best preliminary D + R + SQR: F(N_cav=4) = 0.9993513138326724, but discarded physically.
	- Best physical D + R + SQR: F(12) = 0.9270750608631796, leakage_worst(12) = 0.16854221060407115, active support touches the truncation boundary.
	- Best retained D + R + CPSQR: F(12) = 0.9665477058496347, leakage_worst(12) = 0.06408327089881649, active support 0..7, retained low confidence.
	- GRAPE reference unchanged: replay fidelity 0.9561, open-system process fidelity 0.9009 at 400 ns.
- Regenerated corrected local figures, rewrote README.md, IMPROVEMENTS.md, report.tex, and scripts/reproducibility_notebook.ipynb, and compiled report/report.pdf successfully.