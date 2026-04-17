"""Test model-based synthesis for the cluster target."""
import sys, numpy as np, time
sys.path.insert(0, r'C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation')

from cqed_sim import DispersiveTransmonCavityModel
from cqed_sim.unitary_synthesis import (
    Displacement as SDisp, QubitRotation as SRot, SNAP as SSNAP,
    Subspace, TargetUnitary, UnitarySynthesizer)
from cqed_sim.unitary_synthesis.targets import make_target

TWO_PI = 2*np.pi
NCAV=8; NTR=2
OQ=TWO_PI*6.150e9; OC=TWO_PI*5.241e9
AL=TWO_PI*(-255e6); CH=TWO_PI*(-2.84e6)
CP=TWO_PI*(-21e3);  KR=TWO_PI*(-28e3)

model = DispersiveTransmonCavityModel(
    omega_c=OC, omega_q=OQ, alpha=AL,
    chi=CH, chi_higher=(CP,), kerr=KR, n_cav=NCAV, n_tr=NTR)

U = make_target('cluster', n_match=1)
sub = Subspace.custom(NTR*NCAV, (0,1,NCAV,NCAV+1),
                      ('|g0>','|g1>','|e0>','|e1>'))
tgt = TargetUnitary(U, ignore_global_phase=True)

Ds=48e-9; Rs=16e-9; Ss=200e-9

# Test with model: D-R-SNAP-D-R-SNAP-D-R
print("=== Model-based: D-R-SNAP-D-R-SNAP-D-R (2 SNAP) ===")
gl = [
    SDisp(name='D1', alpha=.2+0j, duration=Ds, optimize_time=False),
    SRot(name='R1', theta=.1, phi=0., duration=Rs, optimize_time=False),
    SSNAP(name='SNAP1', phases=[0.]*NCAV, duration=Ss, optimize_time=False),
    SDisp(name='D2', alpha=.2+0j, duration=Ds, optimize_time=False),
    SRot(name='R2', theta=.1, phi=0., duration=Rs, optimize_time=False),
    SSNAP(name='SNAP2', phases=[0.]*NCAV, duration=Ss, optimize_time=False),
    SDisp(name='D3', alpha=.2+0j, duration=Ds, optimize_time=False),
    SRot(name='R3', theta=.1, phi=0., duration=Rs, optimize_time=False),
]

t0 = time.time()
try:
    syn = UnitarySynthesizer(
        primitives=gl, subspace=sub, target=tgt,
        model=model, seed=42, optimize_times=False, optimizer='powell')
    r = syn.fit(multistart=3, maxiter=400)
    print(f"  obj={r.objective:.6f} success={r.success} time={time.time()-t0:.1f}s")
    if hasattr(r, 'fidelity'):
        print(f"  fidelity={r.fidelity:.6f}")
    print(f"  report fidelity: {r.report.get('metrics', {}).get('fidelity', 'N/A')}")
    # Print fidelity from metrics
    m = r.report.get('metrics', {})
    for k in ['fidelity', 'nominal_fidelity', 'leakage_average']:
        if k in m:
            print(f"  metrics.{k} = {m[k]:.6f}")
except Exception as e:
    print(f"  FAILED: {e}")
    import traceback; traceback.print_exc()
    print(f"  time={time.time()-t0:.1f}s")

# Also try without model but with 3 SNAPs
print("\n=== No model: D-R-SNAP (3 layers) ===")
gl3 = []
for k in range(3):
    gl3.append(SDisp(name=f'D{k+1}', alpha=.2+0j, duration=Ds, optimize_time=False))
    gl3.append(SRot(name=f'R{k+1}', theta=.1, phi=0., duration=Rs, optimize_time=False))
    gl3.append(SSNAP(name=f'SNAP{k+1}', phases=[0.]*NCAV, duration=Ss, optimize_time=False))
gl3.append(SDisp(name='D4', alpha=.2+0j, duration=Ds, optimize_time=False))
gl3.append(SRot(name='R4', theta=.1, phi=0., duration=Rs, optimize_time=False))

t0 = time.time()
syn = UnitarySynthesizer(
    primitives=gl3, subspace=sub, target=tgt,
    seed=42, optimize_times=False, optimizer='powell')
r = syn.fit(multistart=3, maxiter=400)
print(f"  obj={r.objective:.6f} success={r.success} time={time.time()-t0:.1f}s")
if hasattr(r, 'fidelity'):
    print(f"  fidelity={r.fidelity:.6f}")
