"""Test synthesis using entangling primitives (SQR, CPSQR) for the cluster target."""
import sys, numpy as np, time
sys.path.insert(0, r'C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation')

from cqed_sim.unitary_synthesis import (
    Displacement as SDisp, QubitRotation as SRot, SNAP as SSNAP,
    SQR, ConditionalPhaseSQR as CPSQR,
    Subspace, TargetUnitary, UnitarySynthesizer)
from cqed_sim.unitary_synthesis.targets import make_target

NCAV=8; NTR=2
U = make_target('cluster', n_match=1)
sub = Subspace.custom(NTR*NCAV, (0,1,NCAV,NCAV+1),
                      ('|g0>','|g1>','|e0>','|e1>'))
tgt = TargetUnitary(U, ignore_global_phase=True)
Ds=48e-9; Rs=16e-9; Ss=200e-9

# Test 1: CPSQR + SNAP + D + R  (2 layers)
print("=== Test 1: R-SNAP-CPSQR-D-R-SNAP-CPSQR-D (2 selective layers) ===")
gl1 = [
    SRot(name='R1', theta=np.pi/2, phi=0., duration=Rs, optimize_time=False),
    SSNAP(name='SNAP1', phases=[0.]*NCAV, duration=Ss, optimize_time=False),
    CPSQR(name='CP1', duration=Ss, optimize_time=False, phases_n=[0.]*NCAV),
    SDisp(name='D1', alpha=.2+0j, duration=Ds, optimize_time=False),
    SRot(name='R2', theta=np.pi/2, phi=0., duration=Rs, optimize_time=False),
    SSNAP(name='SNAP2', phases=[0.]*NCAV, duration=Ss, optimize_time=False),
    CPSQR(name='CP2', duration=Ss, optimize_time=False, phases_n=[0.]*NCAV),
    SDisp(name='D2', alpha=.2+0j, duration=Ds, optimize_time=False),
    SRot(name='R3', theta=0.1, phi=0., duration=Rs, optimize_time=False),
]
t0 = time.time()
syn = UnitarySynthesizer(primitives=gl1, subspace=sub, target=tgt,
    seed=42, optimize_times=False, optimizer='powell')
r = syn.fit(multistart=5, maxiter=500)
fid1 = r.report.get('metrics',{}).get('fidelity', 1-r.objective)
print(f"  obj={r.objective:.6f} fidelity={fid1:.6f} time={time.time()-t0:.1f}s")

# Test 2: SQR + D + R (library B style)
print("\n=== Test 2: R-SQR-R-D-R-SQR-R (library B style) ===")
gl2 = [
    SRot(name='R1', theta=np.pi/2, phi=0., duration=Rs, optimize_time=False),
    SQR(name='S1', theta_n=[0.1]*NCAV, phi_n=[0.]*NCAV, duration=Ss, optimize_time=False),
    SRot(name='R2', theta=np.pi/2, phi=np.pi/2, duration=Rs, optimize_time=False),
    SDisp(name='D1', alpha=.2+0j, duration=Ds, optimize_time=False),
    SRot(name='R3', theta=np.pi/2, phi=0., duration=Rs, optimize_time=False),
    SQR(name='S2', theta_n=[0.1]*NCAV, phi_n=[0.]*NCAV, duration=Ss, optimize_time=False),
    SRot(name='R4', theta=np.pi/2, phi=np.pi/2, duration=Rs, optimize_time=False),
]
t0 = time.time()
syn = UnitarySynthesizer(primitives=gl2, subspace=sub, target=tgt,
    seed=42, optimize_times=False, optimizer='powell')
r2 = syn.fit(multistart=5, maxiter=500)
fid2 = r2.report.get('metrics',{}).get('fidelity', 1-r2.objective)
print(f"  obj={r2.objective:.6f} fidelity={fid2:.6f} time={time.time()-t0:.1f}s")

# Test 3: SNAP + CPSQR + D + R (using both SNAP for cavity phase + CPSQR for qubit-conditional)
print("\n=== Test 3: D-SNAP-CPSQR-R (x2 layers, simpler structure) ===")
gl3 = [
    SDisp(name='D1', alpha=.2+0j, duration=Ds, optimize_time=False),
    SSNAP(name='SNAP1', phases=[0.]*NCAV, duration=Ss, optimize_time=False),
    CPSQR(name='CP1', duration=Ss, optimize_time=False, phases_n=[0.]*NCAV),
    SRot(name='R1', theta=0.1, phi=0., duration=Rs, optimize_time=False),
    SDisp(name='D2', alpha=.2+0j, duration=Ds, optimize_time=False),
    SSNAP(name='SNAP2', phases=[0.]*NCAV, duration=Ss, optimize_time=False),
    CPSQR(name='CP2', duration=Ss, optimize_time=False, phases_n=[0.]*NCAV),
    SRot(name='R2', theta=0.1, phi=0., duration=Rs, optimize_time=False),
    SDisp(name='D3', alpha=.2+0j, duration=Ds, optimize_time=False),
]
t0 = time.time()
syn = UnitarySynthesizer(primitives=gl3, subspace=sub, target=tgt,
    seed=42, optimize_times=False, optimizer='powell')
r3 = syn.fit(multistart=5, maxiter=500)
fid3 = r3.report.get('metrics',{}).get('fidelity', 1-r3.objective)
print(f"  obj={r3.objective:.6f} fidelity={fid3:.6f} time={time.time()-t0:.1f}s")

# Print best gate parameters
best_r = min([(r,fid1,'T1'),(r2,fid2,'T2'),(r3,fid3,'T3')], key=lambda x:-x[1])
print(f"\nBest: {best_r[2]} with F={best_r[1]:.6f}")
for g in best_r[0].sequence.gates:
    gt = type(g).__name__
    if gt == 'SNAP':
        print(f"  {g.name}: phases={[f'{p:.4f}' for p in g.phases[:3]]}...")
    elif gt == 'ConditionalPhaseSQR':
        print(f"  {g.name}: phases_n={[f'{p:.4f}' for p in g.phases_n[:3]]}...")
    elif gt == 'SQR':
        print(f"  {g.name}: theta={[f'{t:.4f}' for t in g.theta_n[:3]]}...")
    elif gt == 'QubitRotation':
        print(f"  {g.name}: theta={g.theta:.4f} phi={g.phi:.4f}")
    elif gt == 'Displacement':
        print(f"  {g.name}: alpha={g.alpha:.4f}")
