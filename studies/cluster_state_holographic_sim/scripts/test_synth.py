"""Quick diagnostic: test UnitarySynthesizer on the cluster target."""
import sys, numpy as np, time
sys.path.insert(0, r'C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation')

from cqed_sim.unitary_synthesis import (
    Displacement as SDisp, QubitRotation as SRot, SNAP as SSNAP,
    Subspace, TargetUnitary, UnitarySynthesizer)
from cqed_sim.unitary_synthesis.targets import make_target

NCAV=8; NTR=2
U = make_target('cluster', n_match=1)
sub = Subspace.custom(NTR*NCAV, (0,1,NCAV,NCAV+1), ('|g0>','|g1>','|e0>','|e1>'))
tgt = TargetUnitary(U, ignore_global_phase=True)
Ds=48e-9; Rs=16e-9; Ss=200e-9

# Test A: D-SNAP-D (1 SNAP, no rotations) — library A style
print("=== Test A: D-SNAP-D (1 SNAP) ===")
gl_a = [
    SDisp(name='D1', alpha=.2+0j, duration=Ds, optimize_time=False),
    SSNAP(name='SNAP1', phases=[0.]*NCAV, duration=Ss, optimize_time=False),
    SDisp(name='D2', alpha=.2+0j, duration=Ds, optimize_time=False),
]
t0 = time.time()
syn = UnitarySynthesizer(primitives=gl_a, subspace=sub, target=tgt,
    seed=42, optimize_times=False, optimizer='powell')
r = syn.fit(multistart=3, maxiter=300)
print(f"  obj={r.objective:.6f} success={r.success} time={time.time()-t0:.1f}s")
print(f"  report={r.report}")
if hasattr(r, 'fidelity'):
    print(f"  fidelity={r.fidelity:.6f}")

# Test B: D-SNAP-D-SNAP-D (2 SNAP, no rotations)
print("\n=== Test B: D-SNAP-D-SNAP-D (2 SNAP) ===")
gl_b = [
    SDisp(name='D1', alpha=.2+0j, duration=Ds, optimize_time=False),
    SSNAP(name='SNAP1', phases=[0.]*NCAV, duration=Ss, optimize_time=False),
    SDisp(name='D2', alpha=.2+0j, duration=Ds, optimize_time=False),
    SSNAP(name='SNAP2', phases=[0.]*NCAV, duration=Ss, optimize_time=False),
    SDisp(name='D3', alpha=.2+0j, duration=Ds, optimize_time=False),
]
t0 = time.time()
syn = UnitarySynthesizer(primitives=gl_b, subspace=sub, target=tgt,
    seed=42, optimize_times=False, optimizer='powell')
r = syn.fit(multistart=5, maxiter=500)
print(f"  obj={r.objective:.6f} success={r.success} time={time.time()-t0:.1f}s")
print(f"  report={r.report}")
if hasattr(r, 'fidelity'):
    print(f"  fidelity={r.fidelity:.6f}")

# Test C: D-R-SNAP-D-R-SNAP-D-R (2 SNAP with rotations)
print("\n=== Test C: D-R-SNAP-D-R-SNAP-D-R (2 SNAP + rotations) ===")
gl_c = [
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
syn = UnitarySynthesizer(primitives=gl_c, subspace=sub, target=tgt,
    seed=42, optimize_times=False, optimizer='powell')
r = syn.fit(multistart=5, maxiter=500)
print(f"  obj={r.objective:.6f} success={r.success} time={time.time()-t0:.1f}s")
print(f"  report={r.report}")
if hasattr(r, 'fidelity'):
    print(f"  fidelity={r.fidelity:.6f}")

# Test D: D-R-SNAP-D-R-SNAP-D-R with nelder-mead
print("\n=== Test D: Same as C but nelder-mead ===")
t0 = time.time()
syn = UnitarySynthesizer(primitives=gl_c, subspace=sub, target=tgt,
    seed=73, optimize_times=False, optimizer='nelder-mead')
r = syn.fit(multistart=5, maxiter=1000)
print(f"  obj={r.objective:.6f} success={r.success} time={time.time()-t0:.1f}s")
print(f"  report={r.report}")
if hasattr(r, 'fidelity'):
    print(f"  fidelity={r.fidelity:.6f}")

# Print gate sequence from best result
print("\n--- Best gate sequence ---")
for g in r.sequence.gates:
    gt = type(g).__name__
    if gt == 'SNAP':
        print(f"  {g.name}: phases={[f'{p:.3f}' for p in g.phases[:4]]}...")
    elif gt == 'QubitRotation':
        print(f"  {g.name}: theta={g.theta:.4f} phi={g.phi:.4f}")
    elif gt == 'Displacement':
        print(f"  {g.name}: alpha={g.alpha:.4f}")
