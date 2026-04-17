"""Check existing results."""
import json
d = json.load(open(r'c:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cqed_based_study\studies\cluster_state_holographic_sim\data\results.json'))

# Phase 4 SNAP GRAPE
p4 = d.get('p4', {})
for k, v in p4.items():
    if isinstance(v, dict) and 'sweep' in v:
        print(f'{k}: phases(first 2)={v.get("phases", [])[:2]}')
        for s in v.get('sweep', []):
            print(f'  {s["dns"]}ns: F={s["fid"]:.6f} ok={s.get("ok", "?")}')
    else:
        print(f'{k}: {v}')

# Phase 5 timing
p5 = d.get('p5', {})
print('\nPhase 5 timing:')
print(json.dumps(p5, indent=2)[:500])

# Phase 2b
p2b = d.get('p2b', {})
print('\nPhase 2b:', json.dumps(p2b, indent=2)[:200])

# Phase 6 GRAPE
p6 = d.get('p6', {})
print('\nPhase 6 GRAPE:')
for k, v in p6.items():
    print(f'  {k}: {v}')
