import json
d = json.load(open('data/results.json', encoding='utf-8'))

p1 = d['p1']
print('=== PHASE 1 ===')
for k,v in p1.items():
    if isinstance(v, list) and len(v) > 4:
        print(f"  {k}: (array len {len(v)})")
    else:
        print(f"  {k}: {v}")
print()

p2 = d['p2']
print('=== PHASE 2 ===')
print('n sites:', p2['n'])
print('Stabilisers (Ki):')
for s in p2['stab']:
    print(f"  site {s['site']}: K={s['K']}")
print('Pauli:')
for p in p2['pauli']:
    print(f"  site {p['site']}: X={p['X']:.2e} Y={p['Y']:.2e} Z={p['Z']:.2e}")
print('String-order:')
for s in p2['strings']:
    print(f"  ({s['i']},{s['j']}): v={s['v']:.6f}")
print()

p2b = d['p2b']
print('=== PHASE 2b (Holographic) ===')
for k,v in p2b.items():
    if isinstance(v, list) and len(v) > 4:
        print(f"  {k}: (array len {len(v)})")
    else:
        print(f"  {k}: {v}")
print()

p6 = d['p6']
print('=== PHASE 6 (GRAPE) ===')
print(type(p6), p6)
