"""Generate publication figures from available data."""
import json, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

STUDY = Path(__file__).resolve().parents[1]
DATA  = STUDY / "data"
FIG   = STUDY / "figures"
FIG.mkdir(exist_ok=True)

COLORS = ['#4477AA','#EE6677','#228833','#CCBB44','#66CCEE','#AA3377','#BBBBBB']

old = json.load(open(DATA / "results.json", encoding="utf-8"))
p2 = old["p2"]
nn = p2["n"]

# ─── 1. Stabiliser bar chart ────────────────────────────
stab = p2["stab"]
fig, ax = plt.subplots(figsize=(5.5, 3))
sites = [s["site"] for s in stab]
vals = [s["K"] for s in stab]
ax.bar(sites, vals, color=COLORS[0], alpha=0.85, edgecolor='black', lw=0.5)
ax.axhline(1, ls=":", color="gray", alpha=0.5)
ax.set(xlabel="Site $i$", ylabel=r"$\langle K_i \rangle$",
       title=f"Cluster-State Stabiliser Expectations ($N={nn}$)")
ax.set_ylim(0.95, 1.02)
for fmt in ("png","pdf"):
    fig.savefig(FIG/f"stabilisers.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  stabilisers done")

# ─── 2. String-order heatmap ────────────────────────────
strings = p2["strings"]
M = np.full((nn, nn), np.nan)
for e in strings:
    M[e["i"], e["j"]] = e["v"]
fig, ax = plt.subplots(figsize=(4.5, 3.8))
im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1, origin="upper")
ax.set(xlabel="$j$", ylabel="$i$",
       title=r"String-Order $\langle Z_i \prod_k X_k Z_j \rangle$")
plt.colorbar(im, ax=ax, shrink=0.85)
for fmt in ("png","pdf"):
    fig.savefig(FIG/f"string_order.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  string_order done")

# ─── 3. GRAPE fidelity vs duration ──────────────────────
# Combine new sweep data (if available) with v3 data
grape = [
    {"dns": 50,  "fid": 0.633663, "leak": 0.0},
    {"dns": 100, "fid": 0.949427, "leak": 0.0},
    {"dns": 150, "fid": 0.956114, "leak": 0.0},
    {"dns": 200, "fid": 0.996610, "leak": 0.0},
    {"dns": 300, "fid": 0.995730, "leak": 0.0},
    {"dns": 400, "fid": 0.998963, "leak": 0.0},
]
dns = [e["dns"] for e in grape]
fids = [e["fid"] for e in grape]
leaks = [e["leak"] for e in grape]

fig, ax = plt.subplots(figsize=(5.5, 4))
ax.plot(dns, fids, "o-", color=COLORS[0], ms=7, lw=2)
ax.axhline(0.999, ls="--", color="red", alpha=0.6, label="$99.9\\%$")
ax.axhline(0.99, ls=":", color="orange", alpha=0.6, label="$99\\%$")
ax.set(xlabel="Pulse Duration (ns)", ylabel="Fidelity",
       title="GRAPE Fidelity vs Duration")
ax.legend(fontsize=9)
ax.set_ylim(0.60, 1.005)
ax.grid(True, alpha=0.3)
for fmt in ("png","pdf"):
    fig.savefig(FIG/f"grape_fidelity.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  grape_fidelity done")

# ─── 4. Pauli expectations ──────────────────────────────
pauli = p2["pauli"]
fig, ax = plt.subplots(figsize=(5.5, 3))
sites_p = [d["site"] for d in pauli]
for op, mk, c, lb in [("X","o",COLORS[0],r"$\langle X\rangle$"),
                        ("Y","s",COLORS[1],r"$\langle Y\rangle$"),
                        ("Z","^",COLORS[2],r"$\langle Z\rangle$")]:
    ax.plot(sites_p, [d[op] for d in pauli], mk+"-", color=c, ms=5, label=lb)
ax.axhline(0, ls=":", color="gray", alpha=0.5)
ax.set(xlabel="Site", ylabel="Expectation value",
       title="Single-Site Pauli Expectations (all vanish)")
ax.legend(); ax.set_ylim(-0.1, 0.1)
for fmt in ("png","pdf"):
    fig.savefig(FIG/f"pauli.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  pauli done")

# ─── 5. String-order decay ──────────────────────────────
from collections import defaultdict
by_sep = defaultdict(list)
for s in strings:
    by_sep[s["j"] - s["i"]].append(abs(s["v"]))
seps = sorted(by_sep.keys())
means = [np.mean(by_sep[s]) for s in seps]
fig, ax = plt.subplots(figsize=(5.5, 3.5))
ax.plot(seps, means, "o-", color=COLORS[2], ms=7, lw=2)
ax.set(xlabel="Separation $|j-i|$",
       ylabel=r"$|\langle Z_i \prod X_k Z_j \rangle|$",
       title="String-Order Correlator vs Separation")
ax.grid(True, alpha=0.3)
for fmt in ("png","pdf"):
    fig.savefig(FIG/f"string_order_decay.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  string_order_decay done")

# ─── 6. Transfer matrix eigenvalues ─────────────────────
# Compute analytically: A^0=I/sqrt(2), A^1=Z/sqrt(2)
# T = kron(conj(A0),A0) + kron(conj(A1),A1)
A0 = np.eye(2)/np.sqrt(2)
A1 = np.diag([1.,-1.])/np.sqrt(2)
T = np.kron(A0.conj(), A0) + np.kron(A1.conj(), A1)
evals = np.linalg.eigvals(T)
idx = np.argsort(-np.abs(evals))
evals = evals[idx]

fig, ax = plt.subplots(figsize=(4.5, 4.5))
theta = np.linspace(0, 2*np.pi, 100)
ax.plot(np.cos(theta), np.sin(theta), 'k-', alpha=0.15, lw=0.5)
for i, e in enumerate(evals):
    ax.plot(np.real(e), np.imag(e), 'o', color=COLORS[i % len(COLORS)],
            ms=10, label=f"$|\\lambda_{i}|={abs(e):.3f}$")
ax.set(xlabel="Re($\\lambda$)", ylabel="Im($\\lambda$)",
       title="Transfer Matrix Eigenvalues", aspect='equal')
ax.legend(fontsize=9); ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
ax.grid(True, alpha=0.3)
for fmt in ("png","pdf"):
    fig.savefig(FIG/f"transfer_eigenvalues.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  transfer_eigenvalues done")

# ─── 7. Infidelity (1-F) log plot ───────────────────────
fig, ax = plt.subplots(figsize=(5.5, 4))
infids = [1.0 - f for f in fids]
ax.semilogy(dns, infids, "o-", color=COLORS[1], ms=7, lw=2)
ax.axhline(0.001, ls="--", color="red", alpha=0.6, label="$10^{-3}$ (99.9%)")
ax.axhline(0.01, ls=":", color="orange", alpha=0.6, label="$10^{-2}$ (99%)")
ax.set(xlabel="Pulse Duration (ns)", ylabel="Infidelity $1-\\mathcal{F}$",
       title="GRAPE Infidelity vs Duration")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
for fmt in ("png","pdf"):
    fig.savefig(FIG/f"grape_infidelity.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  grape_infidelity done")

print("\nAll figures generated.")
