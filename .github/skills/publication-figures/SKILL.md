---
name: publication-figures
description: "Create consistent, publication-quality matplotlib plots. Use when: generating figures, making plots, styling graphs, creating visualizations, saving figures for LaTeX reports. Provides a shared matplotlib style and colorblind-friendly palettes."
argument-hint: "What to plot, e.g. 'fidelity vs drive amplitude sweep'"
---

# Publication-Quality Figures

## When to Use

- Generating any plot or figure during a study (AGENTS.md Step 3)
- User asks for "publication-quality" or "paper-ready" figures
- Styling existing plots for consistency across studies

## Procedure

### 1. Apply the Style

Load the shared matplotlib style at the top of every plotting script:

```python
import matplotlib.pyplot as plt

plt.style.use('.github/skills/publication-figures/assets/cqed_style.mplstyle')
```

Or if running from a study's `scripts/` directory:

```python
from pathlib import Path

style_path = Path(__file__).resolve().parents[2] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
plt.style.use(str(style_path))
```

### 2. Use Colorblind-Friendly Colors

See the [color palettes reference](./references/color_palettes.md) for approved color cycles.

The default style file already sets a colorblind-friendly cycle (Tol's Bright). For manual color selection, import the palettes:

```python
# Tol's Bright palette
TOL_BRIGHT = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']

# Use for specific line assignments
fig, ax = plt.subplots()
ax.plot(x, y1, color=TOL_BRIGHT[0], label='Fidelity')
ax.plot(x, y2, color=TOL_BRIGHT[1], label='Leakage')
```

### 3. Label Everything

Every figure **must** have:

- **Axis labels** with units (e.g., `$\\omega_d / 2\\pi$ (GHz)`, `Time (ns)`, `Fidelity`)
- **Legend** if multiple datasets are shown
- **Title** only if the caption won't provide context (omit for LaTeX reports)

### 4. Save in Both Formats

Always save both raster and vector formats:

```python
fig.savefig('figures/fidelity_vs_drive.png', dpi=300, bbox_inches='tight')
fig.savefig('figures/fidelity_vs_drive.pdf', bbox_inches='tight')
plt.close(fig)
```

- `.png` at 300 DPI — for README and quick inspection
- `.pdf` — for LaTeX report inclusion (vector, no quality loss)

### 5. Common Plot Types for cQED Studies

| Plot Type | Use Case | Key Settings |
|-----------|----------|-------------|
| Fidelity sweep | OPT, ANA | x = parameter, y = fidelity, horizontal dashed line at target |
| Time evolution | DES, REP | x = time (ns), y = population/expectation value |
| Spectrum | ANA, REP | x = frequency (GHz), y = amplitude or transmission |
| Wigner function | DES, ANA | 2D colormap, use diverging colormap (`RdBu_r`) |
| Convergence | Validation | x = parameter value, y = metric, log scale for y-axis |
| Phase diagram | ANA | 2D colormap, labeled contours for operating points |

### 6. Paired Pulse Plots for OPT/DES Studies

If the study optimizes a pulse or waveform, generate both of these appendix-ready figures:

- **Time domain** — I/Q or amplitude/phase versus time, with units.
- **Frequency domain** — magnitude spectrum versus frequency, with units.

Use consistent styling and save both `.png` and `.pdf` versions. Future agents should be able to verify at a glance that the pulse is physically reasonable and not just numerically high fidelity.

## Rules

- Never use the default matplotlib color cycle — always use the style file or Tol palettes.
- Always include units on axes.
- Always save both `.png` and `.pdf`.
- Always call `plt.close(fig)` after saving to avoid memory leaks in batch scripts.
- Use `bbox_inches='tight'` to avoid clipped labels.
- For waveform-optimization studies, paired time-domain and frequency-domain plots are mandatory, not optional.
