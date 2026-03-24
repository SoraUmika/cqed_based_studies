# Colorblind-Friendly Color Palettes

Reference for all approved color palettes used in cQED study figures.
Based on Paul Tol's color schemes (https://personal.sron.nl/~pault/).

## Tol Bright (Default cycle — 7 colors)

Used as the default `axes.prop_cycle` in the style file.

| Index | Hex       | Name      | Use Case |
|-------|-----------|-----------|----------|
| 0     | `#4477AA` | Blue      | Primary data series |
| 1     | `#EE6677` | Red       | Secondary / comparison |
| 2     | `#228833` | Green     | Third series |
| 3     | `#CCBB44` | Yellow    | Fourth series |
| 4     | `#66CCEE` | Cyan      | Fifth series |
| 5     | `#AA3377` | Purple    | Sixth series |
| 6     | `#BBBBBB` | Grey      | Reference / baseline |

```python
TOL_BRIGHT = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']
```

## Tol Vibrant (6 colors — higher contrast)

For presentations or when maximum distinction is needed.

```python
TOL_VIBRANT = ['#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988']
```

## Tol Muted (10 colors — for many-series plots)

When you need more than 7 distinguishable colors.

```python
TOL_MUTED = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933',
             '#DDCC77', '#CC6677', '#882255', '#AA4499', '#DDDDDD']
```

## Sequential Colormaps

For heatmaps, Wigner functions, and 2D parameter sweeps:

| Data Type | Colormap | Notes |
|-----------|----------|-------|
| Positive-only (population, fidelity) | `viridis` | Perceptually uniform, safe |
| Diverging (Wigner, detuning) | `RdBu_r` | Red-Blue diverging, centered on white |
| Phase | `twilight` | Cyclic colormap for 0→2π data |

```python
# Wigner function example
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
im = ax.pcolormesh(x, p, W, cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
fig.colorbar(im, ax=ax, label='$W(x, p)$')
```

## Quick Copy-Paste

```python
# At the top of any plotting script:
TOL_BRIGHT  = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']
TOL_VIBRANT = ['#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988']
TOL_MUTED   = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933',
               '#DDCC77', '#CC6677', '#882255', '#AA4499', '#DDDDDD']
```
