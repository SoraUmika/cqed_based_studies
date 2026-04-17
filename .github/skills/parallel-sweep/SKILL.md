---
name: parallel-sweep
description: "Parallelize parameter sweeps and independent simulations for cQED studies. Use when: running multi-parameter sweeps, generating figures for multiple configurations, or any embarrassingly parallel simulation workload. Provides patterns for joblib, multiprocessing, and subagent-level parallelization."
argument-hint: "study=studies/<name> — the study with parameter sweeps to parallelize"
---

# Parallel Sweep Orchestration

## When to Use

- Parameter sweeps over multiple dimensions (e.g., chi x kappa x duration grids)
- Running the same simulation with different initial conditions or seeds
- Generating independent figures or data files
- Any computation that takes > 2 minutes and is embarrassingly parallel

## Decision Tree

```
Is the sweep embarrassingly parallel?
  (Each point independent of others, no shared state)
|
+-- YES: How many points?
|   |
|   +-- < 20 points: Use joblib.Parallel (simplest)
|   +-- 20-500 points: Use multiprocessing.Pool (better memory)
|   +-- > 500 points: Use numpy vectorization first, then parallelize remaining loops
|
+-- NO: Sequential with progress tracking
    (But check if inner loops can be vectorized)
```

## Pattern 1: joblib.Parallel (Recommended Default)

Best for: < 100 parameter combinations, moderate per-point cost.

```python
from joblib import Parallel, delayed
import numpy as np
import json
from pathlib import Path
import time

STUDY_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = STUDY_DIR / "data"

def simulate_single_point(chi, kappa, duration):
    """Run one simulation point. Must be self-contained."""
    # Import heavy modules inside the function for spawn safety
    import qutip
    # ... simulation code ...
    return {
        "chi": chi, "kappa": kappa, "duration": duration,
        "fidelity": fidelity, "leakage": leakage
    }

def main():
    # Define sweep grid
    chi_values = np.linspace(-5.0, -1.0, 10)  # MHz
    kappa_values = np.linspace(0.1, 3.0, 8)    # MHz
    duration_values = [100, 200, 400]            # ns

    param_grid = [
        (chi, kappa, dur)
        for chi in chi_values
        for kappa in kappa_values
        for dur in duration_values
    ]

    print(f"Running {len(param_grid)} parameter combinations...")
    t0 = time.time()

    # n_jobs=-1 uses all CPU cores
    results = Parallel(n_jobs=-1, verbose=10)(
        delayed(simulate_single_point)(chi, kappa, dur)
        for chi, kappa, dur in param_grid
    )

    elapsed = time.time() - t0
    print(f"Completed in {elapsed:.1f}s ({elapsed/len(param_grid):.2f}s per point)")

    # Save results
    DATA_DIR.mkdir(exist_ok=True)
    np.savez(DATA_DIR / "sweep_results.npz", results=results)

    # Log timing
    print(f"Wall time: {elapsed:.1f}s for {len(param_grid)} points")

if __name__ == "__main__":
    main()
```

## Pattern 2: multiprocessing.Pool (Better Memory for Large Sweeps)

Best for: > 100 points, or when each point uses significant memory.

```python
import multiprocessing as mp
import numpy as np
from functools import partial

def simulate_point(params, shared_config=None):
    """Worker function. Receives a tuple of parameters."""
    chi, kappa, duration = params
    # Import inside worker for Windows spawn safety
    import qutip
    # ... simulation ...
    return {"chi": chi, "kappa": kappa, "fidelity": fidelity}

def main():
    param_grid = [...]  # Build parameter grid

    # Use spawn context explicitly on Windows
    ctx = mp.get_context("spawn")

    with ctx.Pool(processes=mp.cpu_count()) as pool:
        results = pool.map(simulate_point, param_grid, chunksize=4)

    # Save and log...

if __name__ == "__main__":
    main()
```

**Windows-specific notes:**
- Always use `if __name__ == "__main__":` guard
- Import heavy modules (qutip, cqed_sim) inside worker functions, not at module level
- The `spawn` context is default on Windows and requires picklable arguments

## Pattern 3: Vectorized + Parallel Hybrid

Best for: Large grids where inner calculations can be vectorized.

```python
import numpy as np
from joblib import Parallel, delayed

def vectorized_sweep_slice(chi_array, kappa):
    """Vectorize over chi for a fixed kappa."""
    # Use numpy broadcasting instead of Python loops
    results = np.zeros(len(chi_array))
    for i, chi in enumerate(chi_array):
        # Only the truly non-vectorizable part loops
        results[i] = expensive_simulation(chi, kappa)
    return results

def main():
    chi_values = np.linspace(-5.0, -1.0, 50)
    kappa_values = np.linspace(0.1, 3.0, 20)

    # Parallelize over kappa, vectorize over chi
    all_results = Parallel(n_jobs=-1)(
        delayed(vectorized_sweep_slice)(chi_values, kappa)
        for kappa in kappa_values
    )

    results_2d = np.array(all_results)  # shape: (n_kappa, n_chi)
```

## Timing and Logging

Always log sweep timing in `IMPROVEMENTS.md` under `## Compute & Resource Notes`:

```markdown
## Compute & Resource Notes
- 2D chi-kappa sweep (240 points): 47.3s wall time, 0.20s/point, 8 cores
- Bottleneck: GRAPE optimization at large chi (max 2.8s single point)
- Speedup: 5.9x over sequential (theoretical 8x, overhead from spawn)
```

## Anti-Patterns

| Anti-Pattern | Why Bad | Do Instead |
|-------------|---------|-----------|
| Nested `Parallel` calls | Oversubscription, deadlocks | Flatten to single parallel level |
| Global state in worker functions | Race conditions, wrong results | Pass all data as arguments |
| `import qutip` at module level in workers | Slow spawn, memory bloat | Import inside worker function |
| No progress reporting for long sweeps | Can't tell if stuck or running | Use `verbose=10` or tqdm |
| Saving results inside workers | File I/O contention | Collect results, save once after |
