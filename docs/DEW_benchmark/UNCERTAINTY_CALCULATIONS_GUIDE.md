# Uncertainty Calculations Guide

This guide explains how to use the new C++-backed uncertainty helpers exposed through Python bindings.

Location of the implementation:
- [Reaktoro/Equilibrium/EquilibriumBenchmarkUtils.hpp](Reaktoro/Equilibrium/EquilibriumBenchmarkUtils.hpp)
- [Reaktoro/Equilibrium/EquilibriumBenchmarkUtils.cpp](Reaktoro/Equilibrium/EquilibriumBenchmarkUtils.cpp)
- [Reaktoro/Equilibrium/EquilibriumBenchmarkUtils.py.cxx](Reaktoro/Equilibrium/EquilibriumBenchmarkUtils.py.cxx)

## What Is Available

The module provides three groups of functionality:

1. Mineral database perturbation (in memory, no temp files)
- `perturbMineralDatabaseJSON(base_json, entities, shifts_j_per_mol)`
- `perturbMineralDatabase(base_json, entities, shifts_j_per_mol)`
- `perturbMineralDatabases(base_json, entities, shifts_j_per_mol_samples, num_threads=1)`

2. Species filtering and stoichiometry helpers
- `aqueousSpeciesNamesWithAllowedElements(database, allowed_elements, excluded_species=[])`
- `elementStoichiometryTerms(database, species_names, target_elements)`

3. Residual/interpolation and uncertainty bands
- `interpolateCurveValue(x, y, x_query, atol=1e-8)`
- `computeResidualsInterpolated(curve_x, curve_y, query_x, query_y, atol=1e-8)`
- `computeUncertaintyBand(samples, ci_percent=95.0)`

## Typical Full-Forward Workflow

1. Load the baseline mineral JSON text.
2. Draw Monte Carlo shifts for each entity (for example from your covariance matrix).
3. Build sampled mineral databases with `perturbMineralDatabases(...)`.
4. For each sampled database, build your `ChemicalSystem` and evaluate curves.
5. Stack sampled curve values in a matrix and call `computeUncertaintyBand(...)`.

## Python Example

```python
import json
import numpy as np
from reaktoro import *

# 1) Base mineral JSON
with open("embedded/databases/hollandpowell/tc-ds62-reaktoro.json", "r", encoding="utf-8") as f:
    base_json = f.read()

# 2) Example entities and sampled shifts [nsamples, nentities] in J/mol
entities = ["q"]
nsamples = 300
shifts = np.random.normal(loc=0.0, scale=1500.0, size=(nsamples, len(entities)))

# 3) Build sampled databases in memory (parallelizable)
databases = perturbMineralDatabases(base_json, entities, shifts, num_threads=8)

# 4) Evaluate one scalar result per sample (replace with your full curve evaluation)
values = []
for db in databases:
    # Build sampled system from db + aqueous database, solve, extract value
    # ... user code ...
    values.append(np.random.random())

# 5) Compute uncertainty band over sampled results
samples = np.asarray(values, dtype=float).reshape(-1, 1)  # [nsamples, npoints]
band = computeUncertaintyBand(samples, ci_percent=95.0)
print(band.lower, band.median, band.upper)
```

## Is It Parallelized?

Short answer: yes, partially.

Current parallelized parts:
- `perturbMineralDatabases(..., num_threads=N)` is multithreaded when `N > 1`.
- `EquilibriumSweepSolver` (separate API) supports parallel point evaluation for sweeps through `EquilibriumSweepOptions.num_threads`.

Current non-parallelized part in this utility module:
- `computeUncertaintyBand(...)` is currently single-threaded.

Important practical note:
- If your Python script still loops over sampled systems and solves one-by-one, the total workflow is only partially parallelized.
- To maximize speed, combine:
  1. `perturbMineralDatabases(..., num_threads>1)` for sampled DB creation
  2. C++ sweep APIs (`EquilibriumSweepSolver`) for batched curve evaluation

## Migration Notes from Temp-File Workflow

Previous workflow in benchmark scripts often did:
- deep copy JSON
- write temp file
- read file back into `Database`
- build sampled system

Recommended workflow now:
- perturb JSON in memory and create sampled `Database` objects directly
- avoid temp-file churn and file I/O race/permission issues

## Related Files

- Existing uncertainty benchmark script:
  [DEW_Experimental_Benchmark/Mineral_Solubilities/quartz/quartz_DEW/quartz_solubility_analysis_v2_dew24_uncertainty.py](DEW_Experimental_Benchmark/Mineral_Solubilities/quartz/quartz_DEW/quartz_solubility_analysis_v2_dew24_uncertainty.py)
- Existing sweep solver API:
  [Reaktoro/Equilibrium/EquilibriumSweepSolver.hpp](Reaktoro/Equilibrium/EquilibriumSweepSolver.hpp)
