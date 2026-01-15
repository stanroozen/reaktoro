# DEW Integration Strategy for Reaktoro Minimizer

## Overview

The goal is to integrate the DEW (Deep Earth Water) model into Reaktoro as a proper extension that can be used with the minimizer and `ChemicalSystem` framework, similar to how PHREEQC (port of PREEQC) was integrated.

## Current State

- **DEW exists as:** A separate extension with water models, databases, and test suite
- **Problem:** DEW models are accessed directly but not integrated with Reaktoro's `ActivityModel` interface
- **Root Cause of Quartz Script Failure:** The Python script tried to use `ActivityModel` with water properties directly, but `ActivityModel` expects a callable function that takes `ActivityModelArgs` and returns `ActivityProps`

## Architecture Overview

### Reaktoro's Activity Model System

```
ActivityModelGenerator (SpeciesList) → ActivityModel
                                         ↓
                                      ActivityModel(ActivityModelArgs) → ActivityProps

ActivityModelArgs contains:
  - T (temperature, K)
  - P (pressure, Pa)
  - x (mole fractions)

ActivityProps contains:
  - ln_g (activity coefficients)
  - ln_a (activities)
  - som (state of matter)
  - extra (shared data storage)
```

### Implementation Pattern (from DebyeHuckel)

1. **Generator Function:** Takes `SpeciesList` and returns an `ActivityModel`
2. **Captured Data:** Mineral parameters, species information cached in closures
3. **Lambda Function:** The actual `ActivityModel` that:
   - Takes T, P, x from `ActivityModelArgs`
   - Computes water properties via DEW models
   - Calculates species activities and coefficients
   - Stores intermediate results in `props.extra` for reuse

## Required New Files/Components

### 1. **DEWActivityModel.hpp** (NEW)
Location: `Reaktoro/Extensions/DEW/DEWActivityModel.hpp`

```cpp
// Forward declarations
class DEWDatabase;
struct WaterState;
class WaterThermoModel;
class WaterElectroModel;

// Main generator function
auto ActivityModelDEW(
    const SpeciesList& species,
    const DEWDatabase& database) -> ActivityModelGenerator;

// Optional: Specialized generators
auto ActivityModelDEW(
    const SpeciesList& species,
    const DEWDatabase& database,
    const WaterModelOptions& options) -> ActivityModelGenerator;
```

### 2. **DEWActivityModel.cpp** (NEW)
Location: `Reaktoro/Extensions/DEW/DEWActivityModel.cpp`

**Key Components:**
- `AqueousMixture`-like class that adapts DEW species to Reaktoro species
- Lambda function implementing the activity model calculation
- Integration with `WaterThermoModel`, `WaterElectroModel`, `WaterGibbsModel`
- Caching strategy for water properties in `props.extra`

### 3. **Update DEW.hpp** (MODIFY)
Add export of new `ActivityModelDEW` generator

### 4. **DEWActivityModel.py.cxx** (NEW)
Python bindings for the activity model

### 5. **Tests: DEWActivityModel.test.cxx** (NEW)
Unit tests for the new activity model with:
- Simple aqueous systems
- Quartz saturation calculations
- Comparison with direct water model results

## Implementation Steps

### Phase 1: Core Implementation

1. **Create DEWActivityModel.hpp**
   - Define the `ActivityModelDEW()` generator function
   - Design internal data structures for caching

2. **Create DEWActivityModel.cpp**
   - Implement the activity model calculation
   - Integrate with existing water models (ThermoModel, ElectroModel)
   - Handle species-to-DEW database mapping
   - Calculate ionic strength and activity coefficients

3. **Update CMakeLists.txt**
   - Add new .cpp and .test.cxx files to build

4. **Update DEW.hpp**
   - Include and export `DEWActivityModel.hpp`

### Phase 2: Python Integration

1. **Create DEWActivityModel.py.cxx**
   - Expose `ActivityModelDEW()` to Python
   - Allow Python users to create activity models

2. **Update Extensions/DEW.py.cxx**
   - Include bindings for new C++ functions

### Phase 3: Testing & Validation

1. **Create unit tests** in `DEWActivityModel.test.cxx`
   - Test with simple systems
   - Verify activity calculations

2. **Update quartz_solubility_analysis.py**
   - Use new `ActivityModelDEW` generator
   - Test with experimental data

## Key Design Decisions

### 1. Water Properties Caching

Use `props.extra` dict to cache expensive water property calculations:
```cpp
props.extra["WaterThermoProps"] = std::make_shared<WaterThermoProps>(...);
props.extra["WaterElectroProps"] = std::make_shared<WaterElectroProps>(...);
props.extra["WaterState"] = std::make_shared<WaterState>(...);
```

### 2. Species Mapping

- Map Reaktoro's `SpeciesList` to DEW database species
- Handle missing species gracefully (fallback models)
- Cache species indices and properties

### 3. Activity Calculation

For aqueous mixtures:
- Calculate ionic strength from solution composition
- Compute water properties (density, dielectric constant, etc.) using DEW models
- Apply HKF model with Born term corrections
- Handle temperature and pressure dependencies

### 4. Compatibility

Ensure `ActivityModelDEW` works with:
- `ChemicalSystem` constructor
- `EquilibriumSolver`
- `Minimizer` framework
- Python interface

## Modified Python Script

After implementation, `quartz_solubility_analysis.py` will be:

```python
from reaktoro import *

# Load databases
db = SupcrtDatabase("supcrtbl")
dew_db = DEWDatabase("dew2024-aqueous")

# Create chemical system with DEW activity model
system = ChemicalSystem(
    db,
    AqueousPhase("H2O(aq) H+ OH- SiO2(aq)").setActivityModel(ActivityModelDEW(dew_db)),
    MineralPhases("Quartz")
)

# Now equilibrium calculations will use DEW for water properties
solver = EquilibriumSolver(system)
# ... rest of script works
```

## Expected Behavior After Integration

1. **Direct integration with Reaktoro's minimizer** - No manual water property handling needed
2. **Proper activity calculations** - ActivityModel interface ensures correct format
3. **Scalability** - Works with any number of aqueous species
4. **Extensibility** - Can chain with other activity models (e.g., ionic strength correction)
5. **Performance** - Water properties cached and reused across species calculations

## Testing Strategy

### Unit Tests (DEWActivityModel.test.cxx)
- Activity coefficients at various T, P
- Ionic strength calculations
- Osmotic coefficient tests
- Comparison with direct model results

### Integration Tests (quartz_solubility_analysis.py)
- End-to-end quartz solubility calculations
- Comparison with experimental data from CSV
- Multiple pressure-temperature conditions
- Verification of equilibrium convergence

## Fallback Strategy

If full HKF implementation is complex:
1. Start with simplified activity model (Born term only)
2. Gradually add corrections (ionic strength, salting-out)
3. Integrate full HKF equation later

## Success Criteria

✓ Activity model compiles and links with Reaktoro
✓ Python bindings work
✓ Quartz solubility script runs without errors
✓ Calculated activities match direct water model results
✓ Equilibrium calculations converge
✓ Results agree with experimental data within uncertainty
✓ No crashes or memory leaks
✓ Works with `ChemicalSystem` and `Minimizer`
