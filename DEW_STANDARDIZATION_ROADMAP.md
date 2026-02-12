# DEW Model Standardization Roadmap for Reaktoro

## Executive Summary

This document outlines the required work to fully standardize the Deep Earth Water (DEW) model implementation in Reaktoro, making it seamlessly compatible with all Reaktoro calculation types using standard API calls.

**Current Status:** ⚠️ Partially Integrated
- ✅ DEW thermodynamic properties work
- ✅ Basic equilibrium calculations functional
- ❌ Not standardized for all Reaktoro features
- ❌ API inconsistent with other models (HKF, Pitzer)

**Goal:** Make DEW a first-class citizen in Reaktoro with the same ease of use as HKF/Pitzer models.

---

## PHASE 1: Core API Standardization 🔴 CRITICAL

### 1.1 Database Constructors
**Priority:** CRITICAL | **Effort:** 2-3 days

**Current Problem:**
```python
# This crashes in Jupyter notebooks:
db = DEW2024()

# This works but is inconsistent:
db = DEWDatabase("dew2024-aqueous")
```

**Required Changes:**

#### File: `Reaktoro/Core/DEW.hpp`
```cpp
// Add safe constructor wrapper
class DEW2024 {
public:
    DEW2024() {
        try {
            // Initialize with proper error handling
            this->db_ = DEWDatabase("dew2024-aqueous");
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "DEW2024 initialization failed. "
                "Database files may be missing or corrupted: " +
                std::string(e.what())
            );
        }
    }

    // Ensure proper cleanup
    ~DEW2024() noexcept {
        // Safe cleanup code
    }

private:
    Database db_;
};
```

#### File: `Reaktoro/python/PyDEW.cpp`
```cpp
// Fix Python bindings
void exportDEW(py::module& m) {
    py::class_<DEW2024>(m, "DEW2024")
        .def(py::init([]() {
            try {
                return std::make_unique<DEW2024>();
            } catch (const std::exception& e) {
                throw py::value_error(std::string("DEW2024: ") + e.what());
            }
        }), "Create DEW2024 database with safe initialization")
        .def("__repr__", [](const DEW2024&) {
            return "<DEW2024 Database (dew2024-aqueous)>";
        });
}
```

**Testing:**
- ✅ Works in Python scripts
- ✅ Works in Jupyter notebooks
- ✅ Works in VS Code interactive window
- ✅ Proper error messages on failure

---

### 1.2 Activity Model API Consistency
**Priority:** CRITICAL | **Effort:** 3-5 days

**Current Problem:**
```python
# HKF is simple:
aqueous.setActivityModel(ActivityModelHKF())

# DEW requires configuration:
aqueous.setActivityModel(ActivityModelDEW())  # Uses defaults, not obvious

# Should match this pattern:
aqueous.setActivityModel(ActivityModelPitzer())
```

**Required Changes:**

#### File: `Reaktoro/Thermodynamics/Models/ActivityModelDEW.hpp`
```cpp
/// Options for DEW activity model
struct ActivityModelDEWOptions {
    WaterEosModel eosModel = WaterEosModel::ZhangDuan2005;
    WaterDielectricModel dielectricModel = WaterDielectricModel::PowerFunction;
    WaterGibbsModel gibbsModel = WaterGibbsModel::DewIntegral;
    WaterBornModel bornModel = WaterBornModel::Shock92Dew;

    /// Builder pattern for easy configuration
    auto withEOS(WaterEosModel model) -> ActivityModelDEWOptions& {
        eosModel = model;
        return *this;
    }

    auto withDielectric(WaterDielectricModel model) -> ActivityModelDEWOptions& {
        dielectricModel = model;
        return *this;
    }

    // ... similar for other options
};

/// Standard constructor with sensible defaults
auto ActivityModelDEW() -> ActivityModel;

/// Constructor with custom options
auto ActivityModelDEW(ActivityModelDEWOptions options) -> ActivityModel;
```

**Python API:**
```python
# Simple default usage (like HKF)
aqueous.setActivityModel(ActivityModelDEW())

# Advanced usage with builder pattern
from reaktoro import ActivityModelDEW, ActivityModelDEWOptions

options = (ActivityModelDEWOptions()
    .withEOS(WaterEosModel.ZhangDuan2009)
    .withDielectric(WaterDielectricModel.JohnsonNorton1991))

aqueous.setActivityModel(ActivityModelDEW(options))

# Or inline:
aqueous.setActivityModel(
    ActivityModelDEW()
    .withEOS(WaterEosModel.ZhangDuan2009)
)
```

---

### 1.3 Standard ChemicalSystem Creation
**Priority:** HIGH | **Effort:** 1-2 days

**Goal:** Make DEW system creation identical to HKF/Pitzer pattern

**Standard Pattern (HKF):**
```python
db = SupcrtDatabase("supcrt98")
aqueous = AqueousPhase("H2O(aq) H+ OH- Na+ Cl- SiO2(aq)")
aqueous.setActivityModel(ActivityModelHKF())
system = ChemicalSystem(db, aqueous)
```

**Required DEW Pattern:**
```python
# Should be this simple:
db = DEW2024()  # Fixed in 1.1
aqueous = AqueousPhase("WATER,AQ H+ OH- SiO2_aq")
aqueous.setActivityModel(ActivityModelDEW())  # Fixed in 1.2
system = ChemicalSystem(db, aqueous)
```

**Implementation:**
- Ensure species naming is consistent
- Auto-detect DEW species format
- Provide clear error messages for species mismatches

---

## PHASE 2: Full Reaktoro Feature Support 🟡 HIGH PRIORITY

### 2.1 Equilibrium Calculations
**Priority:** HIGH | **Effort:** 2-3 days

**Required Support:**

#### Basic Equilibrium
```python
# Should work identically to HKF
specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()

solver = EquilibriumSolver(specs)
state = ChemicalState(system)
state.set("H2O(aq)", 1.0, "kg")

conditions = EquilibriumConditions(specs)
conditions.temperature(300, "celsius")
conditions.pressure(1000, "bar")

result = solver.solve(state, conditions)  # Must work with DEW
```

**Testing Matrix:**
- ✅ Temperature range: 25-1000°C
- ✅ Pressure range: 1-100,000 bar
- ✅ All charge balance modes
- ✅ pH specifications
- ✅ Element constraints
- ✅ Activity constraints

#### File: `tests/test_dew_equilibrium.py`
```python
import pytest
from reaktoro import *

def test_dew_equilibrium_basic():
    """Test basic equilibrium with DEW"""
    db = DEW2024()
    aqueous = AqueousPhase("WATER,AQ H+ OH-")
    aqueous.setActivityModel(ActivityModelDEW())
    system = ChemicalSystem(db, aqueous)

    state = ChemicalState(system)
    state.set("H2O(aq)", 1.0, "kg")
    state.temperature(300, "celsius")
    state.pressure(1000, "bar")

    solver = EquilibriumSolver(system)
    result = solver.solve(state)

    assert result.succeeded()
    assert state.temperature() == pytest.approx(300 + 273.15)

def test_dew_equilibrium_high_PT():
    """Test DEW at extreme conditions"""
    # Test at 800°C, 10 kbar
    # ...
```

---

### 2.2 Speciation Calculations
**Priority:** HIGH | **Effort:** 1-2 days

**Current Status:** Partially working but not documented

**Required:**
```python
# Should work exactly like HKF speciation
system = ChemicalSystem(DEW2024(), phases)

# Calculate speciation at high T/P
props = AqueousProps(system)
props.update(state)

# All properties should work:
pH = props.pH()
pe = props.pE()
Eh = props.Eh()
alkalinity = props.alkalinity()
ionicStrength = props.ionicStrength()

# Species activities
activity_H2O = props.speciesActivity("H2O(aq)")
activity_SiO2 = props.speciesActivity("SiO2(aq)")
```

**Testing:** Create test suite for all AqueousProps methods with DEW

---

### 2.3 Solubility Calculations
**Priority:** HIGH | **Effort:** 2-3 days

**Goal:** DEW should work for mineral solubility just like HKF

**Standard Pattern:**
```python
# Create system with DEW + minerals
db_dew = DEW2024()
db_minerals = SUPCRTBL()
db = Database(db_dew.species())
db.addSpecies(db_minerals.species("Quartz"))

aqueous = AqueousPhase("WATER,AQ H+ OH- SiO2_aq")
aqueous.setActivityModel(ActivityModelDEW())
mineral = MineralPhase("Quartz")

system = ChemicalSystem(db, aqueous, mineral)

# Calculate solubility (should be identical to HKF pattern)
state = ChemicalState(system)
state.set("H2O(aq)", 1.0, "kg")
state.set("Quartz", 10.0, "mol")

# Solve at multiple T/P
for T in [200, 300, 400, 500]:
    for P in [500, 1000, 5000]:
        state.temperature(T, "celsius")
        state.pressure(P, "bar")
        solver.solve(state)

        props = AqueousProps(state)
        solubility = props.speciesMolality("SiO2(aq)")
        print(f"T={T}°C, P={P}bar: {solubility:.6f} mol/kg")
```

**Testing:**
- ✅ Quartz solubility vs literature
- ✅ Calcite solubility
- ✅ Other common minerals
- ✅ Comparison with experimental data

---

### 2.4 Reaction Path Modeling
**Priority:** MEDIUM | **Effort:** 3-5 days

**Goal:** DEW should work with EquilibriumPath

**Required:**
```python
# Reaction path with DEW
system = ChemicalSystem(DEW2024(), phases)

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()

path = EquilibriumPath(specs)

# Path should work at high T/P
initial_state = state_at_25C_1bar
final_state = state_at_800C_10kbar

path_result = path.solve(initial_state, final_state)

# Should produce valid path with DEW thermodynamics
```

**Testing:**
- ✅ Temperature paths
- ✅ Pressure paths
- ✅ Fluid-rock interaction paths
- ✅ Evaporation/dilution paths

---

### 2.5 Kinetic Modeling
**Priority:** MEDIUM | **Effort:** 5-7 days

**Goal:** DEW-based equilibrium used in kinetic calculations

**Required:**
```python
# Kinetics with DEW thermodynamics
system = ChemicalSystem(DEW2024(), phases)

kinetics = EquilibriumKinetics(system)
kinetics.add(MineralReaction("Quartz")
    .setRateModel(mineralRateModel)
    .setSurfaceArea(1.0, "m2"))

# Should calculate kinetic path with DEW
solver = KineticsSolver(kinetics)
result = solver.solve(initial_state, dt=1.0)
```

---

### 2.6 Reactive Transport
**Priority:** LOW | **Effort:** 7-10 days

**Goal:** Use DEW in reactive transport simulations

**Required:**
- Integration with transport solvers
- DEW properties in flow calculations
- Multi-phase transport with DEW

---

## PHASE 3: Developer Experience 🟢 MEDIUM PRIORITY

### 3.1 Comprehensive Examples
**Priority:** MEDIUM | **Effort:** 3-5 days

**Create examples for:**
- [ ] Basic equilibrium with DEW
- [ ] Mineral solubility (quartz, calcite, etc.)
- [ ] pH/speciation calculations
- [ ] High T/P geochemistry
- [ ] Fluid-rock interaction
- [ ] Comparison with HKF model
- [ ] Parameter sensitivity analysis

**Example Structure:**
```
examples/
├── dew/
│   ├── 01-basic-equilibrium.py
│   ├── 02-mineral-solubility.py
│   ├── 03-speciation-high-TP.py
│   ├── 04-fluid-rock-interaction.py
│   ├── 05-comparison-hkf-dew.py
│   └── notebooks/
│       ├── quartz-solubility-tutorial.ipynb
│       ├── fluid-inclusions.ipynb
│       └── deep-earth-fluids.ipynb
```

---

### 3.2 Documentation
**Priority:** MEDIUM | **Effort:** 5-7 days

**Required Documentation:**

#### User Guide
- [ ] "Getting Started with DEW" tutorial
- [ ] DEW model theory overview
- [ ] Parameter selection guide
- [ ] When to use DEW vs HKF
- [ ] Limitations and valid ranges

#### API Reference
- [ ] Complete DEW class documentation
- [ ] ActivityModelDEW options
- [ ] Water property models (EOS, dielectric, etc.)
- [ ] Species naming conventions

#### Developer Guide
- [ ] DEW implementation details
- [ ] How to extend DEW
- [ ] Adding new water models
- [ ] Contributing to DEW

---

### 3.3 Testing Suite
**Priority:** HIGH | **Effort:** 5-7 days

**Comprehensive Test Coverage:**

#### Unit Tests
```
tests/
├── test_dew_database.py          # Database loading/initialization
├── test_dew_activity_models.py   # Activity model calculations
├── test_dew_water_properties.py  # Water EOS, dielectric, etc.
├── test_dew_equilibrium.py       # Basic equilibrium
├── test_dew_solubility.py        # Mineral solubility
└── test_dew_integration.py       # Full workflow tests
```

#### Benchmark Tests
```
benchmarks/
├── dew_vs_experimental_data.py   # Compare with literature
├── dew_vs_hkf.py                 # Compare models
└── dew_performance.py            # Speed benchmarks
```

#### Validation Tests
- [ ] Reproduce published DEW results
- [ ] Match experimental quartz solubility
- [ ] Match water density/dielectric data
- [ ] Compare with SUPCRT at overlapping conditions

---

## PHASE 4: Advanced Features 🔵 LOW PRIORITY

### 4.1 Automatic Database Selection
**Priority:** LOW | **Effort:** 2-3 days

**Goal:** Smart database selection based on conditions

```python
# Automatic selection
system = ChemicalSystem.auto(
    species=["H2O(aq)", "SiO2(aq)", "Quartz"],
    temperature=800,  # °C
    pressure=10000    # bar
)
# Automatically uses DEW for high T/P
```

---

### 4.2 Hybrid Models
**Priority:** LOW | **Effort:** 5-7 days

**Goal:** Combine DEW with other models

```python
# Use DEW for water, HKF for other species
aqueous = AqueousPhase(...)
aqueous.setActivityModel(
    ActivityModelHybrid()
        .withDEW(["H2O(aq)"])
        .withHKF(["Na+", "Cl-", "SiO2(aq)"])
)
```

---

### 4.3 Parameter Optimization
**Priority:** LOW | **Effort:** 7-10 days

**Goal:** Fit DEW parameters to experimental data

```python
# Optimize DEW parameters
optimizer = DEWParameterOptimizer()
optimizer.addData(experimental_data)
optimizer.addConstraints(thermodynamic_constraints)

optimized_params = optimizer.fit()
```

---

## PHASE 5: Quality Assurance ✅ ONGOING

### 5.1 Continuous Integration
**Priority:** HIGH | **Effort:** Ongoing

**Required CI/CD:**
- [ ] Run all DEW tests on every commit
- [ ] Test on Windows/Linux/macOS
- [ ] Test with Python 3.8, 3.9, 3.10, 3.11, 3.12
- [ ] Test in Jupyter notebooks
- [ ] Performance regression testing

---

### 5.2 Code Quality
**Priority:** MEDIUM | **Effort:** Ongoing

**Standards:**
- [ ] Consistent code style (clang-format)
- [ ] Complete docstrings
- [ ] Type hints in Python
- [ ] Memory leak checks
- [ ] Static analysis (cppcheck, pylint)

---

### 5.3 Performance Optimization
**Priority:** MEDIUM | **Effort:** 5-7 days

**Optimization Targets:**
- [ ] Cache water properties for repeated T/P
- [ ] Vectorized calculations where possible
- [ ] Parallel equilibrium solving for multiple states
- [ ] Profile and optimize hot paths

**Target Performance:**
- DEW equilibrium solving < 2× HKF time
- Batch calculations use parallel execution
- No memory leaks in long-running calculations

---

## Implementation Priority Matrix

| Phase | Component | Priority | Effort | Impact | Dependencies |
|-------|-----------|----------|--------|--------|--------------|
| 1.1 | Database constructors | 🔴 CRITICAL | 2-3 days | HIGH | None |
| 1.2 | Activity model API | 🔴 CRITICAL | 3-5 days | HIGH | 1.1 |
| 1.3 | ChemicalSystem creation | 🟡 HIGH | 1-2 days | HIGH | 1.1, 1.2 |
| 2.1 | Equilibrium calculations | 🟡 HIGH | 2-3 days | HIGH | 1.1, 1.2 |
| 2.2 | Speciation | 🟡 HIGH | 1-2 days | MEDIUM | 2.1 |
| 2.3 | Solubility | 🟡 HIGH | 2-3 days | HIGH | 2.1 |
| 3.3 | Testing suite | 🟡 HIGH | 5-7 days | CRITICAL | All above |
| 3.1 | Examples | 🟢 MEDIUM | 3-5 days | HIGH | 2.1-2.3 |
| 3.2 | Documentation | 🟢 MEDIUM | 5-7 days | HIGH | All above |
| 2.4 | Reaction paths | 🟢 MEDIUM | 3-5 days | MEDIUM | 2.1 |
| 2.5 | Kinetics | 🟢 MEDIUM | 5-7 days | MEDIUM | 2.1, 2.4 |
| 4.1-4.3 | Advanced features | 🔵 LOW | 15-20 days | LOW | All above |
| 2.6 | Reactive transport | 🔵 LOW | 7-10 days | LOW | All above |

---

## Estimated Timeline

### Sprint 1 (2 weeks): Core Stability
- ✅ Database constructors (1.1)
- ✅ Activity model API (1.2)
- ✅ ChemicalSystem creation (1.3)
- ✅ Basic equilibrium tests

### Sprint 2 (2 weeks): Essential Features
- ✅ Full equilibrium support (2.1)
- ✅ Speciation calculations (2.2)
- ✅ Solubility calculations (2.3)
- ✅ Test suite foundation (3.3)

### Sprint 3 (2 weeks): Polish & Documentation
- ✅ Comprehensive examples (3.1)
- ✅ Complete documentation (3.2)
- ✅ Jupyter notebook tutorials
- ✅ CI/CD setup (5.1)

### Sprint 4+ (Ongoing): Advanced Features
- ⏳ Reaction paths (2.4)
- ⏳ Kinetics (2.5)
- ⏳ Advanced features (4.1-4.3)
- ⏳ Reactive transport (2.6)

**Total Time to Full Standardization:** ~8-10 weeks for core features

---

## Success Metrics

### Developer Experience
- ✅ DEW code identical to HKF/Pitzer pattern
- ✅ Works in Jupyter without crashes
- ✅ Clear error messages
- ✅ < 10 lines for basic setup

### Documentation
- ✅ 10+ working examples
- ✅ Complete API reference
- ✅ Tutorial notebooks
- ✅ User guide sections

### Testing
- ✅ > 80% code coverage
- ✅ All examples run without errors
- ✅ Pass validation against experimental data
- ✅ CI/CD passes on all platforms

### Community
- ✅ Positive user feedback
- ✅ Active usage in publications
- ✅ Community contributions
- ✅ Questions answered in < 24h

---

## Risk Assessment

### High Risk
- **Database initialization crashes**: Requires careful C++/Python binding work
- **Thermodynamic inconsistencies**: Need extensive validation
- **Performance issues at extreme T/P**: May need optimization

### Medium Risk
- **API compatibility**: Breaking changes to existing code
- **Documentation lag**: Features without docs confuse users
- **Platform-specific bugs**: Windows/Linux/macOS differences

### Low Risk
- **Advanced features timeline**: Can be deferred without breaking core functionality
- **Community adoption**: Good examples will drive usage

---

## Current Status Assessment

### ✅ What Works Well
- Basic DEW thermodynamic calculations
- Water property models (multiple EOS options)
- Integration with SUPCRT minerals
- High T/P equilibrium solving

### ⚠️ What Needs Work
- API consistency with other models
- Jupyter notebook stability
- Documentation completeness
- Example coverage

### ❌ What Doesn't Work
- Direct constructors in Jupyter (`DEW2024()`)
- Reaction path with DEW
- Kinetic modeling with DEW
- Some edge cases at extreme conditions

---

## Conclusion

**Current State:** DEW is a powerful model but not fully standardized in Reaktoro.

**Required Work:** ~8-10 weeks of focused development to achieve full standardization.

**Highest Priority:** Fix database constructors and activity model API for consistent user experience.

**Expected Outcome:** DEW becomes as easy to use as HKF/Pitzer for all Reaktoro calculations.

**Recommendation:** Start with Phase 1 (Core API Standardization) immediately, as it blocks user adoption and creates frustration with crashes and inconsistent API patterns.
