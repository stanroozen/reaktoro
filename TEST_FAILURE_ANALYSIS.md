# Test Failure Analysis & Fixes

## Summary of Issues

| Issue | Root Cause | Affected Tests | Fix |
|-------|-----------|---------------|----|
| **Missing Python Bindings** | `.pyd` not recompiled after C++ additions | `test_dew_ph_constraint.py`, `test_dew_sensitivity.py`, `test_kinetics_dew.py`, `test_transport_dew.py` | Rebuild C++ or update test code to work with available bindings |
| **Species Naming Mismatch** | DEW DB uses `WATER,AQ`/`H+` vs test expects `H2O(aq)`/`H+(aq)` | All tests using DEW DB | Update tests to use correct species names from database |
| **Hardcoded Empty Paths** | Tests reference `build-msvc\Reaktoro\Release` (empty) | Multiple smoke tests | Use dynamic path discovery (already fixed for `test_gfsm_handoff.py`) |

---

## Issue 1: Missing Python Bindings

### Problem
Tests call:
```python
params = rkt.ActivityModelParamsDEW()
params.dhModel = rkt.ActivityDHModel.ExtendedDH
```

But current `.pyd` doesn't have:
- `ActivityModelParamsDEW` class
- `ActivityDHModel` enum

### Root Cause
The bindings exist in source (`ActivityModelDEW.py.cxx`) but `.pyd` wasn't recompiled after those changes.

### Source Code Status
✅ **Already in source** — `Reaktoro/Models/ActivityModels/ActivityModelDEW.py.cxx`:
```cpp
py::enum_<ActivityDHModel>(m, "ActivityDHModel")
    .value("ExtendedDH", ActivityDHModel::ExtendedDH)
    .value("Davies",     ActivityDHModel::Davies)
    .export_values();

py::class_<ActivityModelParamsDEW>(m, "ActivityModelParamsDEW")
    .def(py::init<>())
    .def_readwrite("dhModel", &ActivityModelParamsDEW::dhModel)
    ...
```

### Solution: Option A (Recommended)
Update tests to create params without the intermediate step (factory pattern):

```python
# Instead of:
params = rkt.ActivityModelParamsDEW()
params.dhModel = rkt.ActivityDHModel.ExtendedDH
model = rkt.ActivityModelDEW(params)

# Use:
model = rkt.ActivityModelDEW()  # No-arg constructor uses defaults
```

### Solution: Option B
Rebuild the C++ to include these new bindings.

---

## Issue 2: Species Naming Mismatch

### Problem
Tests request species like:
```python
aq = rkt.AqueousPhase("H2O(aq) H+(aq) OH-(aq) Mg+2(aq) MgOH+(aq)")
```

But DEW database has:
- `WATER,AQ` (not `H2O(aq)`)
- `H+` (not `H+(aq)`)
- `OH-` (not `OH-(aq)`)
- `Mg+2` (not `Mg+2(aq)`)
- `MgOH+` (not `MgOH+(aq)`)

### Root Cause
DEW database species naming conventions don't match test expectations. The `(aq)` suffix is not used consistently.

### Solution
Update test species lists to match database naming:

```python
# Before:
aq = rkt.AqueousPhase("H2O(aq) H+(aq) OH-(aq) Mg+2(aq) MgOH+(aq)")

# After:
aq = rkt.AqueousPhase("WATER,AQ H+ OH- Mg+2 MgOH+")
```

### Verified Species Map
```
H2O(aq)     → WATER,AQ
H+(aq)      → H+
OH-(aq)     → OH-
Mg+2(aq)    → Mg+2
MgOH+(aq)   → MgOH+
SiO2(aq)    → SiO2_aq
H2(aq)      → H2_aq
O2(aq)      → O2_aq
```

---

## Issue 3: Hardcoded Paths

### Problem
Tests hardcode `build-msvc\Reaktoro\Release`:
```python
sys.path.insert(0, r"c:\Users\stanroozen\...\build-msvc\Reaktoro\Release")
```

This directory is empty; the `.pyd` is at `temp_build\build-dew\Reaktoro\Release`.

### Solution
Replace with dynamic discovery (already done for `test_gfsm_handoff.py`):

```python
from pathlib import Path

def _setup_path():
    testing_root = Path(__file__).parent.parent.parent
    repo_root = testing_root.parent if testing_root.name == "Testing" else testing_root

    search_dirs = [
        repo_root / "temp_build" / "build-dew" / "Reaktoro" / "Release",
        repo_root / "build-msvc" / "Reaktoro" / "Release",
        repo_root / "build" / "Reaktoro" / "Release",
    ]
    for d in search_dirs:
        if d.exists() and str(d) not in sys.path:
            sys.path.insert(0, str(d))
            break

_setup_path()
```

---

## Recommended Fix Strategy

### Phase 1: Immediate (No Recompile Needed)
1. Update all smoke tests to use correct species names
2. Replace hardcoded paths with dynamic discovery
3. Simplify activity model creation (no params needed for defaults)

### Phase 2: After C++ Recompile
1. Re-enable tests that need `ActivityModelParamsDEW` + `dhModel` binding
2. All 30 smoke tests should pass

---

## Test Files to Fix

| File | Issues | Lines |
|------|--------|-------|
| `test_dew_ph_constraint.py` | Hardcoded path, species names, params binding | 8-11, 30-35, 58-65 |
| `test_dew_sensitivity.py` | Hardcoded path, species names | 14-17, 32-37 |
| `test_build_system.py` | Species names (via imported builder) | ~516 in quartz_solubility_analysis_v2_dew24.py |
| `test_kinetics_dew.py` | Species names | ~40 in _make_system |
| `test_transport_dew.py` | Species names | ~43 in _make_system |
| `test_h2o_aq_dew.py` | Species names | Unknown (need to check) |
| `test_h2o_dummy.py` | Species names | Unknown (need to check) |

---

## Notes

- ✅ `test_gfsm_handoff.py` already has dynamic path discovery
- ✅ `test_dew_water_models.py` passes (doesn't use problematic species/bindings)
- ✅ Conftest.py has auto-discovery function
- ⚠️ Some species are used indirectly through imported builder modules
