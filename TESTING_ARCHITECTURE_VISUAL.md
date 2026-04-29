# Reaktoro Three-Tier Test Architecture - Visual Guide

## High-Level Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                                                                               │
│                  REAKTORO THREE-TIER TEST ARCHITECTURE                       │
│                                                                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                               │
│  ┌───────────────────────────────────────────────────────────────────────┐  │
│  │                                                                       │  │
│  │  TIER 1: CTest/C++ Tests                                            │  │
│  │  ├─ Compiled: Catch2 C++ unit tests (135+ files)                   │  │
│  │  ├─ Command: ctest -L "unit"                                       │  │
│  │  ├─ Time: ~5-10 minutes                                            │  │
│  │  ├─ Status: ✅ ESSENTIAL (PR gating)                               │  │
│  │  └─ Labels: "unit", "core"                                         │  │
│  │                                                                       │  │
│  │  TEST FILES: Reaktoro/**/*.test.cxx                                │  │
│  │  EXECUTABLE: reaktoro-cpptests                                     │  │
│  │  CTest Entry: reaktoro-cpptests                                    │  │
│  │                                                                       │  │
│  └───────────────────────────────────────────────────────────────────────┘  │
│           ↓                                                                   │
│  ┌───────────────────────────────────────────────────────────────────────┐  │
│  │                                                                       │  │
│  │  TIER 2: Python Tests                                              │  │
│  │  ├─ Framework: pytest via CTest                                   │  │
│  │  ├─ Command: ctest -L "python"                                    │  │
│  │  ├─ Time: ~30 seconds - 2 minutes                                 │  │
│  │  ├─ Status: ✅ ESSENTIAL (PR gating)                              │  │
│  │  └─ Labels: "python", "core"                                      │  │
│  │                                                                       │  │
│  │  TEST FILES: tests/test_*.py                                      │  │
│  │  ├─ test_dew_python.py                                            │  │
│  │  ├─ test_standardthermo_dew.py                                    │  │
│  │  ├─ test_water_model_combinations.py                              │  │
│  │  └─ test_none_options.py                                          │  │
│  │                                                                       │  │
│  │  CTest Entries:                                                    │  │
│  │  ├─ pytest-bindings-dew                                           │  │
│  │  ├─ pytest-bindings-standardthermo                                │  │
│  │  ├─ pytest-bindings-water-models                                  │  │
│  │  ├─ pytest-bindings-none-options                                  │  │
│  │  └─ pytest-full-suite (optional)                                  │  │
│  │                                                                       │  │
│  └───────────────────────────────────────────────────────────────────────┘  │
│           ↓                                                                   │
│  ┌───────────────────────────────────────────────────────────────────────┐  │
│  │                                                                       │  │
│  │  TIER 3: External Parity/Regression (OPTIONAL/NIGHTLY)            │  │
│  │  ├─ Framework: Python scripts                                      │  │
│  │  ├─ Command: ctest -L "external"                                  │  │
│  │  ├─ Time: 20 minutes - 2+ hours                                   │  │
│  │  ├─ Status: ⏳ OPTIONAL (nightly/scheduled)                        │  │
│  │  └─ Labels: "external", "regression", "nightly"                   │  │
│  │                                                                       │  │
│  │  ┌─────────────────────────────────────────────────────────────┐  │  │
│  │  │ DEW Experimental Benchmark Suite                           │  │  │
│  │  ├─ Script: DEW_Experimental_Benchmark/regression_suite.py    │  │  │
│  │  ├─ Cases: 15+ tagged regression cases                        │  │  │
│  │  ├─ CTest Entries:                                            │  │  │
│  │  │  ├─ dew-regression-smoke (~20 min)                         │  │  │
│  │  │  └─ dew-regression-full (~2 hours)                         │  │  │
│  │  ├─ Status: ✅ ENABLED                                        │  │  │
│  │  └─ Conditional: Only if regression_suite.py exists           │  │  │
│  │  └─────────────────────────────────────────────────────────────┘  │  │
│  │                                                                       │  │
│  │  ┌─────────────────────────────────────────────────────────────┐  │  │
│  │  │ Perple_X Parity Tests (NOT YET INTEGRATED)                 │  │  │
│  │  ├─ Script: Reaktoro/Extensions/Perple_X/run_regression.py    │  │  │
│  │  ├─ Executable: test_regression.cpp                           │  │  │
│  │  ├─ Baselines: 51+ GFSM CSV files ✅ Committed                │  │  │
│  │  ├─ Status: Commented out (optional integration)              │  │  │
│  │  ├─ To Enable: Uncomment add_subdirectory() in CMakeLists.txt │  │  │
│  │  └─ CTest Entry (when enabled): perplex-parity-gfsm           │  │  │
│  │  └─────────────────────────────────────────────────────────────┘  │  │
│  │                                                                       │  │
│  └───────────────────────────────────────────────────────────────────────┘  │
│                                                                               │
└─────────────────────────────────────────────────────────────────────────────┘
```

## CTest Execution Flow

```
ctest
  │
  ├─→ TIER 1: C++ Unit Tests (5-10 min)
  │   ├─ reaktoro-cpptests [CORE, UNIT] ✅
  │   └─ [Optional] perplex-parity-gfsm [PARITY, CORE] ⏳
  │
  ├─→ TIER 2: Python Tests (30 sec - 2 min)
  │   ├─ pytest-bindings-dew [CORE, PYTHON] ✅
  │   ├─ pytest-bindings-standardthermo [CORE, PYTHON] ✅
  │   ├─ pytest-bindings-water-models [CORE, PYTHON] ✅
  │   ├─ pytest-bindings-none-options [CORE, PYTHON] ✅
  │   └─ pytest-full-suite [OPTIONAL, PYTHON] ⏳
  │
  └─→ TIER 3: External Tests (20 min - 2+ hours)
      ├─ dew-regression-smoke [EXTERNAL, REGRESSION] ✅
      └─ dew-regression-full [EXTERNAL, NIGHTLY] ✅
```

## Label Filtering Matrix

```
Command                              What Runs                    Time
───────────────────────────────────────────────────────────────────────
ctest -L "core"                      Tier 1 + 2 (no external)    ~15 min
ctest -L "core" -E "optional"        Essential tests only         ~12 min
ctest -L "unit"                      C++ tests only              ~5 min
ctest -L "python"                    Python tests only           ~1 min
ctest -L "external"                  External + parity tests     ~2 hours
ctest -L "dew"                       DEW regression only         20 min - 2h
ctest -L "regression"                All regression tests        ~2 hours
ctest                                ALL tests                   ~2 hours
ctest -E "external"                  Skip external tests         ~15 min
```

## Development Workflow

```
┌──────────────────────────────────────────────────────────────────┐
│                                                                  │
│  LOCAL DEVELOPMENT                                              │
│  ┌────────────────────────────────────────────────────────────┐ │
│  │ 1. Make code changes                                       │ │
│  │ 2. Build: cmake --build . --target reaktoro-setuptools    │ │
│  │ 3. Test: ctest -L "core"  ← PR gating tests              │ │
│  │          (C++ unit + Python binding)                      │ │
│  │ 4. Commit & Push                                          │ │
│  └────────────────────────────────────────────────────────────┘ │
│                          ↓                                       │
│  GITHUB CI/CD (PR Check)                                        │
│  ┌────────────────────────────────────────────────────────────┐ │
│  │ 1. Build Reaktoro (10 min)                               │ │
│  │ 2. Run: ctest -L "core" -E "optional"                   │ │
│  │    Result: PR approved or blocked                        │ │
│  └────────────────────────────────────────────────────────────┘ │
│                          ↓                                       │
│  NIGHTLY CI/CD (Full Suite)                                     │
│  ┌────────────────────────────────────────────────────────────┐ │
│  │ 1. Build Reaktoro (10 min)                               │ │
│  │ 2. Run: ctest --verbose (all tests)                      │ │
│  │    ├─ Tier 1: C++ unit tests (5-10 min)                │ │
│  │    ├─ Tier 2: Python tests (1-2 min)                   │ │
│  │    └─ Tier 3: External regression (1-2 hours)          │ │
│  │ 3. Generate report: Pass/fail/timeout stats              │ │
│  └────────────────────────────────────────────────────────────┘ │
│                                                                  │
└──────────────────────────────────────────────────────────────────┘
```

## Test Dependencies and Execution Order

```
Build Phase (cmake --build .)
│
├─→ [ALWAYS] Compile Reaktoro library
│   └─→ Output: libReaktoro.a / Reaktoro.dll
│
├─→ [ALWAYS] Compile C++ tests
│   ├─→ Depends on: Reaktoro library + Catch2
│   └─→ Output: reaktoro-cpptests executable
│
├─→ [ALWAYS] Build Python bindings
│   ├─→ Depends on: Reaktoro library
│   └─→ Output: reaktoro4py Python module
│
└─→ [OPTIONAL] Build Perple_X extension
    ├─→ Depends on: Reaktoro library (if integrated)
    └─→ Output: test_regression executable


Test Phase (ctest -L ...)
│
├─→ [ALWAYS] Tier 1: C++ Unit Tests (reaktoro-cpptests)
│   ├─→ Fast, deterministic
│   └─→ Can run in parallel
│
├─→ [ALWAYS] Tier 2: Python Binding Tests (pytest-*)
│   ├─→ Requires: reaktoro4py module
│   ├─→ Fast, deterministic
│   └─→ Can run in parallel
│
├─→ [OPTIONAL] Tier 1.5: Perple_X Parity (perplex-parity-gfsm)
│   ├─→ Fast, deterministic
│   ├─→ Requires: test_regression executable
│   └─→ Depends on: Committed CSV baseline data
│
└─→ [NIGHTLY] Tier 3: External Regression (dew-regression-*)
    ├─→ Slow, semi-deterministic (may vary with seed/platform)
    └─→ Should run sequentially (resource intensive)
```

## File Organization

```
reaktoro/
├── tests/
│   ├── CMakeLists.txt          ← MODIFIED: Three-tier orchestration
│   ├── pytest.ini.in           ← pytest configuration
│   ├── test_dew_python.py      ← TIER 2: DEWDatabase tests
│   ├── test_standardthermo_dew.py ← TIER 2: StandardThermo tests
│   ├── test_water_model_combinations.py ← TIER 2: Water models
│   └── test_none_options.py    ← TIER 2: Born model tests
│
├── Reaktoro/
│   ├── CMakeLists.txt          ← C++ compilation
│   ├── **/*.test.cxx           ← TIER 1: 135+ C++ unit tests
│   └── Extensions/
│       ├── Perple_X/
│       │   ├── CMakeLists.txt  ← OPTIONAL: Standalone extension
│       │   ├── run_regression.py ← TIER 3 script (optional)
│       │   ├── test_regression.cpp ← TIER 1 executable (optional)
│       │   └── test/gfsm/*.csv ← BASELINE DATA (51+ files) ✅
│       └── DEW/
│           └── tests/          ← DEW test data
│
├── DEW_Experimental_Benchmark/
│   ├── regression_suite.py     ← TIER 3: Main orchestrator
│   ├── test_*.py               ← TIER 3: Individual cases
│   └── regression_results/     ← Output directory
│
├── TESTING_ARCHITECTURE.md     ← NEW: Comprehensive guide
├── TESTING_QUICKREF.md         ← NEW: Quick reference
└── TEST_IMPLEMENTATION_CHECKLIST.md ← NEW: Implementation status
```

## Legend

- ✅ Implemented & working
- ⏳ Optional or not yet integrated
- 🔄 Needs configuration
- ❌ Not available

## Related Documents

1. **TESTING_ARCHITECTURE.md** - Comprehensive guide with CI/CD recommendations
2. **TESTING_QUICKREF.md** - Developer quick reference for CTest commands
3. **TEST_IMPLEMENTATION_CHECKLIST.md** - Implementation status and verification steps

---

**Quick Start**:
```bash
# PR check (local before pushing)
ctest -L "core"

# Developer validation
ctest -L "unit"

# Nightly (all tests)
ctest
```

