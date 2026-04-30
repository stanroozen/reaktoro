# Reaktoro Three-Tier Test Architecture

## Overview

The Reaktoro test suite is organized into three distinct tiers for different purposes, each with specific characteristics, dependencies, and CI/CD usage patterns.

```
┌─────────────────────────────────────────────────────────────────────┐
│                                                                       │
│  TIER 1: CTest/C++ Tests (Deterministic, Core Unit Tests)           │
│  ├─ Catch2 C++ unit tests (135+ test files)                         │
│  ├─ Fast: ~5-10 minutes total                                       │
│  ├─ Deterministic: No external dependencies                         │
│  ├─ Essential: PR gating, pre-commit                                │
│  └─ Labels: "unit", "core"                                          │
│                                                                       │
├─────────────────────────────────────────────────────────────────────┤
│                                                                       │
│  TIER 2: Python Tests (API & Binding Tests)                         │
│  ├─ pytest-based Python binding validation                          │
│  ├─ Covers: test_dew_python.py, test_standardthermo_dew.py, etc.   │
│  ├─ Medium: ~30-60 seconds total                                    │
│  ├─ Deterministic: Requires ReactoroAPI bindings only               │
│  ├─ Purpose: Verify Python layer functionality                      │
│  └─ Labels: "python", "core"                                        │
│                                                                       │
├─────────────────────────────────────────────────────────────────────┤
│                                                                       │
│  TIER 3: External Parity/Regression (Long-running, Optional)        │
│  ├─ DEW Experimental Benchmark (regression_suite.py)                │
│  ├─ Perple_X Parity Tests (run_regression.py) - Optional            │
│  ├─ Slow: 30 minutes to 2+ hours                                    │
│  ├─ Conditional: May require external tools/data                    │
│  ├─ Purpose: Regression & long-term consistency verification        │
│  └─ Labels: "external", "regression", "nightly"                     │
│                                                                       │
└─────────────────────────────────────────────────────────────────────┘
```

## Tier Details

### Tier 1: CTest/C++ Tests

**Purpose**: Core unit and deterministic integration tests essential for PR gating.

**Scope**:
- `reaktoro-cpptests`: Compiled executable containing ~135 Catch2 test files
- Location: `Reaktoro/**/*.test.cxx`
- Compiled into: `reaktoro-cpptests` executable

**Characteristics**:
- ✅ Fast (~5-10 minutes)
- ✅ Deterministic (no external dependencies)
- ✅ Self-contained (use embedded databases)
- ✅ Always enabled

**CTest Labels**: `unit`, `core`

**Run locally**:
```bash
cd build
ctest -L "unit"
```

**CI Usage** (PR gating):
```bash
ctest -L "core"  # Runs all core tests (C++ and Python)
```

---

### Tier 2: Python Tests

**Purpose**: Verify Python API bindings and high-level functionality.

**Scope**:
- `test_dew_python.py`: DEWDatabase, WaterModels, WaterProperties
- `test_standardthermo_dew.py`: StandardThermoModelDEW binding
- `test_water_model_combinations.py`: Water model enum validation
- `test_none_options.py`: Born model "None" behavior documentation
- Location: `tests/test_*.py`
- Executed via: pytest with proper environment setup

**Characteristics**:
- ✅ Medium speed (~30-60 seconds)
- ✅ Deterministic (requires reaktoro4py module only)
- ✅ Self-contained (use embedded databases)
- ✅ Always enabled

**CTest Labels**: `python`, `core`

**Run locally**:
```bash
cd build
ctest -L "python"  # Run all Python binding tests
pytest ../tests/test_dew_python.py -v  # Direct pytest invocation
```

**CI Usage** (PR gating + full checks):
```bash
ctest -L "python"     # Binding tests only
ctest -L "core"       # All core tests (C++ + Python)
```

---

### Tier 3: External Parity/Regression Workflows

**Purpose**: Long-running benchmark and regression tests for nightly CI and release validation.

#### DEW Experimental Benchmark Suite
**Location**: `DEW_Experimental_Benchmark/regression_suite.py`

**Characteristics**:
- Centralized registry of 15+ tagged regression cases (smoke, regression, perplexdew, etc.)
- Comprehensive DEW database comparisons
- Long-running: 30 minutes to 2+ hours depending on subset

**CTest Tests**:
1. **dew-regression-smoke**: Quick subset (smoke tests)
   - Labels: `regression`, `dew`, `external`, `optional`
   - Timeout: 30 minutes
   - Purpose: Quick validation in nightly CI

2. **dew-regression-full**: Complete test suite
   - Labels: `regression`, `dew`, `external`, `nightly`
   - Timeout: 2 hours
   - Purpose: Full regression validation (nightly/release)

#### Perple_X Parity Tests (Optional, Not Yet Integrated)
**Location**: `Reaktoro/Extensions/Perple_X/run_regression.py` + `test_regression.cpp`

**Current Status**: Standalone project, not integrated into main build

**To Enable**:
1. Uncomment `add_subdirectory(${PERPLEX_EXTENSION_DIR})` in `tests/CMakeLists.txt`
2. Ensure baseline GFSM CSV files exist at `Reaktoro/Extensions/Perple_X/test/gfsm/`
3. If external Perple_X tools required, set environment variable: `export PERPLEX_MEEMUM_EXE=/path/to/meemum.exe`

**CTest Tests** (commented, waiting for integration):
```cmake
# perplex-parity-gfsm: Parity verification against committed GFSM baselines
# Labels: parity, core, regression
# Timeout: 5 minutes
```

---

## CTest Label Strategy

### Label Categories

| Label | Meaning | Purpose | Inclusion |
|-------|---------|---------|-----------|
| `core` | Essential for PR gating | CI gates, pre-commit | Always |
| `unit` | Pure unit tests | Fast feedback loop | PR checks |
| `python` | Python API tests | Binding verification | PR checks |
| `integration` | Integration tests | System validation | PR checks or nightly |
| `regression` | Baseline regression tests | Consistency verification | Nightly |
| `parity` | Cross-implementation parity | Implementation agreement | Nightly |
| `external` | Requires external tools/data | Long-running validation | Nightly or scheduled |
| `optional` | Informational, not blocking | Documentation, edge cases | Optional/nightly |
| `nightly` | Long-running tests | Scheduled runs | Nightly only |
| `dew` | DEW-specific tests | DEW feature validation | Nightly |

### Usage Patterns

#### PR Gating (Fast, Essential)
```bash
# Run only tests required for PR approval (~5-10 minutes)
ctest -L "core" -E "optional"
```

This will run:
- ✅ All C++ unit tests (Catch2)
- ✅ All Python binding tests
- ❌ DEW regression (external, not core)
- ❌ Perple_X parity (external, not core)

#### Developer Local Build
```bash
# Run all available tests locally
cd build
cmake --build . --target tests

# Or via CTest (all tests)
ctest --verbose

# Or just C++ tests
ctest -L "unit"

# Or just Python tests
ctest -L "python"
```

#### Nightly/Scheduled CI
```bash
# Run all tests including long-running regression
ctest --verbose
```

This will run:
- ✅ All C++ unit tests
- ✅ All Python binding tests
- ✅ DEW smoke regression (quick validation)
- ✅ DEW full regression (long validation)
- ✅ Perple_X parity (when available)

#### Exclude External Tests (Limited Resources)
```bash
# Skip tests that require external tools
ctest -E "external" --verbose
```

---

## CTest Execution Examples

### Show available tests
```bash
ctest --print-labels       # All labels used in project
ctest --print-test-labels # All tests with labels
```

### Run with output
```bash
ctest -V                   # Verbose (show all output)
ctest -VV                  # Extra verbose (detailed progress)
ctest --output-on-failure  # Only show failed test output
```

### Run specific test
```bash
ctest -R "dew-regression"  # Run tests matching pattern
ctest -N                   # Show tests without running
```

### Parallel execution
```bash
ctest -j 4                 # Run 4 tests in parallel
ctest -j 4 -L "core"      # 4 parallel core tests only
```

### Filter by timeout
```bash
ctest -I 0,0,10            # Run first 10 tests
ctest --stop-time HH:MM    # Stop at specific time
```

---

## Integration Points

### Local Development Workflow

1. **Quick validation** (before git commit):
   ```bash
   cd build && cmake --build . --target tests
   ```
   Runs: C++ unit tests + Python binding tests

2. **Full regression** (before opening PR):
   ```bash
   cd build && ctest -L "core"
   ```
   Runs: All PR gating tests

3. **Nightly/integration** (after merge):
   ```bash
   cd build && ctest  # All tests including external
   ```
   Runs: Everything (slow, ~2+ hours)

### CI/CD Configuration Recommendations

#### GitHub Actions / GitLab CI

**PR Check Job** (blocking):
```yaml
# Fast: ~15 minutes
ctest -L "core" -E "optional"
```

**Nightly Job** (informational):
```yaml
# Slow: ~2+ hours
ctest
```

**Release Build** (blocking + full):
```yaml
# All tests + profiling
ctest -V
ctest --cpack-config CPackConfig.cmake
```

#### Azure Pipelines

```yaml
- script: ctest -L "core" -E "optional" --verbose
  displayName: 'PR Gating Tests (Core + Python)'
  timeoutInMinutes: 15

- script: ctest --verbose
  displayName: 'Nightly Full Test Suite'
  timeoutInMinutes: 180
  condition: eq(variables['Build.Reason'], 'Schedule')
```

---

## Future Enhancements

### Proposed Tier 1 Additions
- [ ] Perple_X regression tests (requires integration of standalone extension)
- [ ] Performance benchmark regression (comparing build-to-build performance)

### Proposed Tier 2 Additions
- [ ] Additional Python API tests for new bindings
- [ ] Jupyter notebook integration tests

### Proposed Tier 3 Enhancements
- [ ] Distributed stress testing (multi-threaded scenarios)
- [ ] Long-term consistency monitoring
- [ ] Optional external Perple_X/PerplexRxn workflow validation

---

## Troubleshooting

### Python tests fail to find reaktoro module
**Symptom**: `ImportError: No module named 'reaktoro'`

**Fix**:
```bash
# Ensure PYTHONPATH is set correctly
export PYTHONPATH=/path/to/build/Reaktoro/:/path/to/build/python/package:$PYTHONPATH
cd build && ctest -L "python" -V
```

### DEW regression tests missing reference data
**Symptom**: `FileNotFoundError: test/gfsm/*.csv`

**Fix**:
Ensure committed CSV baseline files exist:
```bash
ls -la Reaktoro/Extensions/Perple_X/test/gfsm/
# Should show 50+ .csv files
```

### External Perple_X tool not found
**Symptom**: `perplex-parity tests skipped` or `Perple_X baseline data not found`

**Fix** (optional):
```bash
# Set environment before building
export PERPLEX_MEEMUM_EXE=/path/to/perplex/meemum.exe
cmake ..  # Re-run CMake
ctest -L "parity"
```

---

## File Structure Reference

```
reaktoro/
├── tests/
│   ├── CMakeLists.txt           ← Three-tier test orchestration
│   ├── pytest.ini.in            ← pytest configuration
│   ├── test_dew_python.py       ← TIER 2: DEWDatabase tests
│   ├── test_standardthermo_dew.py ← TIER 2: StandardThermoModelDEW tests
│   ├── test_water_model_combinations.py ← TIER 2: Water model validation
│   └── test_none_options.py     ← TIER 2: Born model "None" behavior
│
├── Reaktoro/
│   ├── CMakeLists.txt           ← C++ compilation (includes Catch2 tests)
│   ├── **/*.test.cxx            ← TIER 1: C++ unit tests (~135 files)
│   └── Extensions/
│       ├── Perple_X/
│       │   ├── CMakeLists.txt   ← Standalone Perple_X extension
│       │   ├── run_regression.py ← TIER 3 (optional): Parity workflow
│       │   ├── test_regression.cpp ← TIER 1 (optional): Parity executable
│       │   └── test/gfsm/*.csv  ← Baseline GFSM reference data (51+ files)
│       └── DEW/
│           └── tests/           ← DEW test data
│
└── DEW_Experimental_Benchmark/
    ├── regression_suite.py      ← TIER 3: DEW regression orchestrator
    ├── test_*.py                ← TIER 3: Individual benchmark cases
    └── regression_results/      ← TIER 3: Test output summaries
```

---

## Summary: CI/CD Decision Matrix

| Scenario | CTest Command | Time | Purpose |
|----------|---------------|------|---------|
| **PR Check** | `ctest -L "core"` | ~15 min | Gate PR approval |
| **Pre-commit** | `ctest -L "unit"` | ~5 min | Quick sanity check |
| **Developer** | `ctest` | ~15 min | Full local validation |
| **Nightly** | `ctest` | ~2 hrs | Full regression suite |
| **Release** | `ctest -V` | ~2.5 hrs | All tests + profiling |
| **Skip External** | `ctest -E "external"` | ~15 min | Limited resources |

