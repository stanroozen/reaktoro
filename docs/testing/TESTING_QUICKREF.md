# Quick Reference: Running Reaktoro Tests

## One-Line Commands

### Most Common
```bash
# PR gating (fast, ~15 minutes)
ctest -L "core"

# All tests locally (comprehensive, ~30 minutes with DEW)
ctest

# Just C++ unit tests (fast, ~5 minutes)
ctest -L "unit"

# Just Python binding tests (~1 minute)
ctest -L "python"

# DEW regression only (long, ~1-2 hours)
ctest -L "dew"
```

## Detailed Test Tiers

### Tier 1: C++ Unit Tests (Always Fast)
```bash
ctest -R "reaktoro-cpptests" -V     # Full C++ test suite with output
ctest -R "cpptests" --output-on-failure
```
- Time: ~5-10 minutes
- Coverage: 135+ test files
- Status: Essential for PR

### Tier 2: Python API Tests (Fast)
```bash
ctest -L "python" -V                # All Python binding tests
ctest -R "pytest-bindings" -V       # Just the binding tests
ctest -R "pytest-full" -V           # Full pytest suite
```
- Time: ~30 seconds (bindings) to ~2 minutes (full)
- Coverage: DEW, StandardThermo, water models
- Status: Essential for PR

### Tier 3: External Regression (Slow & Optional)
```bash
ctest -L "external" -V              # All external tests (2+ hours)
ctest -R "dew-smoke" -V             # Quick DEW validation (~20 min)
ctest -R "dew-full" -V              # Full DEW regression (~2 hours)
```
- Time: 20 minutes (smoke) to 2+ hours (full)
- Coverage: DEW benchmark cases
- Status: Nightly/scheduled only

## Filtering & Control

### Run subset of tests
```bash
ctest -R "pattern"                  # Match test name
ctest -L "label"                    # Match label
ctest -L "label" -E "pattern"       # Match label, exclude pattern
```

### Show available tests
```bash
ctest -N                            # List all tests without running
ctest --print-labels                # Show all labels used
```

### Output control
```bash
ctest -V                            # Verbose (show test output)
ctest -VV                           # Extra verbose (very detailed)
ctest --output-on-failure           # Only show failed tests
ctest -O test-output.log            # Save to file
```

### Parallel execution
```bash
ctest -j 4                          # Run 4 tests in parallel
ctest -j 4 -L "core"               # 4 parallel core tests
ctest -j                            # Auto-detect CPU count
```

### Time control
```bash
ctest --timeout 60                  # Global timeout per test (seconds)
ctest -T Test                       # Include timing info
ctest --schedule-random             # Random order (finds race conditions)
```

## Label Cheat Sheet

| Label | Tests Included | Time | Use When |
|-------|---|---|---|
| `core` | C++ + Python | ~15m | PR gating |
| `unit` | C++ unit tests only | ~5m | Pre-commit |
| `python` | Python binding tests | ~1m | Verify bindings |
| `regression` | DEW + parity tests | ~2h | Nightly validation |
| `external` | Long-running workflows | ~2h | Scheduled CI |
| `optional` | Non-blocking tests | varies | Full CI runs |

## Common Workflows

### Before committing locally
```bash
ctest -L "unit"
ctest -L "python"
```

### Before pushing to GitHub
```bash
ctest -L "core" --output-on-failure
```

### Nightly CI job
```bash
ctest --verbose
```

### Skip resource-heavy tests
```bash
ctest -E "external" -L "core"
```

### Debug failing test
```bash
ctest -R "test-name" -V --rerun-failed
ctest --rerun-failed  # Repeat last failures
```

## Environment Setup

### Python tests need reaktoro module
```bash
# In build directory, before running ctest:
cd build
cmake --build . --target reaktoro-setuptools  # Build Python bindings first
ctest -L "python"
```

### DEW tests need databases
```bash
# Automatically loaded from embedded directory
# No extra setup needed!
ctest -L "python"  # Just works
```

### Optional: Perple_X parity tests
```bash
# Set external tool path (if you have Perple_X installed)
export PERPLEX_MEEMUM_EXE=/path/to/perplex/meemum.exe

# Re-run cmake to detect tools
cd build && cmake ..

# Then enable parity tests (currently commented out in CMakeLists.txt)
ctest -L "parity"
```

## Troubleshooting Quick Fixes

| Error | Fix |
|-------|-----|
| `ImportError: No module named 'reaktoro'` | Run `cmake --build . --target reaktoro-setuptools` first |
| `test_regression executable not found` | Uncomment Perple_X integration in tests/CMakeLists.txt |
| `pytest: command not found` | Ensure Python environment has pytest installed |
| `DEW database not found` | Verify embedded databases in `Reaktoro/embedded/databases/` |
| `Test timeout exceeded` | Increase timeout: `set_tests_properties(...  TIMEOUT 3600)` |

## Reference

- **Detailed docs**: See `TESTING_ARCHITECTURE.md`
- **Build system**: `tests/CMakeLists.txt`
- **Python config**: `tests/pytest.ini.in`
- **C++ tests**: `Reaktoro/**/*.test.cxx`
- **Python tests**: `tests/test_*.py`
- **DEW regression**: `DEW_Experimental_Benchmark/regression_suite.py`

