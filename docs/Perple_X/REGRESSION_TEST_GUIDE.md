# Regression Testing Guide for Perple_X Implementation

## Overview

This guide describes how to perform comprehensive regression testing of the Reaktoro Perple_X implementation against the original Perple_X code.

## One-Command Workflow

Use the unified runner:

```bash
python run_regression.py --meemum-exe "C:\\Program Files (x86)\\Perplex\\meemum.exe"
```

Strict GFSM gate (recommended for release parity):

```bash
python run_regression.py \
  --meemum-exe "C:\\Program Files (x86)\\Perplex\\meemum.exe" \
  --require-full-gfsm-pure
```

Optional build integration:

```bash
python run_regression.py \
  --meemum-exe "C:\\Program Files (x86)\\Perplex\\meemum.exe" \
  --build-cmd "cmake --build build --config Release --target test_regression"
```

What it does:

1. Regenerates GFSM reference CSV tables from Perple_X (`generate_gfsm_meemum_references.py`).
2. Optionally builds `test_regression`.
3. Runs `test_regression` to compare Reaktoro against Perple_X references.
4. With `--require-full-gfsm-pure`, enforces presence of all 51 GFSM reference CSVs in
  `test/gfsm/`. Coverage:
  - **Pure H2O** (7): `h2o_mrk`, `h2o_hsmrk`, `h2o_cork`, `h2o_pseos`, `h2o_haar`, `h2o_zd05`, `h2o_zd09`
  - **Pure CO2** (5): `co2_mrk`, `co2_hsmrk`, `co2_cork`, `co2_pseos`, `co2_zd09`
    _(note: `co2_brmrk` / Bottinga & Richet 1981 excluded — EOS diverges for pure CO2)_
  - **Pure CH4** (3): `ch4_mrk`, `ch4_hsmrk`, `ch4_zd09`
  - **Binary H2O+CO2** (11): `h2o_co2_{H2O-eos}_mrk` (7 cases varying H2O) + `h2o_co2_mrk_{CO2-eos}` (4 cases varying CO2)
  - **Binary H2O+CH4** (9): `h2o_ch4_{H2O-eos}_mrk` (7) + `h2o_ch4_mrk_{CH4-eos}` (2)
  - **Binary CO2+CH4** (7): `co2_ch4_{CO2-eos}_mrk` (5) + `co2_ch4_mrk_{CH4-eos}` (2)
  - **Ternary** (3): `h2o_co2_ch4_mrk`, `h2o_co2_ch4_hsmrk`, `h2o_co2_ch4_zd09`
  - **O2/Redox Stress** (6): `redox_o2_excess_mrk`, `redox_o2_excess_hot`, `redox_h2_excess_hot`,
    `redox_co_bias_hot`, `redox_o2_excess_hybrid`, `redox_mixed_hybrid`
  - These are driven by explicit matrix rows in `test/gfsm_case_matrix.csv`.

## Test Strategy

### 1. Reference Data Generation
Generate "ground truth" GFSM results from Perple_X meemum using the type-39 solution model workflow.

### 2. Reaktoro Computation
Run identical calculations using Reaktoro Perple_X implementation.

### 3. Comparison & Validation
Compare results with strict tolerances and identify any deviations.

### Scope Lock (GFSM Only)

This regression suite is intentionally scoped to:

1. GFSM-relevant pure EoS checks (pure H2O/CO2 anchor states used by GFSM pathways)
2. GFSM mixture checks (H2O-CO2 binary and P-T/composition series)
3. DEW-speciation energy parity via HKF Gibbs comparison (`hkf_reference.csv`)

Component-only diagnostics (standalone g-function, Debye-Huckel, Born omega) are kept in
code but not part of the default pass/fail run list.

Additionally, the unified runner is pinned to GFSM/meemum generation and can enforce
full GFSM coverage with `--require-full-gfsm-pure`.

### Architecture Snapshot

The current regression architecture is matrix-driven and generator-backed:

1. Reference generation layer
  - `generate_gfsm_meemum_references.py` regenerates GFSM CSVs from Perple_X meemum.
2. Coverage definition layer
  - `test/gfsm_regression_matrix.csv` defines enabled regression rows.
  - `test/gfsm_case_matrix.csv` adds explicit stress/redox rows.
3. Execution and parity layer
  - `test_regression.cpp` runs EOS gate checks, matrix volume parity, DEW/HKF parity,
    and integration-flow assertions.
4. Orchestration layer
  - `run_regression.py` wires generation and execution, including strict gate mode
    with `--require-full-gfsm-pure`.

Current strict gate footprint is 51 GFSM reference cases:
- 15 pure anchors (H2O/CO2/CH4)
- 27 binary mixture cases
- 3 ternary H2O-CO2-CH4 cases
- 6 O2/redox stress cases

---

## Files

| File | Purpose |
|------|---------|
| `test_regression.cpp` | Main C++ regression test suite |
| `run_regression.py` | Unified GFSM/meemum regression workflow runner |
| `reference_data.json` | Reference data output (generated) |
| `regression_test_report.txt` | Test results report (generated) |
| `test/h2o_co2_mrk.csv` | MRK reference CSV (ifug=0) |
| `test/h2o_co2_hybrid_mrk.csv` | Legacy hybrid-ifug=2 reference CSV (kept for archival comparison only; not part of GFSM regression gate) |
| `test/h2o_co2_hsmrk.csv` | HSMRK H2O-CO2 reference CSV (ifug=1) |
| `test/h2o_co2_cork.csv` | CORK H2O-CO2 reference CSV (ifug=5) |
| `test/gfsm/h2o_*.csv` | Pure H2O meemum-GFSM references (7: mrk/hsmrk/cork/pseos/haar/zd05/zd09) |
| `test/gfsm/co2_*.csv` | Pure CO2 meemum-GFSM references (5: mrk/hsmrk/cork/pseos/zd09; brmrk excluded) |
| `test/gfsm/ch4_*.csv` | Pure CH4 meemum-GFSM references (3: mrk/hsmrk/zd09) |
| `test/gfsm/h2o_co2_*.csv` | Binary H2O+CO2 mixture references (11: systematic EOS variation) |
| `test/gfsm/h2o_ch4_*.csv` | Binary H2O+CH4 mixture references (9: systematic EOS variation) |
| `test/gfsm/co2_ch4_*.csv` | Binary CO2+CH4 mixture references (7: systematic EOS variation) |
| `test/gfsm/h2o_co2_ch4_*.csv` | Ternary H2O+CO2+CH4 mixture references (3: mrk/hsmrk/zd09) |
| `test/gfsm/redox_*.csv` | O2-bearing/redox-stress OHC references (6 explicit compositions) |
| `test/gfsm_case_matrix.csv` | Matrix for extra meemum GFSM cases (O2/redox stress states) |
| `generate_gfsm_meemum_references.py` | meemum-driven GFSM reference generator |

---

## Step 1: Generate Reference Data from Perple_X

### Requirements
- Perple_X installation (version 7.1.6+)
- Compiled `meemum` executable
- Python 3.7+

### Running Reference Data Generation

```bash
# Generate GFSM reference CSVs via meemum
python generate_gfsm_meemum_references.py \
  --meemum-exe "C:\Program Files (x86)\Perplex\meemum.exe" \
  --template-dat "C:\Users\stanroozen\Documents\Projects\Perplex\Perple_X\test\weigang\gfsm_fluid_probe.dat"

# Navigate to test directory
cd Reaktoro/Extensions/Perple_X

# Or invoke the unified GFSM regression workflow
python run_regression.py \
  --meemum-exe "C:\Program Files (x86)\Perplex\meemum.exe" \
  --require-full-gfsm-pure
```

### Reference Data Format

```json
[
  {
    "test_case": "pure_h2o_hsmrk",
    "state": {
      "T": 500.0,
      "P": 1000.0,
      "composition": [1.0],
      "species": [1],
      "ln_fugacity": [-0.123456],
      "volumes": [18.234],
      "total_volume": 18.234,
      "epsilon": 20.45,
      "g_function": 0.0123,
      "adh_factor": -0.456
    }
  }
]
```

---

## Step 2: Update Test Suite with Reference Data

The C++ tests load reference values from the CSV files in `test/`:

```cpp
// CSV lookup example
const auto table = loadCsvTable("test/h2o_co2_mrk.csv");
const auto row = findRow(table, 1000.0, 523.15, 0.5);
const double fH2O = getValue(table, *row, "f(H2O)");
ref.ln_f = {std::log(fH2O)};
```

---

## Step 3: Compile and Run Regression Tests

### Compilation

```bash
# Using CMake
mkdir build && cd build
cmake ..
make test_regression

# Or directly with g++
g++ -std=c++17 -O2 -I../../.. \
    test_regression.cpp \
    PerpleXFluidModel.cpp \
    PerpleXHybridEos.cpp \
    PerpleXMrkMixture.cpp \
    PerpleXPureEos.cpp \
    PerpleXElectrolyte.cpp \
    PerpleXHKF.cpp \
    -o test_regression
```

### Running Tests

```bash
./test_regression

# Output will show:
# ======================================================================
# Perple_X Regression Tests - Reaktoro Implementation
# ======================================================================
#
# Running tests...
#
# Test: Pure H2O HSMRK
# Status: ✅ PASSED
#
# Test: H2O-CO2 Binary (50-50)
# Status: ❌ FAILED
# Errors:
#   - CO2 fugacity coefficient: computed=1.234e-01, reference=1.235e-01, rel_error=8.1e-05 (tol=1.0e-06)
# ...
```

---

## Test Coverage

### Active Test Set (Current Runner)

The default GFSM run executes 11 test groups:

1. GFSM EOS coverage gate
2. Pure H2O MRK (ifug=2 defaults)
3. Pure H2O HSMRK (ifug=1)
4. Pure CO2 MRK (ifug=2 defaults)
5. Pure CO2 CORK (ifug=5)
6. H2O-CO2 binary (50-50)
7. H2O-CO2 composition series
8. Dielectric H2O-CO2 parity
9. DEW speciation energy parity (HKF Gibbs)
10. GFSM full matrix volume parity
11. Integration-flow scenarios

### Conditions and Composition Envelope

- Pressure and temperature ranges are matrix-defined and include hot/redox stress states.
- Composition classes include pure, binary, ternary, and O2/redox stress mixtures.
- Species basis follows the GFSM allowed set used by current generator and solver paths.

---

## Tolerance Levels

| Property | Tolerance | Type | Rationale |
|----------|-----------|------|-----------|
| Fugacity coefficient | 1×10⁻⁶ | Relative | Numerical precision |
| Volume | 1×10⁻⁶ | Relative | Numerical precision |
| Dielectric constant | 1×10⁻⁴ | Absolute | Physical measurement |
| g-function | 1×10⁻⁸ | Absolute | Å precision |
| DH factor | 1×10⁻⁶ | Relative | Activity coefficient |
| HKF Gibbs | 100 | Absolute | J/mol precision |
| ΔG mix proxy | 25 | Absolute | J/mol precision |

### Adjusting Tolerances

If systematic deviations are found, tolerances can be adjusted in `test_regression.cpp`:

```cpp
constexpr double TOL_FUGACITY = 1e-6;     // Adjust if needed
constexpr double TOL_VOLUME = 1e-6;
constexpr double TOL_DIELECTRIC = 1e-4;
// etc.
```

Volume parity additionally uses dynamic tolerance tiers in `test_regression.cpp`
for boundary and extreme-state robustness.

---

## Interpreting Results

### ✅ PASSED
All computed values within tolerance. Implementation verified.

### ❌ FAILED
One or more values exceed tolerance. Investigate:

1. **Systematic error**: All values off by same amount → Check constants
2. **Specific condition error**: Only certain P-T-X fail → Check boundary conditions
3. **Component error**: Only one property fails → Check that component implementation

### Common Issues

| Error Pattern | Likely Cause | Fix |
|---------------|--------------|-----|
| All fugacities off by ~1e-5 | MRK mixing rule constant | Check `a_mix`, `b_mix` formulas |
| Volumes wrong at high P | Compressibility term | Check pressure dependence |
| Dielectric wrong | Species parameters | Verify H&L/H&M coefficients |
| g-function wrong at ρ<0.5 | Region I coefficients | Check polynomial coefficients |
| HKF wrong at high T | Temperature terms | Check `ft = T - theta` |

---

## CI Metrics and Thresholds

### Published Baseline (Current Report)

From `regression_test_report.txt`:

- Total tests: 11
- Passed: 11
- Failed: 0
- Assertions executed: 116
- Assertions passed: 116
- Assertions failed: 0
- Row-level assertions executed: 51
- Row-level assertions passed: 51
- Row-level assertions failed: 0

### CI Gate Policy

Recommended CI thresholds:

1. `Failed == 0`
2. `Assertions Failed == 0`
3. `Row-level Assertions Failed == 0`
4. Strict GFSM coverage mode enabled (`--require-full-gfsm-pure`)
5. Tolerance constants unchanged unless explicitly reviewed

### Strict CI Command

```bash
python run_regression.py \
  --meemum-exe "$PERPLEX_MEEMUM" \
  --require-full-gfsm-pure
```

---

## Automated Testing

### Continuous Integration

Add to CI pipeline:

```yaml
# .github/workflows/test.yml
- name: Run Perple_X Regression Tests
  run: |
    cd build
    ./test_regression
    if [ $? -ne 0 ]; then
      cat regression_test_report.txt
      exit 1
    fi
```

### Nightly Tests

Run comprehensive GFSM regression nightly:

```bash
#!/bin/bash
# nightly_regression.sh

# Generate GFSM references and run tests
python run_regression.py \
  --meemum-exe "$PERPLEX_MEEMUM" \
  --require-full-gfsm-pure

# Email report if failures
if [ $? -ne 0 ]; then
    mail -s "Perple_X Regression Failures" dev@example.com < regression_test_report.txt
fi
```

---

## Advanced Testing

### Benchmarking

Compare performance against Perple_X:

```cpp
#include <chrono>

auto start = std::chrono::high_resolution_clock::now();

// Run 1000 iterations
for (int i = 0; i < 1000; ++i) {
    auto state = model.compute({1, 2}, y, P, T, opts);
}

auto end = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

std::cout << "Average time: " << duration.count() / 1000.0 << " μs\n";
```

### Fuzzing

Test edge cases with random inputs:

```cpp
#include <random>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> T_dist(273.15, 1273.15);
std::uniform_real_distribution<> P_dist(1.0, 10000.0);
std::uniform_real_distribution<> X_dist(0.0, 1.0);

for (int i = 0; i < 10000; ++i) {
    double T = T_dist(gen);
    double P = P_dist(gen);
    double x = X_dist(gen);

    std::array<double, 19> y{};
    y[0] = 1.0 - x;
    y[1] = x;

    try {
        auto state = model.compute({1, 2}, y, P, T, opts);
        // Check for NaN/Inf
        assert(!std::isnan(state.vol));
        assert(!std::isinf(state.vol));
    } catch (...) {
        // Log conditions that cause failures
    }
}
```

---

## Reporting Issues

If regression tests fail consistently:

1. **Document the failure**:
   - P-T-X conditions
   - Expected vs computed values
   - Error magnitude
   - Test name

2. **Check Perple_X version**:
   ```bash
   ./fluids --version
   ```

3. **Verify reference data**:
   - Re-run Perple_X manually
   - Check for Perple_X updates
   - Compare with published results

4. **Create minimal reproducible example**:
   ```cpp
   // Minimal failing case
   PerpleXFluidModel model;
   auto state = model.compute({1}, {1.0}, 1000.0, 500.0, opts);
   std::cout << "Computed: " << state.g[0] << "\n";
   std::cout << "Expected: " << -0.123456 << "\n";
   ```

5. **Submit issue** with:
   - Test output
   - `reference_data.json` excerpt
   - Perple_X version
   - Compiler/OS info

---

## Validation Checklist

Before declaring implementation complete:

- [ ] All 11 active GFSM test groups pass
- [ ] Assertions failed = 0
- [ ] Row-level assertions failed = 0
- [ ] Strict coverage gate (`--require-full-gfsm-pure`) passes
- [ ] No systematic errors > 0.1% in reviewed outputs
- [ ] Rare failures documented and understood
- [ ] Performance within 2× of Perple_X
- [ ] No memory leaks (run with valgrind)
- [ ] Thread-safe (run with thread sanitizer)
- [ ] Documentation complete

---

## Summary

**Regression testing ensures**:
✅ Numerical accuracy to machine precision
✅ Correct implementation of all formulas
✅ Parity with Perple_X across full P-T-X range
✅ Confidence for production use

**Test early, test often!**
