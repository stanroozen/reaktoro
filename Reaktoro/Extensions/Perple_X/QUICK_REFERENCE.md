# Perple_X GFSM Quick Reference Card

## Three Ways to Use Perple_X GFSM in Reaktoro

### 1️⃣ GFSM Mixture Model (Recommended for most users)
```cpp
#include <Reaktoro/Extensions/Perple_X.hpp>

// Create species list (must be exactly the 13 GFSM species)
SpeciesList species;
species.push_back({.name = "H2O", .formula = Formula("H2O"), ...});
species.push_back({.name = "CO2", .formula = Formula("CO2"), ...});
// ... add 11 more species: CH4, H2, CO, H2S, SO2, O2, N2, NH3, HF, C2H6, HCl

// Create GFSM model
auto model = Reaktoro::ActivityModelPerplexGFSM();  // ← Uses MRK baseline

// Create phase
auto phase = Phase(species).setActivityModel(model);

// Evaluate at T, P, composition
ActivityProps props = phase.evaluate(700.0, 1.0e8, x);
// ✅ Gets: ln_g (fugacity coefficients), Vx (molar volume)
```

**Best for**: Fluid mixture thermodynamics, EOS calculations, phase equilibria

---

### 2️⃣ Reference State Properties (For coupling with activity model)
```cpp
#include <Reaktoro/Extensions/Perple_X.hpp>

// Configure reference state for H2O
StandardThermoModelParamsPerplexGFSM h2o_params;
h2o_params.G0 = -235000.0;     // J/mol
h2o_params.H0 = -286000.0;     // J/mol
h2o_params.V0 = 1.89e-5;       // m³/mol
h2o_params.Tmax = 1200.0;      // K

// Create standard thermo model
auto stdthermo = StandardThermoModelPerplexGFSM(h2o_params);

// Evaluate at T, P
StandardThermoProps props = stdthermo(700.0, 1.0e8);
// ✅ Gets: G0, H0, V0 (reference state)
```

**Best for**: Database integration, coupling activity models with reference states

---

### 3️⃣ Individual Pure EOS (For validation/research)
```cpp
#include <Reaktoro/Extensions/Perple_X.hpp>

using namespace Reaktoro::PerpleX;

// Evaluate H2O with HSMRK at P=1000 bar, T=700 K
double vol = 50.0;              // cm³/mol (initial guess)
double P_bar = 1000.0;
double T_K = 700.0;

double ln_fugacity = hsmrkf(vol, 1, P_bar, T_K);  // ← vol updated to solution
// ✅ Gets: ln_fugacity (pure H2O), vol (updated to converged value)
```

**Available for H2O (7 options)**:
- `hsmrkf()` - HSMRK (Kerrick & Jacobs)
- `crkH2O()` - CORK (Holland & Powell)
- `pseos()` - PSEOS
- `brmrk()` - BRMRK
- `haar()` - Haar empirical
- `zhdh2o()` - Zhang-Duan 2005
- `zd09pr()` - Zhang-Duan 2009

**Available for CO2 (6 options)**: `hsmrkf()`, `crkCO2()`, `pseos()`, `brmrk()`, `zd09pr()`, MRK

**Available for CH4 (3 options)**: `hsmrkf()`, `zd09pr()`, MRK

**Available for others (9 species)**: MRK only (CO, H2, H2S, O2, SO2, N2, NH3, HF, C2H6, HCl)

**Best for**: Direct EOS evaluation, validation, research, low-level calculations

---

## Which Interface to Use?

| Use Case | Interface | Function |
|----------|-----------|----------|
| Fluid thermodynamics at (T,P,x) | 1️⃣ Mixture | `ActivityModelPerplexGFSM()` |
| Get mixture molar volume | 1️⃣ Mixture | `ActivityModelPerplexGFSM()` |
| Get fugacity coefficients | 1️⃣ Mixture | `ActivityModelPerplexGFSM()` |
| Setup reference state | 2️⃣ RefState | `StandardThermoModelPerplexGFSM()` |
| Pure species EOS at (T,P) | 3️⃣ PureEOS | `hsmrkf()`, `crkH2O()`, etc. |
| Compare EOS options | 3️⃣ PureEOS | Try multiple functions |
| Validate GFSM results | 3️⃣ PureEOS | Cross-check with pure EOS |

---

## 13 GFSM Species (Model Type 39)

```
H2O, CO2, CH4, H2, CO, H2S, SO2, O2, N2, NH3, HF, C2H6, HCl
```

**MUST use exactly these 13 species** - no subsets, no additions

---

## Unit Conventions

| Quantity | Reaktoro API | Perple_X Internal |
|----------|--------------|-------------------|
| Pressure | **Pa** | bar (converted automatically) |
| Temperature | **K** | K (same) |
| Volume | **m³/mol** | cm³/mol (converted automatically) |
| Energy | **J/mol** | J/mol (same) |

**Good News**: Unit conversions handled automatically in wrapper functions!

---

## Hybrid EOS Selection (Advanced)

```cpp
// Configure which pure EOS to use for each species
ActivityModelParamsPerplexGFSM params;
params.hybridEosOptions.waterEos = HybridEosOptions::WaterEos::CORK;
params.hybridEosOptions.co2Eos = HybridEosOptions::CO2Eos::BRMRK;
params.hybridEosOptions.ch4Eos = HybridEosOptions::CH4Eos::HSMRK;

// All other species stay on MRK (cannot be changed)

auto model = ActivityModelPerplexGFSM(params);
```

**Default**: All species use MRK (Modified Redlich-Kwong)

---

## Common Workflows

### Workflow A: Simple Mixture Calculation
```cpp
auto model = ActivityModelPerplexGFSM();
auto phase = Phase(species).setActivityModel(model);
auto props = phase.evaluate(T, P, x);
// Done! You have fugacity coefficients and molar volume
```

### Workflow B: With Reference State
```cpp
auto model = ActivityModelPerplexGFSM();
auto stdthermo = StandardThermoModelPerplexGFSM(params);
auto phase = Phase(species).setActivityModel(model);
// Use both for complete thermodynamics
```

### Workflow C: Validate with Pure EOS
```cpp
// Get GFSM results
auto model = ActivityModelPerplexGFSM();
auto props_gfsm = /* evaluate */;

// Compare with pure CO2
double vol = 50.0;
double ln_f_pure = PerpleX::crkCO2(P_bar, T_K, vol, /* unused */);
// Check consistency
```

### Workflow D: EOS Comparison Research
```cpp
// Compare different EOS for H2O
double vol = 50.0;
double ln_f_hsmrk = PerpleX::hsmrkf(vol, 1, P_bar, T_K);
vol = 50.0;
double ln_f_cork = 0.0;
PerpleX::crkH2O(P_bar, T_K, vol, ln_f_cork);
// Analyze differences
```

---

## Key Features at a Glance

✅ **13 GFSM species supported** (Model Type 39)
✅ **Three callable interfaces** for different use cases
✅ **Automatic unit conversion** (Pa↔bar, m³/mol↔cm³/mol)
✅ **Hybrid EOS selection** for H2O, CO2, CH4
✅ **Species validation** at runtime
✅ **Matches Reaktoro patterns** (DEW model style)
✅ **Fully documented** with examples
✅ **Production ready** code

---

## Troubleshooting

### "Species 'X' is not a valid Perple_X GFSM fluid species"
**Fix**: Check you're using exactly the 13 allowed species:
H2O, CO2, CH4, H2, CO, H2S, SO2, O2, N2, NH3, HF, C2H6, HCl

### "Temperature X exceeds maximum Y K"
**Fix**: Check `StandardThermoModelParamsPerplexGFSM::Tmax` setting (default 1200 K)

### "Volume not converging in hsmrkf"
**Fix**: Try initial guess closer to expected value, check `PerpleXPureEosOptions` tolerance settings

### Results don't match Perple_X output
**Fix**: Verify (T, P) units are correct, check species order matches Perple_X input

---

## Performance Tips

1. **Reuse ActivityModel**: Create once, use many times
   ```cpp
   auto model = ActivityModelPerplexGFSM();  // Create once
   for (/* many conditions */) {
       auto props = phase.evaluate(T, P, x);  // Reuse model
   }
   ```

2. **Pre-allocate arrays**: For repeated evaluations
   ```cpp
   ActivityProps props;
   for (/* many conditions */) {
       phase.evaluate(props, T, P, x);  // Reuse props structure
   }
   ```

3. **Choose appropriate EOS**: H2O CORK slower than HSMRK but more accurate at high P

---

## Documentation

- **Full API Reference**: See `GFSM_API_REFERENCE.md`
- **Code Examples**: See `GFSM_USAGE_GUIDE.cpp`
- **Implementation Details**: See `IMPLEMENTATION_SUMMARY.md`
- **File Organization**: See `FILE_MANIFEST.md`

---

## Quick Links

| Document | Purpose |
|----------|---------|
| `GFSM_API_REFERENCE.md` | Complete technical specification |
| `GFSM_USAGE_GUIDE.cpp` | Runnable code examples |
| `IMPLEMENTATION_SUMMARY.md` | Architecture and design decisions |
| `FILE_MANIFEST.md` | File organization and integration |
| `DELIVERY_CHECKLIST.md` | What was delivered and verification |

---

## Version Info

**Perple_X**: Type 39 GFSM (Generic Fluid Solution Model)
**Reaktoro**: 1.x compatible
**Status**: ✅ Production Ready

---
