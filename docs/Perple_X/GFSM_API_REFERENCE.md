# Perple_X GFSM Integration for Reaktoro - Complete API Reference

## Overview

This document describes the complete Reaktoro integration for Perple_X GFSM (Generic Fluid Solution Model, Type 39). The integration provides three complementary callable interfaces for different use cases.

## Model Details

**Perple_X Type 39: GFSM (Generic Fluid Solution Model)**

- **Speciation Space**: Explicit (user provides all 13 mole fractions)
- **Framework**: Modified Redlich-Kwong (MRK) baseline for all species
- **Pure EOS Options**: Hybrid alternatives for H2O, CO2, CH4 only
- **13 Allowed Species**:
  - H2O, CO2, CH4, H2, CO, H2S, SO2, O2, N2, NH3, HF, C2H6, HCl
  - (Updated from original 18, restricting to Type 39 allowed species only)

## Callable Interfaces

### 1. ActivityModelPerplexGFSM()

**Location**: `Reaktoro/Extensions/Perple_X/ActivityModelPerplexGFSM.hpp`

**Purpose**: Mixture-level activity model for GFSM fluid phases.

**Signature**:
```cpp
auto ActivityModelPerplexGFSM(
    const ActivityModelParamsPerplexGFSM& params = {}) -> ActivityModelGenerator;
```

**Parameters** (`ActivityModelParamsPerplexGFSM`):
- `hybridEosOptions`: HybridEosOptions for selecting pure EOS (H2O, CO2, CH4)
- `mrkMixOptions`: MRK mixture calculation settings
- `useLowTMrk`: Boolean flag for low-temperature MRK variant
- `enableElectrolyte`: Boolean flag for electrolyte solvent properties
- `pureEosOptions`: PerpleXPureEosOptions (convergence tolerances)

**Usage**:
```cpp
SpeciesList species = {H2O, CO2, CH4, H2, CO, H2S, SO2, O2, N2, NH3, HF, C2H6, HCl};
auto actmodel = ActivityModelPerplexGFSM();
auto phase = Phase(species).setActivityModel(actmodel);

// Evaluate properties at T, P, composition
ActivityProps props = phase.evaluate(T, P, x);
// Provides: ln_g (log activity coefficients), Vx (molar volume), ln_a (log activities)
```

**Features**:
- Returns `ActivityModelGenerator` (callable accepting `SpeciesList`)
- Automatically maps Reaktoro species names to Perple_X indices (1-13)
- Validates that input species are GFSM-allowed
- Converts units internally (Pa→bar, m³/mol→cm³/mol)
- Computes mixture properties in explicit speciation space

### 2. StandardThermoModelPerplexGFSM()

**Location**: `Reaktoro/Extensions/Perple_X/StandardThermoModelPerplexGFSM.hpp`

**Purpose**: Standard thermodynamic properties (reference state) for individual species.

**Signature**:
```cpp
auto StandardThermoModelPerplexGFSM(
    const StandardThermoModelParamsPerplexGFSM& params) -> StandardThermoModel;
```

**Parameters** (`StandardThermoModelParamsPerplexGFSM`):
- `G0`: Standard molar Gibbs energy at reference state (J/mol)
- `H0`: Standard molar enthalpy at reference state (J/mol)
- `V0`: Standard molar volume at reference state (m³/mol)
- `hybridEosOptions`: For context (though per-species model)
- `Tmax`: Maximum applicable temperature (K), default 1200

**Usage**:
```cpp
StandardThermoModelParamsPerplexGFSM h2o_params;
h2o_params.G0 = -235000.0;  // J/mol
h2o_params.H0 = -286000.0;  // J/mol
h2o_params.V0 = 1.89e-5;    // m³/mol
h2o_params.Tmax = 1200.0;   // K

auto stdthermo = StandardThermoModelPerplexGFSM(h2o_params);
StandardThermoProps props = stdthermo(T, P);
// Provides: G0, H0, V0, Cp0 (derivatives VT0, VP0 at zero)
```

**Features**:
- Returns `StandardThermoModel` (callable accepting T, P)
- Provides reference state properties for consistency
- Used alongside ActivityModelPerplexGFSM for complete thermodynamics
- Can be instantiated per species with different parameter sets

### 3. Individual Pure EOS Functions

**Location**: `Reaktoro/Extensions/Perple_X/PerpleXPureEos.hpp`

**Purpose**: Direct pure-species equation of state evaluation (lowest level).

**Available Functions**:

#### Water (H2O, species index 1):
- `hsmrkf()` - HSMRK (Kerrick & Jacobs 1981)
- `crkH2O()` - CORK (Holland & Powell 1990)
- `pseos()` - PSEOS
- `brmrk()` - BRMRK
- `haar()` - Haar empirical
- `zhdh2o()` - Zhang-Duan 2005
- `zd09pr()` - Zhang-Duan 2009 (Peng-Robinson variant)

#### CO2 (species index 2):
- `hsmrkf()` - HSMRK
- `crkCO2()` - CORK
- `pseos()` - PSEOS
- `brmrk()` - BRMRK
- `zd09pr()` - Zhang-Duan 2009

#### CH4 (species index 4):
- `hsmrkf()` - HSMRK
- `zd09pr()` - Zhang-Duan 2009

#### All other species (indices 3, 5, 6, 7, 8, 10, 11, 16, 17, 18):
- **MRK only** (no alternatives in Perple_X data)

**Generic Signatures**:
```cpp
// HSMRK variant
double hsmrkf(double& vol, int species, double pressureBar, double temperatureK,
              const PerpleXPureEosOptions& options = {});

// CORK variants
void crkH2O(double P_bar, double T_K, double& vol, double& lnfug);
void crkCO2(double P_bar, double T_K, double& vol, double& lnfug);

// Others follow similar patterns
```

**Usage Example**:
```cpp
using namespace Reaktoro::PerpleX;

// Evaluate H2O with HSMRK
double volume = 50.0;  // cm³/mol initial guess
int species = 1;       // H2O
double P = 1000.0;     // bar
double T = 700.0;      // K
PerpleXPureEosOptions opts;

double ln_f = hsmrkf(volume, species, P, T, opts);
// volume updated to converged value, ln_f is ln(fugacity coefficient)
```

**Hybrid EOS Selector**:
```cpp
auto options = makePerpleXHybridEosOptions();
options.waterEos = HybridEosOptions::WaterEos::CORK;
options.co2Eos = HybridEosOptions::CO2Eos::BRMRK;
options.ch4Eos = HybridEosOptions::CH4Eos::HSMRK;

HybridEosResult result = hybEos(species, ln_f_mrk, g_mrk, v_mrk, P, T, options);
```

## Hierarchical Architecture

```
Reaktoro User Code
       │
       ├─→ ActivityModelPerplexGFSM()      [Mixture level]
       │   └─→ Returns ActivityModelGenerator
       │       └─→ Accepts SpeciesList
       │           └─→ Returns ActivityModel
       │               └─→ Computes ln_g, Vx at (T,P,x)
       │
       ├─→ StandardThermoModelPerplexGFSM() [Species reference state]
       │   └─→ Returns StandardThermoModel
       │       └─→ Computes G0, H0, V0 at (T,P)
       │
       └─→ Individual Pure EOS              [Direct calculation]
           (hsmrkf, crkH2O, crkCO2, ...)
           └─→ Compute ln_f, volume at (T,P,species)
```

## Technical Specifications

### Thermodynamic Quantities Computed

| Quantity | ActivityModel | StandardThermoModel | Pure EOS |
|----------|---------------|-------------------|----------|
| ln(fugacity coefficient) | ✓ (ln_g) | ✗ | ✓ |
| ln(activity) | ✓ (ln_a) | ✗ | implicit |
| Molar volume | ✓ (Vx) | ✓ (V0) | ✓ |
| Gibbs energy | ✗ | ✓ (G0 only) | ✗ |
| Enthalpy | ✗ | ✓ (H0 only) | ✗ |
| Cp | ✗ | ✓ (0 in current impl) | ✗ |

### Unit Conversions

| Quantity | Input | Internal | Output |
|----------|-------|----------|--------|
| Pressure | Pa | bar | Pa |
| Temperature | K | K | K |
| Volume | m³/mol | cm³/mol | m³/mol |
| Energy | J/mol | J/mol | J/mol |

### Species Index Mapping

| Name | Formula | Index | EOS Options | MRK only? |
|------|---------|-------|-------------|-----------|
| H2O | H₂O | 1 | 7 | No |
| CO2 | CO₂ | 2 | 6 | No |
| CO | CO | 3 | 1 | Yes |
| CH4 | CH₄ | 4 | 3 | No |
| H2 | H₂ | 5 | 1 | Yes |
| H2S | H₂S | 6 | 1 | Yes |
| O2 | O₂ | 7 | 1 | Yes |
| SO2 | SO₂ | 8 | 1 | Yes |
| N2 | N₂ | 10 | 1 | Yes |
| NH3 | NH₃ | 11 | 1 | Yes |
| HF | HF | 17 | 1 | Yes |
| C2H6 | C₂H₆ | 16 | 1 | Yes |
| HCl | HCl | 18 | 1 | Yes |

*Note: Indices skip 9, 12, 13, 14, 15, 19 (unused in Type 39)*

## File Structure

```
Reaktoro/Extensions/Perple_X/
├── Perple_X.hpp                          [Main extension header]
├── PerpleXSpecies.hpp                    [Species enum (13 species)]
├── PerpleXMrkParameters.hpp              [MRK parameters]
├── PerpleXMrkMixture.hpp                 [MRK mixture rules]
├── PerpleXMrkPure.hpp                    [MRK pure species]
├── PerpleXHybridEos.hpp                  [Hybrid EOS selector]
├── PerpleXFluidModel.hpp/cpp             [Legacy interface]
├── PerpleXGFSMModel.hpp/cpp              [GFSM core compute]
├── PerpleXPureEos.hpp                    [Individual pure EOS functions]
├── ActivityModelPerplexGFSM.hpp/cpp      [NEW: GFSM activity model]
├── StandardThermoModelPerplexGFSM.hpp/cpp [NEW: Standard thermo model]
└── GFSM_USAGE_GUIDE.cpp                  [NEW: Usage examples]
```

## Validation & Testing Recommendations

1. **ActivityModelPerplexGFSM**:
   - Validate fugacity coefficients against Perple_X output at reference P-T
   - Check molar volume consistency with pure species volumes

2. **StandardThermoModelPerplexGFSM**:
   - Verify reference state energies match input literature values
   - Check temperature sensitivity against HKF-style models

3. **Individual Pure EOS**:
   - Compare with published EOS papers (Kerrick & Jacobs, Holland & Powell, etc.)
   - Verify convergence behavior (Newton-Raphson iterations)

## Backward Compatibility

- All existing code using PerpleXFluidModel remains unchanged
- New callables are additions only (no deletions)
- PerpleXSpecies.hpp modified (18→13 species), affects:
  - `PerpleXSpecies::speciesCount()` returns 13
  - Any iteration over all species must account for 13 only

## Future Extensions

Potential enhancements:
1. Electrolyte module (aq) coupling for dissolved species
2. Gibbs-Helmholtz integration for enthalpy/entropy from fugacity
3. Dielectric constant evolution for activity coefficient corrections
4. Phase stability analysis (spinodal/binodal)

## References

- Perple_X Documentation: http://www.perplex.ethz.ch/
- Kerrick & Jacobs (1981) "A modified equation of state for CO2-rich fluid mixtures" GCA 45:629-641
- Holland & Powell (1990) "An improved and extended internally consistent thermodynamic dataset for phases of interest to high pressure, high-temperature metamorphism and core formation" JMG 8:89-124
- Zhang & Duan (2005, 2009) EOS papers in Journal of Supercritical Fluids
