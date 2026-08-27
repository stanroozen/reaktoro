# GFSM Usage Guide

## Using Perple_X GFSM Model in Reaktoro

This demonstrates the three callable interfaces:

1. `ActivityModelPerplexGFSM()` — GFSM fluid mixture model (explicit speciation)
2. `StandardThermoModelPerplexGFSM(params)` — Standard thermo per species
3. Individual pure EOS functions — `hsmrkf`, `crkH2O`, `crkCO2`, `pseos`, `brmrk`, `haar`, `zhdh2o`, `zd09pr`

---

## Interface 1: GFSM Model as ActivityModel (Mixture in Explicit Speciation)

GFSM requires the user to specify mole fractions of all 13 fluid species:
`H2O`, `CO2`, `CH4`, `H2`, `CO`, `H2S`, `SO2`, `O2`, `N2`, `NH3`, `HF`, `C2H6`, `HCl`

```cpp
#include <Reaktoro/Extensions/Perple_X.hpp>

auto Example_ActivityModelGFSM()
{
    // Step 1: Create species list (must include exactly the 13 GFSM species)
    Reaktoro::SpeciesList species;
    species.push_back({
        .name = "H2O",
        .formula = Reaktoro::Formula("H2O"),
        .aggregateState = Reaktoro::AggregateState::Gas,
        .charge = 0
    });
    species.push_back({
        .name = "CO2",
        .formula = Reaktoro::Formula("CO2"),
        .aggregateState = Reaktoro::AggregateState::Gas,
        .charge = 0
    });
    // ... add remaining 11 species: CH4, H2, CO, H2S, SO2, O2, N2, NH3, HF, C2H6, HCl

    // Step 2: Create GFSM activity model with default options
    auto actmodel = Reaktoro::ActivityModelPerplexGFSM();

    // Step 3: Create Phase with the model
    auto phase = Reaktoro::Phase(species).setActivityModel(actmodel);

    // Step 4: Evaluate properties at T, P, composition
    double T = 700.0;    // Kelvin
    double P = 1.0e8;    // Pascal
    Reaktoro::ArrayXr x = Reaktoro::ArrayXr::Constant(13, 1.0/13.0);  // uniform mole fractions

    // Activity model computes:
    // - ln(fugacity coefficient) for each species
    // - Partial molar volumes
    // - Mixture molar volume
    // (Gibbs energy and enthalpy derivatives not yet implemented in GFSM)
}
```

---

## Interface 2: Standard Thermodynamic Model (Per-Species Properties)

Use `StandardThermoModelPerplexGFSM` when you need reference state properties for individual species.

```cpp
auto Example_StandardThermoModel()
{
    // Configure parameters for H2O
    Reaktoro::StandardThermoModelParamsPerplexGFSM h2o_params;
    h2o_params.G0 = -235000.0;  // Reference Gibbs energy (J/mol) at 298K, 1bar
    h2o_params.H0 = -286000.0;  // Reference enthalpy (J/mol)
    h2o_params.V0 = 1.89e-5;    // Reference molar volume (m³/mol)
    h2o_params.Tmax = 1200.0;   // Maximum temperature (K)

    // Create standard thermo model
    auto stdthermo = Reaktoro::StandardThermoModelPerplexGFSM(h2o_params);

    // Evaluate at T, P
    double T = 700.0;  // K
    double P = 1.0e8;  // Pa
    auto props = stdthermo(T, P);

    // props contains: G0, H0, V0, Cp0, VT0, VP0
    // These are used as baseline reference properties for the species
}
```

---

## Interface 3: Individual Pure EOS Functions (Direct Callable)

Exposed through `PerpleXPureEos.hpp`. Allows direct computation of ln(fugacity) and molar volume for pure species at any P-T point.

### Available functions

| Function   | EOS              | Species              |
|------------|------------------|----------------------|
| `hsmrkf`   | HSMRK            | H2O, CO2, CH4        |
| `crkH2O`   | CORK             | H2O only             |
| `crkCO2`   | CORK             | CO2 only             |
| `pseos`    | PSEOS            | H2O, CO2             |
| `brmrk`    | BRMRK            | H2O, CO2             |
| `haar`     | Haar empirical   | H2O                  |
| `zhdh2o`   | Zhang-Duan 2005  | H2O                  |
| `zd09pr`   | Zhang-Duan 2009  | All species          |

```cpp
auto Example_PureEOS()
{
    using namespace Reaktoro::PerpleX;

    // Example: Evaluate H2O using HSMRK
    {
        double volume_cm3_per_mol = 50.0;  // Initial guess
        int species_idx = 1;  // H2O
        double P_bar = 1000.0;
        double T_K = 700.0;
        PerpleXPureEosOptions opts;

        double ln_fugacity = hsmrkf(volume_cm3_per_mol, species_idx, P_bar, T_K, opts);
        // volume_cm3_per_mol is updated to converged value
    }

    // Example: Evaluate H2O using CORK
    {
        double P_bar = 1000.0;
        double T_K = 700.0;
        double volume_cm3_per_mol = 0.0;
        double ln_fugacity = 0.0;

        crkH2O(P_bar, T_K, volume_cm3_per_mol, ln_fugacity);
        // Both volume and ln_fugacity are computed
    }

    // Example: Evaluate CO2 using CORK
    {
        double P_bar = 1000.0;
        double T_K = 700.0;
        double volume_cm3_per_mol = 0.0;
        double ln_fugacity = 0.0;

        crkCO2(P_bar, T_K, volume_cm3_per_mol, ln_fugacity);
    }

    // Example: Using hybrid EOS selector to switch between alternatives
    {
        auto options = makePerpleXHybridEosOptions();
        // options.waterEos: MRK, HSMRK, CORK, PSEOS, HAAR, ZDHH2O, ZD09
        // options.co2Eos:   MRK, HSMRK, CORK, BRMRK, PSEOS, ZD09
        // options.ch4Eos:   MRK, HSMRK, ZD09

        HybridEosResult result = hybEos(
            species_idx,   // 1, 2, or 4 for H2O, CO2, CH4
            ln_f_mrk,      // baseline MRK ln(f)
            g_mrk,         // baseline MRK fugacity coefficient
            v_mrk,         // baseline MRK volume
            P_bar,
            T_K,
            options
        );
        // result.ln_f, result.g, result.v updated with hybrid choice
    }
}
```

---

## Complete Workflow: All Three Interfaces Together

```cpp
auto Example_Complete_Workflow()
{
    // Part A: Setup GFSM mixture (13 species)
    Reaktoro::SpeciesList species;
    // ... populate species ...

    auto gfsm_actmodel = Reaktoro::ActivityModelPerplexGFSM();
    auto gfsm_phase = Reaktoro::Phase(species).setActivityModel(gfsm_actmodel);

    // Part B: Setup standard thermo for reference state
    Reaktoro::StandardThermoModelParamsPerplexGFSM h2o_params;
    auto h2o_stdthermo = Reaktoro::StandardThermoModelPerplexGFSM(h2o_params);

    // Part C: Direct evaluation of pure EOS for validation
    double P_bar = 1000.0;
    double T_K = 700.0;
    double vol_cm3 = 50.0;
    double ln_f_hsmrk = Reaktoro::PerpleX::hsmrkf(vol_cm3, 1, P_bar, T_K);

    // Now you have:
    // 1. GFSM mixture at any composition (mixture level)
    // 2. Standard thermo baseline (species level)
    // 3. Pure EOS evaluation (validation/coupling level)
}
```

---

## Key Design Decisions

### 1. Why GFSM uses explicit speciation
- Direct specification of all 12 independent mole fractions
- No implicit speciation (unlike aqueous models)
- Each species evaluated independently at pure-species level, then combined without a mixing law

### 2. Why three callables
- `ActivityModelPerplexGFSM()` — Standard Reaktoro integration (mixture level)
- `StandardThermoModelPerplexGFSM()` — Reference state per species
- Individual pure EOS — Direct calculation for validation/research

### 3. Hybrid EOS selection
| Species | Available options                        |
|---------|------------------------------------------|
| H2O     | HSMRK, CORK, PSEOS, Haar, ZD05, ZD09, MRK |
| CO2     | HSMRK, CORK, BRMRK, PSEOS, ZD09, MRK    |
| CH4     | HSMRK, ZD09, MRK                         |
| Others  | Fixed to MRK (no alternatives in Perple_X data) |

### 4. Unit conventions
| Scope              | Pressure | Volume       |
|--------------------|----------|--------------|
| Reaktoro interface | Pa       | m³/mol       |
| Perple_X internal  | bar      | cm³/mol      |

Conversions are handled internally by wrapper functions.
