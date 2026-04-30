# DEW Activity Model Implementation Guide

## Files to Create/Modify

### 1. `Reaktoro/Extensions/DEW/DEWActivityModel.hpp` (NEW)

```cpp
// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

// Forward declarations
class DEWDatabase;
struct WaterModelOptions;

/// Create an activity model for aqueous species using the DEW (Deep Earth Water) model.
///
/// The DEW activity model calculates activities and activity coefficients for aqueous species
/// using the HKF (Helgeson-Kirkham-Flowers) model with DEW-specific electrostatic calculations
/// and temperature-pressure dependent equations of state for water.
///
/// @param species The species in the aqueous phase
/// @param database The DEW database with thermodynamic data for aqueous species
/// @return An ActivityModelGenerator that creates the activity model
///
/// @see DEWDatabase, WaterThermoModel, WaterElectroModel
/// @ingroup Extensions
auto ActivityModelDEW(const SpeciesList& species, const DEWDatabase& database) -> ActivityModelGenerator;

/// Create an activity model for aqueous species using the DEW model with custom options.
///
/// @param species The species in the aqueous phase
/// @param database The DEW database with thermodynamic data for aqueous species
/// @param options Custom options for water property calculations
/// @return An ActivityModelGenerator that creates the activity model
auto ActivityModelDEW(
    const SpeciesList& species,
    const DEWDatabase& database,
    const WaterModelOptions& options) -> ActivityModelGenerator;

} // namespace Reaktoro
```

### 2. `Reaktoro/Extensions/DEW/DEWActivityModel.cpp` (NEW)

Key implementation points (pseudo-code structure):

```cpp
#include "DEWActivityModel.hpp"

// Include all necessary DEW components
#include <Reaktoro/Extensions/DEW/DEWDatabase.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoModel.hpp>
#include <Reaktoro/Extensions/DEW/WaterElectroModel.hpp>
#include <Reaktoro/Extensions/DEW/WaterGibbsModel.hpp>
#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>
#include <Reaktoro/Core/AqueousMixture.hpp>

namespace Reaktoro {

namespace {

/// Helper function to get HKF Born coefficient (omega) for a species
auto calcBornCoefficient(const Species& species, real effectiveRadius) -> real
{
    const auto eta = 1.66027e5; // Born parameter (J*K/mol)
    const auto charge = species.charge();
    const auto omegaH = 0.5387e+05; // Born coefficient of H+
    return eta * charge * charge / effectiveRadius - charge * omegaH;
}

/// Helper function to get effective ionic radius
auto effectiveIonicRadius(const Species& species) -> real
{
    // Similar to ActivityModelHKF implementation
    // Use charge-based estimation if specific radius not available
    const auto Zi = species.charge();
    if(Zi == -1) return 1.81;      // based on Cl-
    if(Zi == -2) return 3.00;      // based on CO3--, SO4--
    if(Zi == +1) return 2.31;      // based on NH4+
    if(Zi == +2) return 2.80;      // based on +2 species average
    if(Zi == +3) return 3.60;      // based on +3 species average
    return 0.0; // Will use formula-based defaults
}

} // namespace

auto ActivityModelDEW(const SpeciesList& species, const DEWDatabase& database) -> ActivityModelGenerator
{
    WaterModelOptions options;
    return ActivityModelDEW(species, database, options);
}

auto ActivityModelDEW(
    const SpeciesList& species,
    const DEWDatabase& database,
    const WaterModelOptions& options) -> ActivityModelGenerator
{
    return [species, database, options](const SpeciesList& specieslist)
    {
        // Create aqueous mixture wrapper
        AqueousMixture mixture(specieslist);

        // Cache indices and properties
        const auto iwater = mixture.indexWater();
        const auto icharged = mixture.indicesCharged();
        const auto ineutral = mixture.indicesNeutral();

        // Extract DEW thermo and electro models from database
        // (These will be retrieved based on species in the database)

        // Initialize water models
        WaterThermoModel thermoModel(options);
        WaterElectroModel electroModel(options);

        // Shared pointers for caching
        auto stateptr = std::make_shared<AqueousMixtureState>();
        auto mixtureptr = std::make_shared<AqueousMixture>(mixture);

        // Return the actual activity model function
        ActivityModel model = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
        {
            const auto& [T, P, x] = args;

            // 1. Calculate water properties using DEW models
            WaterState wstate = thermoModel(T, P);
            auto wthermoprops = wstate.thermoProps();  // rho, kappa, etc.
            auto welectroprops = electroModel(T, P);   // epsilon, DH_A, DH_B, Born Q, etc.

            // 2. Evaluate aqueous mixture state
            auto const& mixstate = *stateptr = mixture.state(T, P, x);
            const auto& I = mixstate.Is;       // ionic strength
            const auto& m = mixstate.m;        // molalities
            const auto& ms = mixstate.ms;      // stoichiometric molalities

            // 3. Cache water properties in extras
            props.extra["WaterThermoProps"] = std::make_shared<decltype(wthermoprops)>(wthermoprops);
            props.extra["WaterElectroProps"] = std::make_shared<decltype(welectroprops)>(welectroprops);
            props.extra["AqueousMixtureState"] = stateptr;
            props.extra["AqueousMixture"] = mixtureptr;

            // 4. Calculate HKF activity coefficients
            const auto sqrtI = sqrt(I);
            const auto xw = x[iwater];
            const auto ln_xw = log(xw);

            // 4a. Neutral species: ln_g = b * I
            for(auto idx : ineutral)
            {
                const auto b = 0.1; // PHREEQC default
                props.ln_g[idx] = std::log(10.0) * b * I;
            }

            // 4b. Charged species: use HKF equation with DEW water properties
            real osmotic_coeff = 0.0;
            const auto A = welectroprops.DH_A;  // Debye-Hückel parameter A
            const auto B = welectroprops.DH_B;  // Debye-Hückel parameter B
            const auto eta = 1.66027e5;         // Born parameter

            for(auto idx : icharged)
            {
                const auto& sp = mixture.species(idx);
                const auto charge = sp.charge();
                const auto eff_radius = effectiveIonicRadius(sp);

                // Born coefficient
                const auto omega = calcBornCoefficient(sp, eff_radius);

                // Debye-Hückel ion size parameter
                const auto a = 2.0 * (eff_radius + 1.81 * abs(charge)) / (abs(charge) + 1.0);

                // Lambda parameter
                const auto lambda = 1.0 + a * B * sqrtI;

                // Activity coefficient (HKF equation)
                const auto log10_g = -(A * charge * charge * sqrtI) / lambda
                                   + std::log10(xw)
                                   + omega * I; // Simplified Born term

                props.ln_g[idx] = log10_g * std::log(10.0);
            }

            // 5. Calculate activities (ln_a = ln_g + ln_m)
            props.ln_a = props.ln_g + m.log();

            // 6. Activity of water
            if(xw != 1.0)
                props.ln_a[iwater] = std::log(10.0) * mixture.water().molarMass() * osmotic_coeff;
            else
                props.ln_a[iwater] = ln_xw;

            props.ln_g[iwater] = props.ln_a[iwater] - ln_xw;

            // 7. Set state of matter
            props.som = StateOfMatter::Liquid;
        };

        return model;
    };
}

} // namespace Reaktoro
```

### 3. Update `Reaktoro/Extensions/DEW.hpp`

Add the export after existing DEW imports:

```cpp
// DEW Activity Models
#include <Reaktoro/Extensions/DEW/DEWActivityModel.hpp>
```

### 4. Update `Reaktoro/Extensions/DEW/CMakeLists.txt`

Add to the source files list:
```cmake
DEWActivityModel.cpp
DEWActivityModel.hpp
```

And to test files:
```cmake
DEWActivityModel.test.cxx
```

### 5. Create `Reaktoro/Extensions/DEW/DEWActivityModel.test.cxx` (NEW)

```cpp
#include <catch2/catch.hpp>

#include <Reaktoro/Extensions/DEW.hpp>

using namespace Reaktoro;

TEST_CASE("DEWActivityModel", "[DEW]")
{
    // Test 1: Simple aqueous system
    // Test 2: Quartz saturation
    // Test 3: Activity coefficient calculations
    // Test 4: Water property caching
}
```

### 6. Create `Reaktoro/Extensions/DEW/DEWActivityModel.py.cxx` (NEW)

```cpp
// Pybind11 bindings for ActivityModelDEW
#include <pybind11/pybind11.h>
#include <Reaktoro/Extensions/DEW/DEWActivityModel.hpp>

namespace py = pybind11;

PYBIND11_MODULE(reaktoro4py, m)
{
    m.def("ActivityModelDEW",
        py::overload_cast<const Reaktoro::SpeciesList&, const Reaktoro::DEWDatabase&>
        (&Reaktoro::ActivityModelDEW),
        py::arg("species"), py::arg("database"));

    m.def("ActivityModelDEW",
        py::overload_cast<const Reaktoro::SpeciesList&,
                          const Reaktoro::DEWDatabase&,
                          const Reaktoro::WaterModelOptions&>
        (&Reaktoro::ActivityModelDEW),
        py::arg("species"), py::arg("database"), py::arg("options"));
}
```

## Integration with Python Script

After implementation, `quartz_solubility_analysis.py` becomes:

```python
from reaktoro import *

# Create system with DEW activity model
db = SupcrtDatabase("supcrtbl")
dew_db = DEWDatabase("dew2024-aqueous")

# Define aqueous phase with DEW activity model
aqueous = AqueousPhase("H2O(aq) H+ OH- SiO2(aq)")
aqueous.setActivityModel(ActivityModelDEW(aqueous.species(), dew_db))

# Create system
system = ChemicalSystem(db, aqueous, MineralPhases("Quartz"))

# Use solver as normal
solver = EquilibriumSolver(system)

for T_C in temperatures:
    for P_kbar in pressures:
        state = ChemicalState(system)
        state.temperature(T_C, "celsius")
        state.pressure(P_kbar * 100, "bar")
        # ... set initial composition

        result = solver.solve(state)

        # Extract SiO2(aq) molality
        props = EquilibriumProps(state)
        siO2_molality = props.speciesMolality("SiO2(aq)")
```

## Key Design Features

1. **Lazy Evaluation**: Water properties calculated only when needed
2. **Caching**: Results stored in `props.extra` for potential reuse
3. **Modularity**: Can chain with other activity models via `chain()`
4. **Compatibility**: Works seamlessly with `ChemicalSystem` and `Minimizer`
5. **Temperature/Pressure Dependent**: Full T-P dependency via DEW models

## Testing Checklist

- [ ] Code compiles without errors
- [ ] Python bindings work
- [ ] Simple aqueous systems calculate correctly
- [ ] Water properties match direct model calculations
- [ ] Ionic strength calculations are correct
- [ ] Activity coefficients are physically reasonable
- [ ] Quartz solubility script runs to completion
- [ ] Results converge with equilibrium solver
- [ ] No memory leaks or crashes
- [ ] Results agree with experimental data
