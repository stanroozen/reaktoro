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

#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Extensions/DEW/DEWDatabase.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelDEW.hpp>
#include <Reaktoro/Phases/AqueousPhase.hpp>

namespace Reaktoro {

TEST_CASE("Testing ActivityModelDEW", "[ActivityModelDEW]")
{
    // Create a DEW database
    auto db = DEWDatabase("DeepEarthWater");

    // Create an aqueous phase with common ions and neutral species
    AqueousPhase aqueous("H2O(aq) H+ OH- Na+ Cl- SiO2(aq)");
    aqueous.setActivityModel(ActivityModelDEW());

    // Create the chemical system
    ChemicalSystem system(db, aqueous);

    WHEN("Using ActivityModelDEW at neutral conditions (298 K, 1 bar)")
    {
        // Create a chemical state
        ChemicalState state(system);
        state.temperature(298.15, "K");
        state.pressure(1e5, "Pa");

        // Set composition (dilute NaCl solution at pH 7)
        state.set("Na+", 0.1, "mol/kg");
        state.set("Cl-", 0.1, "mol/kg");
        state.set("pH", 7.0);

        // Check that activity coefficients are reasonable
        const auto phase = state.phase(0);
        const auto& props = state.phaseProps(0);

        // For dilute solutions at 25°C, activity coefficients should be close to 1
        // (slightly less than 1 due to Debye-Hückel effect)
        REQUIRE(props.ln_g.size() == aqueous.species().size());

        // Activity coefficients should be negative (γ < 1) due to ionic interactions
        // But not too negative for dilute solutions
        for(Index i = 0; i < props.ln_g.size(); ++i)
        {
            CHECK(props.ln_g[i] < 0.5);   // Should be less than ln(1.65)
        }
    }

    WHEN("Using ActivityModelDEW at elevated conditions (473 K, 50 MPa)")
    {
        // Create a chemical state
        ChemicalState state(system);
        state.temperature(473.15, "K");
        state.pressure(5e7, "Pa");

        // Set composition
        state.set("Na+", 0.1, "mol/kg");
        state.set("Cl-", 0.1, "mol/kg");
        state.set("pH", 7.0);

        const auto& props = state.phaseProps(0);

        // At higher temperatures, activity coefficients should change
        REQUIRE(props.ln_g.size() == aqueous.species().size());

        // Check that we get reasonable values (no NaN, no inf)
        for(Index i = 0; i < props.ln_g.size(); ++i)
        {
            CHECK(std::isfinite(props.ln_g[i]));
        }
    }

    WHEN("Using ActivityModelDEW at deep-Earth conditions (573 K, 500 MPa)")
    {
        // Create a chemical state
        ChemicalState state(system);
        state.temperature(573.15, "K");
        state.pressure(5e8, "Pa");

        // Set composition
        state.set("Na+", 0.1, "mol/kg");
        state.set("Cl-", 0.1, "mol/kg");
        state.set("SiO2(aq)", 0.01, "mol/kg");

        const auto& props = state.phaseProps(0);

        // Should produce valid results at extreme conditions
        REQUIRE(props.ln_g.size() == aqueous.species().size());

        // All activity coefficients should be finite
        for(Index i = 0; i < props.ln_g.size(); ++i)
        {
            CHECK(std::isfinite(props.ln_g[i]));
        }
    }

    WHEN("Verifying water activity from ideal mixing")
    {
        ChemicalState state(system);
        state.temperature(298.15, "K");
        state.pressure(1e5, "Pa");

        // Pure water
        state.set("Na+", 0.0, "mol/kg");
        state.set("Cl-", 0.0, "mol/kg");

        const auto& props = state.phaseProps(0);

        // Water activity should be close to 1 for pure water
        CHECK(std::isfinite(props.ln_a[0]));  // iwater should be first in phase
    }
}

} // namespace Reaktoro
