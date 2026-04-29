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
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Extensions/DEW/DEWDatabase.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelDEW.hpp>

namespace Reaktoro {

TEST_CASE("Testing ActivityModelDEW", "[ActivityModelDEW]")
{
    // Create a DEW database
    auto db = DEWDatabase("DeepEarthWater");

    // Create an aqueous phase with common ions and neutral species
    AqueousPhase aqueous("H2O(aq) H+(aq) OH-(aq) Na+(aq) Cl-(aq) SiO2(aq)");
    aqueous.setActivityModel(ActivityModelDEW());

    // Create the chemical system
    ChemicalSystem system(db, aqueous);

    const auto setReferenceComposition = [](ChemicalState& state)
    {
        state.set("H2O(aq)", 55.5, "mol");
        state.set("H+(aq)", 1.0e-7, "mol");
        state.set("OH-(aq)", 1.0e-7, "mol");
        state.set("Na+(aq)", 0.1, "mol");
        state.set("Cl-(aq)", 0.1, "mol");
    };

    WHEN("Using ActivityModelDEW at neutral conditions (298 K, 1 bar)")
    {
        // Create a chemical state
        ChemicalState state(system);
        state.temperature(298.15, "K");
        state.pressure(1e5, "Pa");

        // Set composition (dilute NaCl solution near neutral conditions)
        setReferenceComposition(state);

        const auto& props = state.props().phaseProps("AqueousPhase");
        const auto lng = props.speciesActivityCoefficientsLn();

        REQUIRE(lng.size() == aqueous.species().size());

        for(Index i = 0; i < lng.size(); ++i)
            CHECK(std::isfinite(static_cast<double>(lng[i])));
    }

    WHEN("Using ActivityModelDEW at elevated conditions (473 K, 50 MPa)")
    {
        // Create a chemical state
        ChemicalState state(system);
        state.temperature(473.15, "K");
        state.pressure(5e7, "Pa");

        setReferenceComposition(state);

        const auto& props = state.props().phaseProps("AqueousPhase");
        const auto lng = props.speciesActivityCoefficientsLn();

        REQUIRE(lng.size() == aqueous.species().size());

        for(Index i = 0; i < lng.size(); ++i)
            CHECK(std::isfinite(static_cast<double>(lng[i])));
    }

    WHEN("Using ActivityModelDEW at deep-Earth conditions (573 K, 500 MPa)")
    {
        // Create a chemical state
        ChemicalState state(system);
        state.temperature(573.15, "K");
        state.pressure(5e8, "Pa");

        setReferenceComposition(state);
        state.set("SiO2(aq)", 0.01, "mol");

        const auto& props = state.props().phaseProps("AqueousPhase");
        const auto lng = props.speciesActivityCoefficientsLn();

        REQUIRE(lng.size() == aqueous.species().size());

        for(Index i = 0; i < lng.size(); ++i)
            CHECK(std::isfinite(static_cast<double>(lng[i])));
    }

    WHEN("Verifying water activity from ideal mixing")
    {
        ChemicalState state(system);
        state.temperature(298.15, "K");
        state.pressure(1e5, "Pa");

        // Pure water
        state.set("H2O(aq)", 55.5, "mol");

        const auto& props = state.props().phaseProps("AqueousPhase");
        const auto lna = props.speciesActivitiesLn();

        // Water activity should be close to 1 for pure water
        CHECK(std::isfinite(static_cast<double>(lna[0])));
    }
}

} // namespace Reaktoro
