#pragma once

#include <stdexcept>

#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Species.hpp>

namespace Reaktoro::TestUtils {

inline auto reorderPilotMixedConditionsSpeciesAmounts(
    ChemicalSystem const& system,
    SpeciesList const& phaseSpecies,
    ArrayXrConstRef logicalAmounts,
    const char* context) -> ArrayXr
{
    const Strings aqueousSpecies = {"H2O(aq)", "H+(aq)", "OH-(aq)", "Na+(aq)", "Cl-(aq)"};
    const auto expectedSize = static_cast<Index>(aqueousSpecies.size() + phaseSpecies.size() + 1);
    if(logicalAmounts.size() != expectedSize)
        throw std::runtime_error(context);

    ArrayXr reordered = ArrayXr::Zero(system.species().size());
    Index offset = 0;

    for(const auto& name : aqueousSpecies)
        reordered[system.species().indexWithName(name)] = logicalAmounts[offset++];

    for(Index i = 0; i < phaseSpecies.size(); ++i)
        reordered[system.species().indexWithName(phaseSpecies[i].name())] = logicalAmounts[offset++];

    reordered[system.species().indexWithName("SiO2(s)")] = logicalAmounts[offset];
    return reordered;
}

} // namespace Reaktoro::TestUtils