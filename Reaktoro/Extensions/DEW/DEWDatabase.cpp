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

#include "DEWDatabase.hpp"

// C++ includes
#include <fstream>
#include <sstream>

// yaml-cpp includes
#include <yaml-cpp/yaml.h>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/Embedded.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelDEW.hpp>
#include <Reaktoro/Serialization/Models/StandardThermoModels.hpp>

namespace Reaktoro {
namespace {

/// Load species from a YAML database content string
auto loadSpeciesFromYAML(const String& yaml_content) -> SpeciesList
{
    SpeciesList species_list;
    Set<String> inserted_elements;
    ElementList elements;

    try {
        YAML::Node root = YAML::Load(yaml_content);

        if(!root["Species"] || !root["Species"].IsMap())
        {
            errorif(true, "DEW database: missing or invalid 'Species' node in YAML");
        }

        const YAML::Node& speciesNode = root["Species"];

        for(auto it = speciesNode.begin(); it != speciesNode.end(); ++it)
        {
            const String speciesName = it->first.as<String>();
            const YAML::Node& spec = it->second;

            // Extract basic species information
            String name = spec["Name"] ? spec["Name"].as<String>() : speciesName;
            String formula = spec["Formula"] ? spec["Formula"].as<String>() : speciesName;
            double charge = spec["Charge"] ? spec["Charge"].as<double>() : 0.0;

            // Determine aggregate state
            String aggStateStr = spec["AggregateState"] ? spec["AggregateState"].as<String>() : "Aqueous";
            AggregateState aggState = AggregateState::Aqueous;
            if(aggStateStr == "Gas" || aggStateStr == "Gaseous")
                aggState = AggregateState::Gas;
            else if(aggStateStr == "Solid" || aggStateStr == "Mineral")
                aggState = AggregateState::Solid;
            else if(aggStateStr == "Liquid")
                aggState = AggregateState::Liquid;

            // Parse elements from formula
            // Handle DEW-specific suffix ",aq" (e.g., "MgO,aq")
            // Replace comma with underscore for species name, strip for formula parsing
            String species_name_suffix = "";
            String clean_formula = formula;
            auto comma_pos = clean_formula.find(",aq");
            if(comma_pos != String::npos)
            {
                species_name_suffix = "_aq";  // Will be used in species name
                clean_formula = clean_formula.substr(0, comma_pos);  // Strip for parsing
            }

            ChemicalFormula chem_formula(clean_formula);
            Pairs<Element, double> element_pairs;

            for(const auto& [symbol, coeff] : chem_formula.elements())
            {
                // Add element if not already present
                if(inserted_elements.find(symbol) == inserted_elements.end())
                {
                    elements.append(Element(symbol));  // Use Element(symbol) constructor - looks up periodic table
                    inserted_elements.insert(symbol);
                }

                // Find the element in our list
                auto idx = elements.find(symbol);
                element_pairs.emplace_back(elements[idx], coeff);
            }

            // Create the species object
            // Use modified name with underscore instead of comma (MgO,aq -> MgO_aq)
            String modified_name = name;
            auto comma_in_name = modified_name.find(",aq");
            if(comma_in_name != String::npos)
                modified_name = modified_name.substr(0, comma_in_name) + "_aq";

            Species species;
            species = species.withName(modified_name)
                             .withFormula(clean_formula)  // Use cleaned formula without ,aq
                             .withElements(element_pairs)
                             .withCharge(charge)
                             .withAggregateState(aggState);

            // Extract HKF parameters if present
            if(spec["StandardThermoModel"] && spec["StandardThermoModel"]["HKF"])
            {
                const YAML::Node& hkf = spec["StandardThermoModel"]["HKF"];

                // Create DEW model parameters (uses same HKF parameters but with DEW water models)
                StandardThermoModelParamsDEW params;
                params.Gf = hkf["Gf"] ? hkf["Gf"].as<double>() : 0.0;
                params.Hf = hkf["Hf"] ? hkf["Hf"].as<double>() : 0.0;
                params.Sr = hkf["Sr"] ? hkf["Sr"].as<double>() : 0.0;
                params.a1 = hkf["a1"] ? hkf["a1"].as<double>() : 0.0;
                params.a2 = hkf["a2"] ? hkf["a2"].as<double>() : 0.0;
                params.a3 = hkf["a3"] ? hkf["a3"].as<double>() : 0.0;
                params.a4 = hkf["a4"] ? hkf["a4"].as<double>() : 0.0;
                params.c1 = hkf["c1"] ? hkf["c1"].as<double>() : 0.0;
                params.c2 = hkf["c2"] ? hkf["c2"].as<double>() : 0.0;
                // wref: the YAML value is expected to be in J/mol already (generation scripts
                // convert from Excel/SUPCRT conventions into SI units). Do not apply any
                // additional ad-hoc scaling here; use the YAML-provided value directly.
                params.wref = hkf["wref"] ? hkf["wref"].as<double>() : 0.0;
                params.charge = charge;
                params.Tmax = hkf["Tmax"] ? hkf["Tmax"].as<double>() : 1000.0; // Default max T
                // waterOptions will use default DEW settings with high-precision integration

                // Attach DEW standard thermo model
                species = species.withStandardThermoModel(StandardThermoModelDEW(params));

                // Store HKF params as attached data for later use
                Data hkf_data;
                hkf_data["Gf"] = params.Gf;
                hkf_data["Hf"] = params.Hf;
                hkf_data["Sr"] = params.Sr;
                hkf_data["a1"] = params.a1;
                hkf_data["a2"] = params.a2;
                hkf_data["a3"] = params.a3;
                hkf_data["a4"] = params.a4;
                hkf_data["c1"] = params.c1;
                hkf_data["c2"] = params.c2;
                hkf_data["wref"] = params.wref;
                hkf_data["charge"] = params.charge;
                species = species.withAttachedData(hkf_data);
            }

            // Add comment/tags if present
            if(spec["Comment"])
            {
                String comment = spec["Comment"].as<String>();
                species = species.withTags({comment});
            }

            species_list.append(species);
        }
    }
    catch(const YAML::Exception& e)
    {
        errorif(true, "DEW database YAML parsing error: ", e.what());
    }
    catch(const std::exception& e)
    {
        errorif(true, "DEW database loading error: ", e.what());
    }

    return species_list;
}

/// Get the path to an embedded DEW database file
auto embeddedDatabasePath(const String& name) -> String
{
    const Strings supported = {"dew2024-aqueous", "dew2019-aqueous", "dew2024-gas", "dew2019-gas"};

    errorif(!contains(supported, name),
        "Could not find an embedded DEW database file with name `", name, "`. ",
        "The supported names are: ", join(supported, ", "), ".");

    return "databases/DEW/" + name + ".yaml";
}

} // namespace

DEWDatabase::DEWDatabase()
: Database()
{
}

DEWDatabase::DEWDatabase(const String& name)
: DEWDatabase()
{
    const auto path = embeddedDatabasePath(name);
    const auto contents = Embedded::get(path);
    const auto species_list = loadSpeciesFromYAML(contents);

    for(const auto& species : species_list)
        addSpecies(species);
}

auto DEWDatabase::load(const String& database) -> DEWDatabase&
{
    String contents;

    // Check if database is a file path or direct contents
    if(database.find('\n') != String::npos || database.find("Species:") != String::npos)
    {
        // Looks like YAML contents
        contents = database;
    }
    else
    {
        // Try to load as file
        std::ifstream file(database);
        errorif(!file.is_open(),
            "Could not open DEW database file at path `", database, "`.");

        std::stringstream buffer;
        buffer << file.rdbuf();
        contents = buffer.str();
    }

    const auto species_list = loadSpeciesFromYAML(contents);

    for(const auto& species : species_list)
        addSpecies(species);

    return *this;
}

auto DEWDatabase::withName(const String& name) -> DEWDatabase
{
    return DEWDatabase(name);
}

auto DEWDatabase::fromFile(const String& path) -> DEWDatabase
{
    DEWDatabase db;
    db.load(path);
    return db;
}

auto DEWDatabase::fromContents(const String& contents) -> DEWDatabase
{
    DEWDatabase db;
    db.load(contents);
    return db;
}

auto DEWDatabase::contents(const String& name) -> String
{
    const auto path = embeddedDatabasePath(name);
    return Embedded::get(path);
}

auto DEWDatabase::namesEmbeddedDatabases() -> Strings
{
    return {"dew2024-aqueous", "dew2019-aqueous", "dew2024-gas", "dew2019-gas"};
}

} // namespace Reaktoro
