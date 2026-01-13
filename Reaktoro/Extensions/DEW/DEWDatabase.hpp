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
#include <Reaktoro/Core/Database.hpp>

namespace Reaktoro {

/// The class used to store and retrieve data of chemical species from DEW databases.
///
/// DEW (Deep Earth Water) databases contain thermodynamic data for aqueous species
/// and minerals at high temperatures (25-1000°C) and pressures (1-60 kbar), suitable
/// for modeling geothermal systems, metamorphic fluids, and deep crustal processes.
///
/// The DEW databases use the HKF (Helgeson-Kirkham-Flowers) model with parameters
/// calibrated for high-temperature/high-pressure conditions using the Zhang-Duan
/// water equations of state and DEW-specific electrostatic models.
///
/// @see Database
/// @ingroup Extensions
class DEWDatabase : public Database
{
public:
    /// Construct a default DEWDatabase object.
    DEWDatabase();

    /// Construct a DEWDatabase object using an embedded DEW database.
    /// If `name` does not correspond to one of the following names, an exception is thrown:
    /// - dew2024-aqueous
    /// - dew2019-aqueous
    /// - dew2024-gas
    /// - dew2019-gas
    /// @param name The name of the embedded DEW database file
    DEWDatabase(const String& name);

    /// Extend this DEWDatabase object with contents in given database file.
    /// This method supports either a path to a YAML database file, including its
    /// file name, or a multi-line string containing the database contents in YAML format.
    /// @param database The path to the database file or its contents as a string
    auto load(const String& database) -> DEWDatabase&;

    /// Return a DEWDatabase object constructed with an embedded database file.
    /// If `name` does not correspond to one of the following names, an exception is thrown:
    /// - dew2024-aqueous
    /// - dew2019-aqueous
    /// - dew2024-gas
    /// - dew2019-gas
    /// @param name The name of the embedded DEW database file
    static auto withName(const String& name) -> DEWDatabase;

    /// Return a DEWDatabase object constructed with a given local YAML file.
    /// @param path The path, including file name, to the database file.
    /// @warning An exception is thrown if `path` does not point to a valid database file.
    static auto fromFile(const String& path) -> DEWDatabase;

    /// Return a DEWDatabase object constructed with given database text contents in YAML format.
    /// @param contents The contents of the database as a YAML string.
    static auto fromContents(const String& contents) -> DEWDatabase;

    /// Return the contents of an embedded DEW database as a string.
    /// @param name The name of the embedded DEW database file
    /// @see withName
    static auto contents(const String& name) -> String;

    /// Return the names of the currently supported embedded DEW databases.
    static auto namesEmbeddedDatabases() -> Strings;
};

} // namespace Reaktoro
