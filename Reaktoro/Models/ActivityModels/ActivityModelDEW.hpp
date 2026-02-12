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
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

/// Return the activity model for aqueous electrolyte phases based on the DEW (Deep Earth Water) model.
///
/// This model implements the Helgeson–Kirkham–Flowers (HKF) extended Debye–Hückel formulation
/// with electrostatic coefficients (A and B parameters) evaluated dynamically from water density
/// and dielectric constant using the DEW equation of state, valid up to 1000°C and 60 kbar.
///
/// **References:**
///   - Helgeson, H. C., Kirkham, D. H., Flowers, G. C. (1981). Theoretical
///     prediction of the thermodynamic behavior of aqueous electrolytes at
///     high pressures and temperatures: IV. Calculation of activity
///     coefficients, osmotic coefficients, and apparent molal and standard and
///     relative partial molal properties to 600°C. American Journal of
///     Science, 281(10), 1249–1516.
///   - Huang, F., Sverjensky, D. A. (2019). Extended Deep Earth Water (DEW)
///     Model for the depths of the Earth's crust and upper mantle.
///     Geochimica et Cosmochimica Acta, 3, 260, 149–161.
///
/// **Key Features:**
///   - Dynamic A(T,P) and B(T,P) from water density and dielectric constant
///   - Extended-term parameter b_c,k set to zero at deep-Earth conditions
///   - Explicit finite-water correction term
///   - Valid from ambient to extreme conditions (0.1 to 1000°C, 1 bar to 60 kbar)
///
/// @ingroup Thermodynamics
auto ActivityModelDEW() -> ActivityModelGenerator;

} // namespace Reaktoro
