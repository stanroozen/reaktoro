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
#include <Reaktoro/Core/StandardThermoModel.hpp>
#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>

namespace Reaktoro {

/// The parameters in the DEW model for calculating standard thermodynamic properties of aqueous solutes.
/// Uses the same HKF parameters but with DEW water models for computing Born functions and water properties.
struct StandardThermoModelParamsDEW
{
    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in J/mol).
    real Gf;

    /// The apparent standard molal enthalpy of formation of the species from its elements (in J/mol).
    real Hf;

    /// The standard molal entropy of the species at reference temperature and pressure (in J/(mol·K)).
    real Sr;

    /// The coefficient `a1` of the HKF equation of state of the aqueous solute (in J/(mol·Pa)).
    real a1;

    /// The coefficient `a2` of the HKF equation of state of the aqueous solute (in J/mol).
    real a2;

    /// The coefficient `a3` of the HKF equation of state of the aqueous solute (in (J·K)/(mol·Pa)).
    real a3;

    /// The coefficient `a4` of the HKF equation of state of the aqueous solute (in (J·K)/mol).
    real a4;

    /// The coefficient `c1` of the HKF equation of state of the aqueous solute (in J/(mol·K)).
    real c1;

    /// The coefficient `c2` of the HKF equation of state of the aqueous solute (in (J·K)/mol).
    real c2;

    /// The conventional Born coefficient of the aqueous solute at reference temperature 298.15 K and pressure 1 bar (in J/mol).
    real wref;

    /// The electrical charge of the aqueous solute.
    real charge;

    /// The maximum temperature at which the DEW model can be applied for the substance (optional, in K).
    real Tmax;

    /// The water model options to use for DEW calculations (EOS, dielectric, Born, Gibbs models).
    /// Defaults to DEW preset (ZhangDuan2005, PowerFunction, DewIntegral, Shock92Dew).
    WaterModelOptions waterOptions = makeWaterModelOptionsDEW();
};

/// Return a function that calculates thermodynamic properties of an aqueous solute using the DEW model.
/// This model uses the same HKF equation structure but with DEW-specific water models for:
/// - Density and derivatives (ZhangDuan EOS)
/// - Dielectric constant (PowerFunction or JohnsonNorton)
/// - Born solvation functions (Shock92)
/// - Water Gibbs energy (∫V dP integral)
auto StandardThermoModelDEW(const StandardThermoModelParamsDEW& params) -> StandardThermoModel;

} // namespace Reaktoro
