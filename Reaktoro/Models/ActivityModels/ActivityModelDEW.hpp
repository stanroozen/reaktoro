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
#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>

/// ============================================================
/// ActivityModelDEW — capability matrix
/// ============================================================
///
/// Supported workflows
///   ✔  Pure aqueous chemistry with DEW database
///   ✔  T/P constraint: temperature() + pressure()
///   ✔  pH constraint: specs.pH()
///   ✔  Fugacity constraint: specs.fugacity()
///   ✔  EquilibriumSensitivity (dndw)
///   ✔  KineticsSolver (model-agnostic kinetics path)
///   ✔  Operator-splitting reactive transport (Python level)
///   ✔  Davies Debye-Hückel (dhModel = Davies)
///   ✔  Extended Debye-Hückel with ion-size a (dhModel = ExtendedDH, default)
///   ✔  Configurable water submodels (EOS, dielectric, Gibbs, Born)
///   ✔  Configurable Psat handling and density tolerance
///   ✔  Extended DH b-dot term (bExtended)
///   ✔  Exports AqueousMixtureState into props.extra
///
/// Not supported in this model
///   ✘  Perple_X GFSM gas-phase water-activity handoff
///   ✘  GFSM standard-state conflict detection
///   ✘  Perple_X hybrid EOS selection for gas phase
///
/// Constructor signatures (C++ and Python)
///   ActivityModelDEW()
///   ActivityModelDEW(ActivityModelParamsDEW const& params)
///
/// Debye-Hückel default
///   dhModel = ActivityDHModel::ExtendedDH
///
/// See also: ActivityModelPerplexDEW.hpp for the Perple_X-backed variant
///           that shares the same ActivityDHModel enum and dhModel field.
/// ============================================================

namespace Reaktoro {

/// Debye-Hückel variant used by the DEW and PerplexDEW aqueous activity models.
///
/// Shared by ActivityModelParamsDEW and ActivityModelParamsPerplexDEW so the
/// same field name (`dhModel`) selects the variant on both backends.
enum class ActivityDHModel
{
    /// Extended Debye-Hückel with species-specific ionic radii (aᵢ):
    ///   log₁₀(γᵢ) = −A zᵢ² √I / (1 + aᵢ B √I) + b I + Cc
    /// A and B are computed from the current water density and dielectric constant.
    /// This is the default for both DEW and PerplexDEW.
    ExtendedDH,

    /// Classic Davies approximation (no ionic-radius parameter):
    ///   log₁₀(γᵢ) = −A zᵢ² (√I / (1 + √I) − 0.3 I)
    /// A is computed from the current water density and dielectric constant.
    /// Valid to ~0.5 mol/kg ionic strength.
    Davies
};

/// Options for the aqueous DEW activity model.
struct ActivityModelParamsDEW
{
    /// Debye-Hückel variant. Default: ExtendedDH (full HKF ionic-radius form).
    /// Set to ActivityDHModel::Davies to use the ionic-radius-free Davies equation.
    ActivityDHModel dhModel = ActivityDHModel::ExtendedDH;

	/// Water-property submodels used to evaluate A(T,P) and B(T,P).
	WaterModelOptions waterOptions = makeWaterModelOptionsDEW();

	/// Extended Debye-Huckel correction term b_c,k.
	///
	/// The DEW default is zero at deep-Earth conditions.
	real bExtended = 0.0;
};

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

/// Return the activity model for aqueous electrolyte phases based on the DEW (Deep Earth Water) model.
///
/// This overload permits callers to configure the underlying DEW water-property models.
auto ActivityModelDEW(const ActivityModelParamsDEW& params) -> ActivityModelGenerator;

} // namespace Reaktoro
