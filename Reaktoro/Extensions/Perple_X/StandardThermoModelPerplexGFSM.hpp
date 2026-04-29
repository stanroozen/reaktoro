#pragma once

// Reaktoro includes
#include <Reaktoro/Core/StandardThermoModel.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXGFSMModel.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXPureEos.hpp>

namespace Reaktoro {

/// Parameters for the Perple_X GFSM standard thermodynamic model.
/// Provides standard state properties for individual species in GFSM fluid mixtures.
struct StandardThermoModelParamsPerplexGFSM
{
    /// Perple_X species index for the target pure fluid species (1..18)
    /// Common GFSM species: 1=H2O, 2=CO2, 3=CO, 4=CH4, 5=H2, 6=H2S, ...
    int speciesIndex = 1;

    /// The standard molar Gibbs energy of formation at reference state (J/mol)
    real G0 = 0.0;

    /// The standard molar enthalpy of formation at reference state (J/mol)
    real H0 = 0.0;

    /// The standard molar volume at reference state (m³/mol)
    real V0 = 0.0;

    /// Hybrid EOS options for H2O, CO2, CH4
    PerpleX::HybridEosOptions hybridEosOptions = PerpleX::makePerpleXHybridEosOptions();

    /// MRK mixture options
    PerpleX::MrkMixOptions mrkMixOptions{};

    /// Pure EOS options
    PerpleX::PerpleXPureEosOptions pureEosOptions{};

    /// Whether to use low-temperature MRK variant
    bool useLowTMrk = false;

    /// Maximum temperature applicable to this model (K), optional
    real Tmax = 1200.0;
};

/// Return a standard thermodynamic model for a Perple_X GFSM fluid species.
///
/// This model evaluates standard state properties (G°, H°, V°) for individual
/// species in the context of explicit GFSM fluid calculations. Properties are
/// computed using the same EOS framework as ActivityModelPerplexGFSM.
///
/// **Important**: This model evaluates properties of a *single species*, not the
/// mixture. Use with care when coupling to activity models.
///
/// @param params The parameters controlling EOS selection and calculation options
///
/// @ingroup Thermodynamics
/// @see ActivityModelPerplexGFSM
auto StandardThermoModelPerplexGFSM(
    const StandardThermoModelParamsPerplexGFSM& params) -> StandardThermoModel;

} // namespace Reaktoro
