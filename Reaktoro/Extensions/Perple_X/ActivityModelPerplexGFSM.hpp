#pragma once

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXGFSMModel.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXHybridEos.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXPureEos.hpp>

namespace Reaktoro {

/// Parameters for the Perple_X GFSM (Generic Fluid Solution Model) activity model.
/// Uses explicit speciation space with optional hybrid pure EOS for H2O, CO2, CH4.
struct ActivityModelParamsPerplexGFSM
{
    /// Hybrid EOS options for H2O, CO2, CH4 pure-species substitution.
    /// Defaults to ZhangDuan09 for H2O, CO2, and CH4, matching the recommended
    /// Perple_X COH-Fluid+ configuration (iopt(25)=7, iopt(26)=7, iopt(27)=7).
    PerpleX::HybridEosOptions hybridEosOptions = PerpleX::makePerplexCOHFluidPlusEosOptions();

    /// MRK mixture calculation options
    PerpleX::MrkMixOptions mrkMixOptions{};

    /// Whether to use low-temperature MRK variant
    bool useLowTMrk = false;

    /// Whether to enable electrolyte solvent properties
    bool enableElectrolyte = false;

    /// Pure EOS options for convergence/iteration control
    PerpleX::PerpleXPureEosOptions pureEosOptions{};
};

/// Return the activity model for Perple_X GFSM (ifug=39) fluid phases.
///
/// The Generic Fluid Solution Model (GFSM) is an explicit solution model that
/// operates in SPECIATION SPACE. Users directly specify mole fractions of all
/// 13 fluid species (H2O, CO2, CH4, H2, CO, H2S, SO2, O2, N2, NH3, HF, C2H6, HCl).
///
/// Key Features:
/// - Explicit speciation space (all 12 species as independent variables)
/// - MRK mixture foundation for all species
/// - Optional hybrid pure EOS for H2O, CO2, CH4 (7, 6, 3 alternatives each)
/// - H2S, SO2, H2, CO, N2, NH3, HF, C2H6, HCl remain on MRK
/// - Dielectric/electrolyte support when enabled
///
/// Usage Example:
/// ~~~cpp
/// SpeciesList species = {H2O(g), CO2(g), CH4(g), ... }; // 13 species
/// auto actmodel = ActivityModelPerplexGFSM();
/// auto phase = Phase(species).setActivityModel(actmodel);
/// ~~~
///
/// @ingroup Thermodynamics
/// @see StandardThermoModelPerplexGFSM, PerpleX::GFSMFluidModel, PerpleXPureEos
auto ActivityModelPerplexGFSM(
    const ActivityModelParamsPerplexGFSM& params = {}) -> ActivityModelGenerator;

} // namespace Reaktoro
