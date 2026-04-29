#pragma once

#include <array>
#include <vector>

#include "PerpleXHybridEos.hpp"
#include "PerpleXMrkMixture.hpp"
#include "PerpleXMrkPure.hpp"
#include "PerpleXElectrolyte.hpp"

namespace Reaktoro::PerpleX {

/// ============================================================================
/// GFSM - Generic Fluid Solution Model (Perple_X ifug=39)
/// ============================================================================
///
/// EXPLICIT FLUID SOLUTION MODEL SOLVED IN SPECIATION SPACE
///
/// ⚠️  CRITICAL ARCHITECTURAL DISTINCTION FROM COMPOSITION-SPACE MODELS ⚠️
///
/// This is NOT a binary model. This operates in SPECIATION SPACE, not COMPOSITION SPACE.
///
/// Composition-space Binary Models (ifug=0-5):
///   - Solved in COMPOSITION SPACE (user provides X_CO2 as single parameter)
///   - Speciation is solved internally (model solves for all 12 species)
///   - Example: User sets X_CO2=0.3, model internally determines H2O, CH4, etc.
///   - Framework: All 12 species combined via MRK mixing rules
///
/// EXPLICIT GFSM (ifug=39, THIS MODEL):
///   - Solved in SPECIATION SPACE (user provides ALL 12 mole fractions explicitly)
///   - Speciation is EXTERNAL/EXPLICIT (user directly specifies every species)
///   - Example: User sets Xn_CO2=0.3, Xn_H2O=0.6, Xn_CH4=0.1, etc. (ALL 12)
///   - Framework: MRK mixture foundation + selective pure EOS replacements
///
/// Key Characteristics of EXPLICIT GFSM:
/// - Input: Temperature, Pressure, ALL 12 species mole fractions (user provides)
/// - Speciation: NO INTERNAL SOLVING (speciation completely specified by user)
/// - Properties: Direct function of input mole fractions (no internal equilibrium solver)
/// - Flexibility: Can use any subset of 12 species or all 12
/// - Pure EOS: User selects alternatives for H2O, CO2, CH4; others fixed to MRK
///
/// Three-Step Algorithm:
/// 1. MRK Foundation: Combine all 12 species via MRK mixing rules
/// 2. EXPLICIT Hybrid Correction: For H2O/CO2/CH4, replace pure EOS if needed
/// 3. Return EXPLICIT Properties: Direct function of user-specified speciation
///
/// Pure EOS Options Available:
/// - H2O (7 total): MRK, HSMRK, CORK, PSEOS, Haar, ZhangDuan05, ZhangDuan09
/// - CO2 (6 total): MRK, HSMRK, CORK, BRMRK, PSEOS, ZhangDuan09
/// - CH4 (3 total): MRK, HSMRK, ZhangDuan09
/// - Others (9 fixed): H2S, SO2, H2, CO, N2, NH3, HF, C2H6, HCl → ALWAYS MRK
///
/// DO NOT CONFUSE:
/// - Binary (composition space) ≠ GFSM (speciation space)
/// - User provides composition (X_CO2) → binary model solves speciation internally
/// - User provides speciation (Xn_all) → GFSM computes properties explicitly
/// - These are fundamentally different model architectures with different solution strategies
///
/// References:
/// - Perple_X solution_model.dat lines 11920-12030 (COH-Fluid definition)
/// - Perple_X rkparm() subroutine: MRK parameters for all 12 species
/// - Perple_X hybeos() subroutine: Pure EOS selection logic
///

struct GFSMFluidOptions
{
    /// Pure EOS selection for hybrid framework (GFSM only)
    HybridEosOptions hybridEosOptions{};

    /// MRK mixture options (applies to all 12 species initially)
    MrkMixOptions mrkMixOptions{};

    /// Species indices that use hybrid pure EOS (typically H2O, CO2, CH4)
    /// For all other species, MRK pure EOS is used (fixed)
    std::vector<int> hybridSpeciesIndices;

    /// Enable electrolyte solvent properties (for aqueous solutions)
    bool enableElectrolyte = false;

    /// Use low-temperature MRK variant (for very low T)
    bool useLowTMrk = false;
};

struct GFSMFluidState
{
    /// Fugacity coefficients for all 12 species (from MRK mixture)
    std::array<double, 19> g_mrk{};

    /// Partial molar volumes for all 12 species (from MRK mixture)
    std::array<double, 19> v_mrk{};

    /// Log fugacities for all 12 species
    std::array<double, 19> ln_f{};

    /// Hybrid-corrected Gibbs contributions (only for hybrid species)
    std::array<double, 19> g_hybrid{};

    /// Hybrid partial molar volumes (only for hybrid species)
    std::array<double, 19> v_hybrid{};

    /// Total molar volume (MRK mixture basis)
    double molarVolume = 0.0;

    /// Volume fraction sum for hybrid species
    double hybridVolume = 0.0;

    /// Volume fractions of hybrid species
    std::array<double, 19> volumeFractions{};

    /// Dielectric constant (when electrolyte enabled)
    DielectricState dielectric{};

    /// Hybrid solvent Gibbs contribution (when electrolyte enabled)
    double solventGibbs = 0.0;

    /// Root-finding state for continuation (for future Newton-Raphson iterations)
    MrkRootState rootState{};
};

/// ============================================================================
/// GFSM Fluid Model Implementation
/// ============================================================================
///
/// Computes thermodynamic properties using the GFSM (Generic Fluid Solution
/// Model, Perple_X ifug=39) in the speciation space.
///
/// This is an EXPLICIT fluid solution model that computes properties from
/// direct mole fractions of up to 12 fluid species.
///
class GFSMFluidModel
{
public:
    /// Compute GFSM fluid properties
    ///
    /// Arguments:
    ///   - species: indices of active fluid species (e.g., [1,2] for H2O+CO2)
    ///   - y: mole fractions of all 19 potential species (only indices in species are used)
    ///   - pressureBar: fluid pressure in bar
    ///   - temperatureK: absolute temperature in Kelvin
    ///   - options: GFSM configuration (hybrid EOS choices, etc.)
    ///
    /// Returns:
    ///   GFSMFluidState with all thermodynamic properties computed
    ///
    GFSMFluidState compute(
        const std::vector<int>& species,
        const std::array<double, 19>& y,
        double pressureBar,
        double temperatureK,
        const GFSMFluidOptions& options) const;
};

} // namespace Reaktoro::PerpleX
