#pragma once

#include <array>
#include <functional>
#include <vector>

namespace Reaktoro::PerpleX {

/// ============================================================================
/// GFSM PURE-SPECIES EOS SELECTION FRAMEWORK (EXPLICIT SPECIATION SPACE)
/// ============================================================================
///
/// This header defines the pure EOS selection mechanism used by the GFSM
/// (Generic Fluid Solution Model, Perple_X Type 39 solution model).
///
/// GFSM Architecture (Explicit Speciation Space):
/// -----------------------------------------------
/// - User specifies: Mole fractions of all 12 species (Xn_CO2, Xn_H2O, ... Xn_HCl)
/// - GFSM applies pure species EOS to each species independently
/// - H2O, CO2, CH4: Can use alternative EOS (HSMRK, CORK, PSEOS, Haar, etc.)
/// - Other 9 species: Fixed to MRK (no switching)
/// - Result: Direct property calculation from explicit speciation
///
/// Pure EOS Options by Species:
////
/// H2O (species index 1):
///   - MRK (default): Modified Redlich-Kwong (fast, good accuracy)
///   - HSMRK: Kerrick & Jacobs hard-sphere MRK (corrects volume)
///   - CORK: Holland & Powell CORK model (high P/T)
///   - PSEOS: Pseff pseudo-Einstein model
///   - Haar: Haar et al. 1979 (very accurate for water)
///   - ZhangDuan05: Zhang & Duan 2005 correlation
///   - ZhangDuan09: Zhang & Duan 2009 correlation
///
/// CO2 (species index 2):
///   - MRK (default): Modified Redlich-Kwong
///   - HSMRK: Hard-sphere MRK
///   - CORK: Holland & Powell CORK
///   - BRMRK: Bottinga & Richet MRK variant (CO2-specific)
///   - PSEOS: Pseff model
///   - ZhangDuan09: Zhang & Duan 2009
///
/// CH4 (species index 3):
///   - MRK (default): Modified Redlich-Kwong
///   - HSMRK: Hard-sphere MRK
///   - ZhangDuan09: Zhang & Duan 2009
///
/// All Other 9 Species (H2S, SO2, H2, CO, N2, NH3, HF, C2H6, HCl):
///   - FIXED to MRK (no switching)
///
/// NOTE: This framework is GFSM-only (explicit speciation space).

struct HybridEosOptions
{
    /// Pure EOS options for H2O (7 available)
    enum class WaterEos {
        Mrk,           ///< Modified Redlich-Kwong (default, fast)
        Hsmrk,         ///< Kerrick & Jacobs hard-sphere MRK
        Cork,          ///< Holland & Powell CORK (high P/T)
        Pseos,         ///< Pseff pseudo-Einstein
        Haar,          ///< Haar et al. 1979 (most accurate)
        ZhangDuan05,   ///< Zhang & Duan 2005 correlation
        ZhangDuan09    ///< Zhang & Duan 2009 correlation
    };

    /// Pure EOS options for CO2 (6 available)
    enum class CO2Eos {
        Mrk,           ///< Modified Redlich-Kwong (default)
        Hsmrk,         ///< Hard-sphere MRK
        Cork,          ///< Holland & Powell CORK
        Brmrk,         ///< Bottinga & Richet MRK variant
        Pseos,         ///< Pseff model
        ZhangDuan09    ///< Zhang & Duan 2009
    };

    /// Pure EOS options for CH4 (3 available)
    enum class CH4Eos {
        Mrk,           ///< Modified Redlich-Kwong (default)
        Hsmrk,         ///< Hard-sphere MRK
        ZhangDuan09    ///< Zhang & Duan 2009
    };

    using PureEosCallback = std::function<double(int /*species*/, double& /*volume*/, double /*pressureBar*/, double /*temperatureK*/)>;

    WaterEos water = WaterEos::Mrk;     ///< Default: MRK for H2O
    CO2Eos co2 = CO2Eos::Mrk;           ///< Default: MRK for CO2
    CH4Eos ch4 = CH4Eos::Mrk;           ///< Default: MRK for CH4

    /// Callback functions for alternative pure EOS evaluations
    /// These are set by caller to provide implementations of HSMRK, CORK, etc.
    PureEosCallback hsmrk;      ///< HSMRK pure EOS implementation
    PureEosCallback cork;       ///< CORK pure EOS implementation
    PureEosCallback brmrk;      ///< BRMRK pure EOS implementation
    PureEosCallback pseos;      ///< PSEOS pure EOS implementation
    PureEosCallback haar;       ///< Haar pure EOS implementation
    PureEosCallback zhangDuan05; ///< Zhang & Duan 2005 implementation
    PureEosCallback zhangDuan09; ///< Zhang & Duan 2009 implementation
};

struct HybridEosResult
{
    /// Pure MRK fugacity coefficients (Perple_X ln_f array)
    std::array<double, 19> ln_f{};

    /// Pure MRK Gibbs contributions (Perple_X g array)
    std::array<double, 19> g{};

    /// Pure MRK partial molar volumes (Perple_X v array)
    std::array<double, 19> v{};

    /// Hybrid Gibbs contributions (for H2O, CO2, CH4 only)
    std::array<double, 19> gh{};

    /// Hybrid partial molar volumes (for H2O, CO2, CH4 only)
    std::array<double, 19> vh{};

    /// MRK volumes before hybrid substitution (backup from Perple_X vmrk0)
    std::array<double, 19> vmrk0{};

    /// Hybrid partial molar volumes (corrected Perple_X vhyb)
    std::array<double, 19> vhyb{};

    /// Volume fractions for hybrid species (Perple_X vf)
    std::array<double, 19> vf{};

    /// MRK Gibbs before hybrid substitution (backup from Perple_X gmrk0)
    std::array<double, 19> gmrk0{};

    /// Total molar volume (MRK mixture basis)
    double vol = 0.0;

    /// Total hybrid volume for selected species (sum of yf[i]*vhyb[i])
    double hyvol = 0.0;
};

/// ============================================================================
/// Apply GFSM hybrid pure-species EoS corrections
/// ============================================================================
///
/// This function implements the GFSM framework by replacing MRK pure-species
/// properties with alternative pure EOS for H2O, CO2, and CH4.
///
/// Purpose:
/// - Takes MRK pure-species results (universal default)
/// - For selected species (H2O, CO2, CH4), substitutes alternative pure EOS
/// - Other 9 species remain as MRK (no substitution)
/// - Returns hybrid-corrected properties for explicit speciation space calculation
///
/// Corresponds to Perple_X hybeos() subroutine (flib.f lines ~7950-8100)
///
/// Arguments:
///   - species: indices of active species
///   - mrk_ln_f: pure MRK log fugacity coefficients (Perple_X ln_f)
///   - mrk_g: pure MRK Gibbs coefficients (Perple_X g)
///   - mrk_v: pure MRK partial molar volumes (Perple_X v)
///   - pressureBar: fluid pressure
///   - temperatureK: absolute temperature
///   - options: GFSM hybrid pure EOS selections
///
/// Returns:
///   HybridEosResult with hybrid-corrected properties for GFSM
///
HybridEosResult hybEos(const std::vector<int>& species,
                       const std::array<double, 19>& mrk_ln_f,
                       const std::array<double, 19>& mrk_g,
                       const std::array<double, 19>& mrk_v,
                       const std::array<double, 19>& mrk_mix_v,
                       double pressureBar,
                       double temperatureK,
                       const HybridEosOptions& options);

} // namespace Reaktoro::PerpleX
