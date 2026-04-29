#include "PerpleXHybridEos.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace Reaktoro::PerpleX {

/// ============================================================================
/// GFSM HYBRID PURE-EOS SELECTION FRAMEWORK (EXPLICIT SPECIATION SPACE)
/// ============================================================================
///
/// NOTE: This is PURE-SPECIES EOS selection for GFSM (Type 39).
///
/// EXPLICIT GFSM (this file):
/// - Composition space formulation: User sets Xn_CO2, Xn_H2O, ... Xn_HCl (all 12)
/// - Speciation EXPLICIT: All 12 mole fractions specified directly by user
/// - Pure EOS: Each species gets pure-species EOS evaluation (independent)
/// - Properties computed directly from these 12 pure EOS evaluations
/// - Flexible: Can use any subset of 12 species or all 12
/// - Pure EOS for H2O, CO2, CH4 can be MRK, HSMRK, CORK, PSEOS, Haar, ZD05, ZD09
/// - Other 9 species: Fixed to MRK
///
/// This file implements PURE-SPECIES EOS SELECTION for EXPLICIT GFSM:
/// - H2O pure species: Can be MRK, HSMRK, CORK, PSEOS, Haar, ZhangDuan05, ZhangDuan09
/// - CO2 pure species: Can be MRK, HSMRK, CORK, BRMRK, PSEOS, ZhangDuan09
/// - CH4 pure species: Can be MRK, HSMRK, ZhangDuan09
/// - Other 9 species: FIXED to MRK
///
/// Result: Properties in EXPLICIT speciation space (direct function of 12 mole fractions)
///

namespace {

/// ============================================================================
/// GFSM Hybrid Pure EOS Evaluation
/// ============================================================================
///
/// These helper functions evaluate alternative pure EOS formulations as part
/// of the GFSM framework. They replace MRK values for H2O, CO2, and CH4.
///

inline double evalPureEos(const HybridEosOptions::PureEosCallback& callback,
                          int species,
                          double& volume,
                          double pressureBar,
                          double temperatureK,
                          const char* name)
{
    if(callback)
        return callback(species, volume, pressureBar, temperatureK);

    throw std::runtime_error(std::string("Hybrid EoS callback not provided for ") + name);
}

} // namespace

/// ============================================================================
/// hybEos - GFSM Pure-Species EOS Correction Framework
/// ============================================================================
///
/// This function implements the Perple_X hybeos() subroutine, which applies
/// the GFSM framework to replace MRK pure-species properties with alternative
/// pure EOS choices for H2O, CO2, and CH4.
///
/// GFSM Architecture:
/// 1. User specifies all 12 mole fractions explicitly
/// 2. For H2O, CO2, CH4: optionally replace MRK with HSMRK, CORK, PSEOS, Haar, etc.
/// 3. Other 9 species: remain as MRK (fixed)
/// 4. Combine all species in explicit speciation space
///
/// This is an EXPLICIT fluid solution model because:
/// - Properties are computed directly from mole fractions
/// - Each species contributes independently to final properties
/// - No hidden solving for composition is needed
/// - No hardcoded H2O-CO2 coupling
///
HybridEosResult hybEos(const std::vector<int>& species,
                       const std::array<double, 19>& mrk_ln_f,
                       const std::array<double, 19>& mrk_g,
                       const std::array<double, 19>& mrk_v,
                       const std::array<double, 19>& mrk_mix_v,
                       double pressureBar,
                       double temperatureK,
                       const HybridEosOptions& options)
{
    HybridEosResult result{};

    // Initialize with MRK values (default for all species, and final for 9 species)
    result.ln_f = mrk_ln_f;
    result.g = mrk_g;
    result.v = mrk_v;

    // Store original MRK values (Perple_X common blocks vmrk0, gmrk0)
    result.vmrk0 = mrk_v;
    result.gmrk0 = mrk_g;

    for(int i = 0; i < static_cast<int>(species.size()); ++i)

    {
        const int j = species[i];

        result.vh[j] = -mrk_v[j];
        result.gh[j] = mrk_g[j];

        double f = mrk_ln_f[j];
        double v = mrk_v[j];

        // ====================================================================
        // H2O PURE EOS SELECTION (species index 1)
        // ====================================================================
        // Default: MRK
        // Alternatives: HSMRK, CORK, PSEOS, Haar, ZhangDuan05, ZhangDuan09
        //
        if(j == 1)
        {
            switch(options.water)
            {
            case HybridEosOptions::WaterEos::Mrk:
                // MRK is default, no change needed
                break;
            case HybridEosOptions::WaterEos::Hsmrk:
                f = evalPureEos(options.hsmrk, j, v, pressureBar, temperatureK, "HSMRK");
                break;
            case HybridEosOptions::WaterEos::Cork:
                f = evalPureEos(options.cork, j, v, pressureBar, temperatureK, "CORK");
                break;
            case HybridEosOptions::WaterEos::Pseos:
                f = evalPureEos(options.pseos, j, v, pressureBar, temperatureK, "PSEOS");
                break;
            case HybridEosOptions::WaterEos::Haar:
                f = evalPureEos(options.haar, j, v, pressureBar, temperatureK, "Haar");
                break;
            case HybridEosOptions::WaterEos::ZhangDuan05:
                f = evalPureEos(options.zhangDuan05, j, v, pressureBar, temperatureK, "ZhangDuan05");
                break;
            case HybridEosOptions::WaterEos::ZhangDuan09:
                f = evalPureEos(options.zhangDuan09, j, v, pressureBar, temperatureK, "ZhangDuan09");
                break;
            }

            result.ln_f[j] = f;
            result.v[j] = v;
            result.g[j] = std::exp(f) / pressureBar;
            result.vol = v;
        }
        // ====================================================================
        // CO2 PURE EOS SELECTION (species index 2)
        // ====================================================================
        // Default: MRK
        // Alternatives: HSMRK, CORK, BRMRK, PSEOS, ZhangDuan09
        //
        else if(j == 2)
        {
            switch(options.co2)
            {
            case HybridEosOptions::CO2Eos::Mrk:
                // MRK is default, no change needed
                break;
            case HybridEosOptions::CO2Eos::Hsmrk:
                f = evalPureEos(options.hsmrk, j, v, pressureBar, temperatureK, "HSMRK");
                break;
            case HybridEosOptions::CO2Eos::Cork:
                f = evalPureEos(options.cork, j, v, pressureBar, temperatureK, "CORK");
                break;
            case HybridEosOptions::CO2Eos::Brmrk:
                f = evalPureEos(options.brmrk, j, v, pressureBar, temperatureK, "BRMRK");
                break;
            case HybridEosOptions::CO2Eos::Pseos:
                f = evalPureEos(options.pseos, j, v, pressureBar, temperatureK, "PSEOS");
                break;
            case HybridEosOptions::CO2Eos::ZhangDuan09:
                f = evalPureEos(options.zhangDuan09, j, v, pressureBar, temperatureK, "ZhangDuan09");
                break;
            }

            result.ln_f[j] = f;
            result.v[j] = v;
            result.g[j] = std::exp(f) / pressureBar;
            result.vol = v;
        }
        else if(j == 4)
        {
            switch(options.ch4)
            {
            case HybridEosOptions::CH4Eos::Mrk:
                break;
            case HybridEosOptions::CH4Eos::Hsmrk:
                f = evalPureEos(options.hsmrk, j, v, pressureBar, temperatureK, "HSMRK");
                break;
            case HybridEosOptions::CH4Eos::ZhangDuan09:
                f = evalPureEos(options.zhangDuan09, j, v, pressureBar, temperatureK, "ZhangDuan09");
                break;
            }

            result.ln_f[j] = f;
            result.v[j] = v;
            result.g[j] = std::exp(f) / pressureBar;
        }

        result.vh[j] += result.v[j];
        result.gh[j] = result.g[j] / result.gh[j];
    }

    // Compute hybrid volumes and volume fractions (rlib.f lines 11692-11707)
    // vhyb(i) = dvhy(i) + v(i), where dvhy(i) = v_pure_hybrid(i) - v_pure_mrk(i)
    // and v(i) here is the MRK mixture partial molar volume.
    result.hyvol = 0.0;
    for(int i = 0; i < static_cast<int>(species.size()); ++i)
    {
        const int j = species[i];
        result.vhyb[j] = result.vh[j] + mrk_mix_v[j];
        // hyvol = sum(yf[i] * vhyb[i]) - need mole fractions from caller
        // For now just store volumes, fractions computed in fluid model
    }

    return result;
}

} // namespace Reaktoro::PerpleX
