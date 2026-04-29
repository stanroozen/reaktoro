#include "PerpleXGFSMModel.hpp"

#include <cmath>
#include <stdexcept>

namespace Reaktoro::PerpleX {

/// ============================================================================
/// GFSM FLUID MODEL - EXPLICIT FLUID SOLUTION MODEL IN SPECIATION SPACE
/// ============================================================================
///
/// CRITICAL DISTINCTION FROM COMPOSITION-SPACE MODELS:
/// ==========================================
///
/// Composition-space Binary Models (ifug=0-5):
/// - Solved in COMPOSITION SPACE (user specifies e.g., X_CO2)
/// - Speciation is solved internally (handled by EOS mechanism)
/// - Fixed binary structure (only one degree of freedom)
/// - All 12 species participate via built-in mixing laws
/// - Example: CO2-H2O binary at X_CO2=0.5 internally solves for all species
///
/// EXPLICIT GFSM (ifug=39) - THIS IMPLEMENTATION:
/// - Solved in SPECIATION SPACE (user specifies all 12 mole fractions directly)
/// - Speciation is EXPLICIT (user directly specifies Xn_CO2, Xn_H2O, etc.)
/// - Full flexibility (all 12 species as independent variables)
/// - Properties computed as EXPLICIT functions of speciation
/// - Example: User specifies Xn_CO2=0.3, Xn_H2O=0.6, Xn_CH4=0.1, etc. directly
///
/// This file implements GFSM computation in speciation space.
///

GFSMFluidState GFSMFluidModel::compute(
    const std::vector<int>& species,
    const std::array<double, 19>& y,
    double pressureBar,
    double temperatureK,
    const GFSMFluidOptions& options) const
{
    GFSMFluidState state{};

    // ========================================================================
    // STEP 1: MRK Mixture Foundation
    // ========================================================================
    //
    // IMPORTANT: Although we are solving in SPECIATION SPACE (EXPLICIT),
    // we begin with an MRK mixing step that combines all 12 species.
    //
    // Compute MRK properties treating all species via built-in mixing rules
    // (as if this were a composition-space model).
    //
    // Then in STEP 2: REPLACE pure H2O/CO2/CH4 properties with alternatives,
    // creating an EXPLICIT hybrid model in speciation space.

    MrkMixResult mixResult = options.useLowTMrk
        ? loMrkMix(species, y, pressureBar, temperatureK, options.mrkMixOptions)
        : mrkMix(species, y, pressureBar, temperatureK, options.mrkMixOptions);

    state.g_mrk = mixResult.g;
    state.v_mrk = mixResult.v;
    state.molarVolume = mixResult.vol;
    state.ln_f = mixResult.ln_f;

    std::array<double, 19> pureMrkG{};
    bool havePureMrk = false;

    // ========================================================================
    // STEP 2: EXPLICIT Hybrid Correction (H2O, CO2, CH4 pure EOS substitution)
    // ========================================================================
    //
    // NOW we make GFSM truly EXPLICIT in speciation space:
    //
    // For H2O, CO2, and CH4 (3 species only):
    // - REPLACE the MRK pure EOS properties with explicit alternatives:
    //   * H2O (7 options): HSMRK, CORK, PSEOS, Haar, ZhangDuan05, ZhangDuan09
    //   * CO2 (6 options): HSMRK, CORK, BRMRK, PSEOS, ZhangDuan09
    //   * CH4 (3 options): HSMRK, ZhangDuan09
    //
    // For other 9 species (H2S, SO2, H2, CO, N2, NH3, HF, C2H6, HCl):
    // - REMAIN with MRK (FIXED, no alternatives)
    //
    // Result: An EXPLICIT hybrid model that computes properties directly
    // from the user-specified mole fractions in speciation space.
    // All 12 species are handled explicitly (not solved internally).

     if(options.hybridEosOptions.water != HybridEosOptions::WaterEos::Mrk ||
         options.hybridEosOptions.co2 != HybridEosOptions::CO2Eos::Mrk ||
         options.hybridEosOptions.ch4 != HybridEosOptions::CH4Eos::Mrk)
    {
        // At least one species needs hybrid substitution
        // Compute pure MRK values first
        const auto mrkPureState = mrkPure(species, pressureBar, temperatureK);
        pureMrkG = mrkPureState.g;
        havePureMrk = true;

        // Apply hybrid pure EOS substitution
        const auto hybridResult = hybEos(
            options.hybridSpeciesIndices,
            mrkPureState.ln_f,
            mrkPureState.g,
            mrkPureState.v,
            mixResult.v,
            pressureBar,
            temperatureK,
            options.hybridEosOptions);

        state.g_hybrid = hybridResult.gh;
        state.v_hybrid = hybridResult.vhyb;
        state.hybridVolume = hybridResult.hyvol;

        // Replace fugacity coefficients for hybrid species
        applyHybridFugacity(state.g_mrk, options.hybridSpeciesIndices, hybridResult.gh);

        // Update partial molar volumes with hybrid values
        state.hybridVolume = 0.0;
        for(const int j : options.hybridSpeciesIndices)
        {
            const double yj = y[j] > 0.0 ? y[j] : options.mrkMixOptions.minY;
            state.ln_f[j] = std::log(state.g_mrk[j] * pressureBar * yj);
            state.hybridVolume += yj * state.v_hybrid[j];

            // Update partial molar volume with hybrid value
            state.v_mrk[j] = state.v_hybrid[j];
        }

        // Recompute total volume from updated partial molar volumes.
        state.molarVolume = 0.0;
        for(const int j : species)
            state.molarVolume += y[j] * state.v_mrk[j];

        // Compute volume fractions for hybrid species
        if (state.hybridVolume > 0.0)
        {
            for(const int j : options.hybridSpeciesIndices)
            {
                const double yj = y[j] > 0.0 ? y[j] : options.mrkMixOptions.minY;
                state.volumeFractions[j] = yj * state.v_hybrid[j] / state.hybridVolume;
            }
        }
    }

    // ========================================================================
    // STEP 3: Optionally enable electrolyte solvent properties
    // ========================================================================
    // If electrolyte model is enabled, compute dielectric constant and
    // solvent Gibbs contribution for aqueous solutions.

    if (options.enableElectrolyte && !species.empty())
    {
        const int nSolvent = static_cast<int>(species.size());

        // Perple_X uses vhyb(i)=v(i)+dvhy(i): for non-hybrid species dvhy=0.
        std::array<double, 19> solventVolumes = state.v_mrk;
        for(const int j : options.hybridSpeciesIndices)
            solventVolumes[j] = state.v_hybrid[j] > 0.0 ? state.v_hybrid[j] : state.v_mrk[j];

        if(!havePureMrk)
        {
            const auto mrkPureState = mrkPure(species, pressureBar, temperatureK);
            pureMrkG = mrkPureState.g;
            havePureMrk = true;
        }

        state.dielectric = computeSolventState(
            y,
            solventVolumes,
            species,
            nSolvent,
            pressureBar,
            temperatureK);

        state.solventGibbs = computeHybridSolventGibbs(
            y,
            state.g_mrk,
            pureMrkG,
            species,
            nSolvent,
            temperatureK);
    }

    return state;
}

} // namespace Reaktoro::PerpleX
