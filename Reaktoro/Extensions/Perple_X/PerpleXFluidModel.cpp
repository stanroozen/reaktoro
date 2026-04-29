#include "PerpleXFluidModel.hpp"

#include <cmath>
#include <stdexcept>

namespace Reaktoro::PerpleX {

PerpleXFluidState PerpleXFluidModel::compute(const std::vector<int>& species,
                                             const std::array<double, 19>& y,
                                             double pressureBar,
                                             double temperatureK,
                                             const PerpleXFluidOptions& options,
                                             MrkRootState* rootState) const
{
    PerpleXFluidState state{};

    MrkMixResult mix = options.useLowTMrk
        ? loMrkMix(species, y, pressureBar, temperatureK, options.mixOptions)
        : mrkMix(species, y, pressureBar, temperatureK, options.mixOptions, rootState);

    state.ln_f = mix.ln_f;
    state.g = mix.g;
    state.v = mix.v;
    state.vol = mix.vol;

    if(!options.hybridSpecies.empty())
    {
        for(const int j : options.hybridSpecies)
        {
            if(j != 1 && j != 2 && j != 4)
            {
                throw std::invalid_argument("PerpleXFluidModel hybridSpecies is limited to GFSM-callable species {1(H2O),2(CO2),4(CH4)}");
            }
        }

        const auto mrkPureState = mrkPure(species, pressureBar, temperatureK);

        const auto hybrid = hybEos(options.hybridSpecies,
                                   mrkPureState.ln_f,
                                   mrkPureState.g,
                                   mrkPureState.v,
                                   mix.v,
                                   pressureBar,
                                   temperatureK,
                                   options.hybridOptions);

        state.gh = hybrid.gh;
        state.vh = hybrid.vh;
        state.vhyb = hybrid.vhyb;
        state.hyvol = hybrid.hyvol;

        applyHybridFugacity(state.g, options.hybridSpecies, hybrid.gh);

        // Compute hybrid volumes with composition
        state.hyvol = 0.0;
        for(const int j : options.hybridSpecies)
        {
            const double yj = y[j] > 0.0 ? y[j] : options.mixOptions.minY;
            state.ln_f[j] = std::log(state.g[j] * pressureBar * yj);
            state.hyvol += yj * state.vhyb[j];

            // Update partial molar volumes with hybrid values
            state.v[j] = state.vhyb[j];
        }

        // Recompute total volume from updated partial molar volumes.
        state.vol = 0.0;
        for(const int j : species)
            state.vol += y[j] * state.v[j];

        // Compute volume fractions
        if (state.hyvol > 0.0)
        {
            for(const int j : options.hybridSpecies)
            {
                const double yj = y[j] > 0.0 ? y[j] : options.mixOptions.minY;
                state.vf[j] = yj * state.vhyb[j] / state.hyvol;
            }
        }

        // Compute electrolyte properties if enabled
        if (options.enableElectrolyte)
        {
            const int nSolvent = static_cast<int>(species.size());

            // Perple_X uses vhyb(i)=v(i)+dvhy(i): for non-hybrid species dvhy=0.
            std::array<double, 19> solventVolumes = state.v;
            for(const int j : options.hybridSpecies)
                solventVolumes[j] = state.vhyb[j] > 0.0 ? state.vhyb[j] : state.v[j];

            state.dielectric = computeSolventState(y, solventVolumes, species,
                                                   nSolvent, pressureBar, temperatureK);

            state.gsolv = computeHybridSolventGibbs(y, state.g, mrkPureState.g, species, nSolvent, temperatureK);
        }
    }

    return state;
}

} // namespace Reaktoro::PerpleX
