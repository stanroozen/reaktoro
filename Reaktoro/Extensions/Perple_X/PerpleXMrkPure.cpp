#include "PerpleXMrkPure.hpp"

#include <cmath>

#include "PerpleXMrkParameters.hpp"

/// ============================================================================
/// MRK PURE-SPECIES PROPERTIES - FOUNDATION FOR GFSM
/// ============================================================================
///
/// GFSM (Generic Fluid Solution Model) operates in speciation space:
/// - User specifies all 12 mole fractions explicitly
/// - Properties are computed directly from those inputs
///
/// Three-Step GFSM Computation Process:
///
/// STEP 1 (BASELINE - This file):
/// - Compute pure-species MRK properties for all 12 species
/// - mrkPure(): Individual species MRK EOS
/// - loMrkMix(): Low-temperature MRK variant
///
/// STEP 2 (HYBRID CORRECTION - PerpleXHybridEos):
/// - For H2O: Can REPLACE with HSMRK, CORK, PSEOS, Haar, ZhangDuan05, ZhangDuan09
/// - For CO2: Can REPLACE with HSMRK, CORK, BRMRK, PSEOS, ZhangDuan09
/// - For CH4: Can REPLACE with HSMRK, ZhangDuan09
/// - For other 9: STAY with MRK (FIXED)
///
/// STEP 3 (PROPERTIES):
/// - Result is a direct function of the 12-species speciation vector
///
/// These pure-species MRK properties provide the MRK baseline used by GFSM.
/// GFSM then applies optional pure EOS substitutions for key species.

namespace Reaktoro::PerpleX {
namespace {

constexpr double rkR = 83.1441; // cm3·bar/(mol·K)

} // namespace

MrkPureResult mrkPure(const std::vector<int>& species,
                      double pressureBar,
                      double temperatureK)
{
    MrkPureResult result{};

    const double dsqrtt = std::sqrt(temperatureK);
    const double rt = rkR * temperatureK;

    const auto params = mrkParameters(temperatureK);

    for(int k = 0; k < static_cast<int>(species.size()); ++k)
    {
        const int i = species[k];
        const double aij = params.a[i];
        const double bx = params.b[i];

        const double c1 = -rt / pressureBar;
        const double c3 = -aij * bx / pressureBar / dsqrtt;
        const double c2 = c1 * bx + aij / dsqrtt / pressureBar - bx * bx;

        const auto roots = roots3(c1, c2, c3);

        double vol = 0.0;

        if(roots.iroots == 3 && roots.ineg == 0 && roots.vmin > bx)
        {
            const double v1 = roots.vmin;
            const double v2 = roots.vmax;
            const double pdv = pressureBar * (v2 - v1)
                - std::log((v2 - bx) / (v1 - bx)) * rt
                - std::log((v2 + bx) / (bx + v1) * v1 / v2) * aij / bx / dsqrtt;

            vol = (pdv > 0.0) ? v1 : v2;
        }
        else if(roots.iroots == 3)
        {
            vol = roots.vmax;
        }
        else
        {
            vol = roots.roots[roots.ipos - 1];
        }

        const double vpb = vol + bx;
        const double vmb = vol - bx;
        const double d2 = std::log(vpb / vol);

        result.v[i] = vol;
        result.ln_f[i] = bx / vmb - (1.0 / vpb + d2 / bx) * aij / rt / dsqrtt + std::log(rt / vmb);
        result.g[i] = std::exp(result.ln_f[i]) / pressureBar;

        result.vol = vol;
    }

    return result;
}

MrkMixResult loMrkMix(const std::vector<int>& species,
                      const std::array<double, 19>& y,
                      double pressureBar,
                      double temperatureK,
                      const MrkMixOptions& options)
{
    MrkMixResult result{};

    const double t2 = temperatureK * temperatureK;
    const double dsqrtt = std::sqrt(temperatureK);
    const double rt = rkR * temperatureK;

    auto params = mrkParameters(temperatureK);

    params.a[1] = 0.3930568949e9 - 0.1273025840e7 * temperatureK
        + 2049.978752 * t2 - 1.122350458 * t2 * temperatureK;
    params.a[2] = 0.9293554e8 - 0.8213073e5 * temperatureK + 0.2129e2 * t2;

    const double ch = std::exp(-11.218 + 6032.0 / temperatureK - 2782.0e3 / t2 + 4.708e8 / t2 / temperatureK)
        * 6912.824964 * t2 * dsqrtt + 79267647.0;

    std::array<double, 19> aj2{};
    double bx = 0.0;

    for(int k = 0; k < static_cast<int>(species.size()); ++k)
    {
        const int i = species[k];
        aj2[i] = 0.0;
        bx += params.b[i] * y[i];
    }

    double aij = 0.0;

    for(int k = 0; k < static_cast<int>(species.size()); ++k)
    {
        const int i = species[k];

        for(int l = 0; l < static_cast<int>(species.size()); ++l)
        {
            const int j = species[l];

            if((i == 1 && j == 2) || (i == 2 && j == 1))
            {
                aij += y[i] * y[j] * ch / 2.0;
                aj2[i] += y[j] * ch;
            }
            else
            {
                const double aij12 = y[j] * std::sqrt(params.a[i] * params.a[j]);
                aij += y[i] * aij12;
                aj2[i] += 2.0 * aij12;
            }
        }
    }

    const double c1 = -rt / pressureBar;
    const double c3 = -aij * bx / pressureBar / dsqrtt;
    const double c2 = c1 * bx + aij / dsqrtt / pressureBar - bx * bx;

    result.roots = roots3(c1, c2, c3);

    if(result.roots.iroots == 3)
        result.vol = result.roots.vmax;
    else
        result.vol = result.roots.roots[0];

    const double d2 = std::log((result.vol + bx) / result.vol);
    const double d1 = rt * dsqrtt * bx;
    const double d3 = d2 - bx / (bx + result.vol);
    const double d4 = result.vol - bx;
    const double d5 = aij * d3 / d1 / bx;
    const double d6 = std::log(rt / d4);

    for(int k = 0; k < static_cast<int>(species.size()); ++k)
    {
        const int l = species[k];

        if(y[l] <= 0.0)
        {
            result.ln_f[l] = 0.0;
            result.g[l] = 1.0;
            continue;
        }

        const double yl = (y[l] > options.minY) ? y[l] : options.minY;
        result.ln_f[l] = std::log(yl) + params.b[l] / d4 - aj2[l] / d1 * d2 + params.b[l] * d5 + d6;
        result.g[l] = std::exp(result.ln_f[l]) / pressureBar / yl;
    }

    return result;
}

} // namespace Reaktoro::PerpleX
