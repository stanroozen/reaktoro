#include "PerpleXMrkMixture.hpp"

#include <cmath>

namespace Reaktoro::PerpleX {

/// ============================================================================
/// MRK MIXTURE IMPLEMENTATION
/// ============================================================================
///
/// This file implements MRK mixture mixing rules used as the baseline for
/// GFSM calculations. GFSM uses these results and then applies optional
/// pure-species substitutions for H2O/CO2/CH4.
///
/// This implementation focuses on mixing-rule behavior in (P, T) space and
/// provides the initial mixture properties used by the GFSM workflow.
///
namespace {

constexpr double rkR = 83.1441; // cm3·bar/(mol·K)

inline double harmonicMean(double a, double b)
{
    return 2.0 / (1.0 / a + 1.0 / b);
}

} // namespace

/// ============================================================================
/// roots3 - Cubic Equation Solver
/// ============================================================================
///
/// Solves the cubic equation: x³ + a1*x² + a2*x + a3 = 0
///
/// Used in MRK mixture calculations to find volume roots from the
/// cubic EOS formulation. This corresponds to Perple_X roots3()
/// subroutine (flib.f).
///
/// Arguments:
///   - a1, a2, a3: cubic equation coefficients
///
/// Returns:
///   Roots3Result with up to 3 real roots and diagnostic info
///
Roots3Result roots3(double a1, double a2, double a3)
{
    Roots3Result res{};

    const double qq = (a1 * a1 - 3.0 * a2) / 9.0;
    const double rr = (a1 * (2.0 * a1 * a1 - 9.0 * a2) + 27.0 * a3) / 54.0;
    const double a5 = a1 / 3.0;

    const double dif = qq * qq * qq - rr * rr;

    if(dif >= 0.0)
    {
        const double phi = (dif > 0.0) ? std::acos(rr / std::pow(qq, 1.5)) : 0.0;
        const double a4 = -2.0 * std::sqrt(qq);
        const double a6 = phi / 3.0;

        double dphi = 0.0;
        res.vmin = 1.0e9;
        res.vmax = -1.0e9;
        res.ineg = 0;

        for(int i = 0; i < 3; ++i)
        {
            const double v = a4 * std::cos(a6 + dphi) - a5;
            if(v > res.vmax) res.vmax = v;
            if(v < res.vmin) res.vmin = v;
            if(v <= 0.0)
                res.ineg += 1;
            else
                res.ipos = i + 1;

            res.roots[i] = v;
            dphi += 2.094395102497915;
        }

        res.iroots = 3;
    }
    else
    {
        const double a7 = std::pow(std::sqrt(-dif) + std::abs(rr), 1.0 / 3.0);
        res.roots[0] = -rr / std::abs(rr) * (a7 + qq / a7) - a5;
        res.iroots = 1;
        res.ineg = 0;
        res.ipos = 1;
    }

    return res;
}

/// ============================================================================
/// mrkMix - MRK Mixture Model (Perple_X mrkmix)
/// ============================================================================
///
/// Computes MRK mixture fugacity coefficients and partial molar volumes
/// for all species using MRK mixing rules.
///
/// Framework:
/// - Uses built-in mixing rules (geometric mean for cross-coefficients)
/// - All species contributions combined via mixing law
/// - Cubic volume equation solved for each P-T point
/// - Returns fugacities and volumes for the MRK mixture
///
/// This function provides the BASIS for GFSM:
/// - GFSM uses mrkMix() to get initial fugacities/volumes for all 12 species
/// - Then optionally replaces H2O, CO2, CH4 pure EOS with hybrid alternatives
/// - Result: explicit speciation-space model (GFSM) built from MRK foundation
///
/// Corresponds to Perple_X mrkmix() subroutine (flib.f lines ~5500-6000)
///
MrkMixResult mrkMix(const std::vector<int>& species,
                    const std::array<double, 19>& y,
                    double pressureBar,
                    double temperatureK,
                    const MrkMixOptions& options,
                    MrkRootState* rootState)
{
    MrkMixResult result{};

    const double dsqrtt = std::sqrt(temperatureK);
    const double rt = rkR * temperatureK;

    const auto params = mrkParameters(temperatureK);

    std::array<double, 19> aj2{};
    double bx = 0.0;

    for(int k = 0; k < static_cast<int>(species.size()); ++k)
    {
        const int i = species[k];
        const double yi = (y[i] < 0.0) ? 0.0 : y[i];
        aj2[i] = 0.0;
        bx += params.b[i] * yi;
    }

    const double ch = std::exp(-11.218 + (6032.0 + (-2782000.0 + 4.708e8 / temperatureK) / temperatureK) / temperatureK)
        * 6912.824964 * temperatureK * temperatureK * dsqrtt + 79267647.0;

    double aij = 0.0;

    for(int k = 0; k < static_cast<int>(species.size()); ++k)
    {
        const int i = species[k];
        const double yi = (y[i] < 0.0) ? 0.0 : y[i];

        for(int l = 0; l < static_cast<int>(species.size()); ++l)
        {
            const int j = species[l];
            const double yj = (y[j] < 0.0) ? 0.0 : y[j];

            if((i == 1 && j == 2) || (i == 2 && j == 1))
            {
                aij += yi * yj * ch / 2.0;
                aj2[i] += yj * ch;
            }
            else
            {
                double ax = 0.0;
                if((i == 14 && j == 15) || (i == 15 && j == 14))
                {
                    ax = harmonicMean(params.a[i], params.a[j]);
                }
                else if(options.iavg == 1)
                {
                    ax = std::sqrt(params.a[i] * params.a[j]);
                }
                else if(options.iavg == 2)
                {
                    ax = (params.a[i] + params.a[j]) / 2.0;
                }
                else
                {
                    ax = harmonicMean(params.a[i], params.a[j]);
                }

                aij += yi * yj * ax;
                aj2[i] += 2.0 * yj * ax;
            }
        }
    }

    const double c1 = -rt / pressureBar;
    const double c3 = -aij * bx / pressureBar / dsqrtt;
    const double c2 = c1 * bx + aij / dsqrtt / pressureBar - bx * bx;

    result.roots = roots3(c1, c2, c3);

    const auto& roots = result.roots;
    double vol = 0.0;

    if(rootState && rootState->sroot)
    {
        if(rootState->irt == 3 && roots.iroots == 3 && roots.ineg == 0 && roots.vmin > bx)
        {
            vol = rootState->max ? roots.vmax : roots.vmin;
        }
        else if(roots.iroots == 3 && rootState->irt == 3)
        {
            vol = roots.vmax;
        }
        else
        {
            double dv = 1.0e99;
            int jrt = 0;
            for(int i = 0; i < roots.iroots; ++i)
            {
                const double ev = roots.roots[i];
                if(ev < 0.0) continue;
                const double diff = std::abs(ev - rootState->vrt);
                if(diff < dv)
                {
                    jrt = i;
                    dv = diff;
                }
            }

            vol = (dv == 1.0e99) ? roots.roots[0] : roots.roots[jrt];
        }
    }
    else if(roots.iroots == 3 && roots.ineg == 0 && roots.vmin > bx)
    {
        const double pdv = pressureBar * (roots.vmax - roots.vmin)
            - std::log((roots.vmax - bx) / (roots.vmin - bx)) * rt
            - std::log((roots.vmax + bx) / (bx + roots.vmin) * roots.vmin / roots.vmax) * aij / bx / dsqrtt;

        if(pdv > 0.0)
        {
            vol = roots.vmin;
            if(rootState) rootState->max = false;
        }
        else
        {
            vol = roots.vmax;
            if(rootState) rootState->max = true;
        }
    }
    else if(roots.iroots == 3)
    {
        vol = roots.vmax;
    }
    else
    {
        vol = roots.roots[roots.ipos - 1];
    }

    if(rootState && !rootState->sroot)
    {
        rootState->irt = roots.iroots;
        rootState->vrt = vol;
    }

    result.vol = vol;

    const double vpb = vol + bx;
    const double vmb = vol - bx;
    const double d1 = rt * dsqrtt * bx;
    const double d2 = std::log(vpb / vol) / d1;
    const double d3 = aij * d2 / bx - aij / vpb / d1 + 1.0 / vmb;
    const double d6 = std::log(rt / vmb);
    const double d4 = vmb * vmb / vpb / (rt * dsqrtt) / vol;
    const double d5 = aij * d4 * (1.0 / vol + 1.0 / vpb) - 1.0;
    const double d7 = -d4 * aij / vpb;

    for(int k = 0; k < static_cast<int>(species.size()); ++k)
    {
        const int l = species[k];
        const double yl = (y[l] > options.minY) ? y[l] : options.minY;
        const double f = std::log(yl) + params.b[l] * d3 - aj2[l] * d2 + d6;
        result.ln_f[l] = f;
        result.g[l] = std::exp(f) / pressureBar / yl;
        result.v[l] = (d4 * aj2[l] - params.b[l] - vmb + d7 * params.b[l]) / d5;
    }

    return result;
}

void applyHybridFugacity(std::array<double, 19>& g,
                         const std::vector<int>& hybridSpecies,
                         const std::array<double, 19>& gh)
{
    for(const int j : hybridSpecies)
        g[j] = gh[j] * g[j];
}

} // namespace Reaktoro::PerpleX
