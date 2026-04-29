#include "PerpleXHKF.hpp"
#include "PerpleXMrkPure.hpp"

#ifndef PERPLEX_STANDALONE
#  include <Reaktoro/Extensions/DEW/WaterEosZhangDuan2005.hpp>
#endif

#include <cmath>
#include <stdexcept>
#include <vector>

namespace Reaktoro::PerpleX {

double calculateBornOmega(double charge,
                          double bornRadius,
                          double omega0,
                          double gf)
{
    using namespace HKFConstants;

    if (std::abs(charge) < 1e-10)
    {
        // Neutral species - omega is constant
        return omega0;
    }

    // Ionic species - Born function
    // omega = eta * z * (z/(r_e + |z|*gf) - 1/(3.082 + gf))
    // Based on rlib.f lines 3014-3015

    double z = charge;
    double absZ = std::abs(z);

    double omega = ETA * z * (z / (bornRadius + absZ * gf)
                            - 1.0 / (NEUTRAL_RADIUS + gf));

    return omega;
}

HKFState computeHKFGibbs(const HKFParams& params,
                         double pressureBar,
                         double temperatureK,
                         double epsilon,
                         double gf)
{
    using namespace HKFConstants;

    HKFState state;

    // Calculate Born omega at P-T
    state.omega = calculateBornOmega(params.charge,
                                     params.bornRadius,
                                     params.omega0,
                                     gf);

    // Temperature and pressure terms
    double ft = temperatureK - THETA;
    double fp = std::log(PSI + pressureBar);

    // HKF equation (rlib.f lines 3026-3031)
    // G = b9 + (b8 + b12*ln(ft) + b13*ln(T))*T + b11*ft
    //   + a1*P + a2*ln(psi+P) + (a3*P + a4*ln(psi+P) + b10)/ft
    //   + omega*(1/ε - 1) - omega0/ε0

    state.G = params.b9
            + (params.b8 + params.b12 * std::log(ft) + params.b13 * std::log(temperatureK)) * temperatureK
            + params.b11 * ft
            + params.a1 * pressureBar
            + params.a2 * fp
            + (params.a3 * pressureBar + params.a4 * fp + params.b10) / ft
            + state.omega * (1.0 / epsilon - 1.0)
            - params.omega0 / EPSILON0;

    state.bornTerm = state.omega * (1.0 / epsilon - 1.0) - params.omega0 / EPSILON0;

    return state;
}

HKFParams preprocessHKFParams(const HKFParams& params,
                              double Tr,
                              double Pr)
{
    using namespace HKFConstants;

    HKFParams processed = params;

    // Pre-compute coefficients as done in Perple_X
    // These formulas are from the comments in rlib.f lines 2960-2978
    //
    // IMPORTANT: yr is the Born Y reference constant (dZ/dT at Tr, Pr where Z=-1/epsilon).
    // Aligned to DEW/SUPCRT92: StandardThermoModelDEW.cpp Yr = -5.795424563e-05.
    // (Perple_X tlib.f uses -5.79865e-5, a truncated approximation of the same J&N1991 value.)
    // It is NOT 1/(theta-Tr) = -0.01426 — that would cause ~Megajoule errors at high T.
    constexpr double yr = -5.795424563e-5;  // Born Y reference constant [K^-1], SUPCRT92/DEW/J&N1991

    // b8 = -S0 + c1*ln(Tr) + c1 + omega0*yr + ln(Tr/(Tr-theta))*c2/theta^2
    processed.b8 = -params.S0
                 + params.c1 * std::log(Tr)
                 + params.c1
                 + params.omega0 * yr
                 + std::log(Tr / (Tr - THETA)) * params.c2 / (THETA * THETA);

    // b9 = (-omega0*yr - c1 + S0)*Tr + omega0 - a1*Pr - a2*ln(psi+Pr) + G0 + c2/theta
    double fpRef = std::log(PSI + Pr);
    processed.b9 = (-params.omega0 * yr - params.c1 + params.S0) * Tr
                 + params.omega0
                 - params.a1 * Pr
                 - params.a2 * fpRef
                 + params.G0
                 + params.c2 / THETA;

    // b10 = -a3*Pr - a4*ln(psi+Pr)
    processed.b10 = -params.a3 * Pr - params.a4 * fpRef;

    // b11 = -c2/(Tr-theta)/theta
    processed.b11 = -params.c2 / ((Tr - THETA) * THETA);

    // b12 = c2/theta^2
    processed.b12 = params.c2 / (THETA * THETA);

    // b13 = -c1 - c2/theta^2
    processed.b13 = -params.c1 - params.c2 / (THETA * THETA);

    // Born radius calculation (from thermo(19) assignment in tlib.f line ~11286):
    //   b9 = 5d9 * eta * z^2 / (1.622323167d9 * eta * z + 5d9 * omega0)
    // Equivalent to: re = eta*z^2 / (eta*z/3.082 + omega0)
    // NOTE: rlib.f comment mistakenly shows 5d10/1622323167 — tlib.f is the source of truth.
    if (std::abs(params.charge) > 1e-10)
    {
        double q = params.charge;
        double q2 = q * q;
        double numerator = 5e9 * ETA * q2;
        double denominator = 1.622323167e9 * ETA * q + 5e9 * params.omega0;
        processed.bornRadius = numerator / denominator;
        processed.chargeSquared = q2;
    }
    else
    {
        processed.bornRadius = NEUTRAL_RADIUS;
        processed.chargeSquared = 0.0;
    }

    return processed;
}

double waterDensity(double pressureBar, double temperatureK)
{
    // Simple water density model for g-function calculation
    // For accurate results, should use IAPWS-95 or similar
    // This is a simplified approximation sufficient for g-function

    // Critical point parameters
    constexpr double Tc = 647.096;  // K
    constexpr double rhoc = 0.322;  // g/cm³

    // Simple corresponding states approximation
    double Tr = temperatureK / Tc;

    if (Tr < 1.0)
    {
        // Liquid branch - approximate correlation
        // rho ≈ rhoc * (1 + 1.9*(1-Tr)^(1/3))
        double tau = 1.0 - Tr;
        double rho = rhoc * (1.0 + 1.9 * std::pow(tau, 1.0/3.0));

        // Pressure correction (simplified)
        // drho/dP ≈ kappa * rho where kappa ≈ 4.5e-5 /bar
        double kappa = 4.5e-5;
        rho *= (1.0 + kappa * pressureBar);

        return rho;
    }
    else
    {
        // Supercritical - use ideal gas approximation with real gas correction
        constexpr double R = 83.14472;  // cm³·bar/(mol·K)
        constexpr double M = 18.01528;  // g/mol

        // Compressibility factor (approximate)
        double Z = 1.0 + pressureBar / (100.0 * temperatureK);

        double rho = pressureBar * M / (Z * R * temperatureK);
        return rho;
    }
}

DielectricState getWaterSolventState(double pressureBar,
                                     double temperatureK,
                                     double& waterVolume)
{
    using namespace HKFConstants;

    DielectricState state;

    constexpr double M_H2O = 18.01528;  // g/mol

    // Prefer Zhang-Duan 2005 EOS in the full Reaktoro build.
    // Standalone Perple_X regression builds use the internal fallback density model.
    double rho = 0.0;
#ifndef PERPLEX_STANDALONE
    {
        const auto wtp = waterThermoPropsZhangDuan2005(temperatureK, pressureBar * 1e5);
        const double rho_kgm3 = static_cast<double>(wtp.D);

        if (std::isfinite(rho_kgm3) && rho_kgm3 > 0.0)
        {
            rho = rho_kgm3 * 1e-3;           // kg/m³ -> g/cm³
            waterVolume = M_H2O / rho;       // cm³/mol
        }
        else
        {
            rho = waterDensity(pressureBar, temperatureK);
            waterVolume = M_H2O / rho;
        }
    }
#else
    rho = waterDensity(pressureBar, temperatureK);
    waterVolume = M_H2O / rho;
#endif

    state.vsolv = waterVolume;
    state.msol = M_H2O * 1e-3;  // kg/mol

    // Dielectric constant for pure water
    double v_jbar = waterVolume / 10.0;  // Convert to J/bar
    state.epsilon = epsh2o(v_jbar, temperatureK);

    // Debye-Hückel factor
    // adh = cdh / sqrt(v_h2o/10 * (ε*T)^3)
    // CDH = -5661800.47810 (from slvnt0)
    constexpr double CDH_SLVNT0 = -5661800.47810;
    state.adh = CDH_SLVNT0 / std::sqrt(v_jbar * std::pow(state.epsilon * temperatureK, 3));

    // Shock g-function (CGS solvent density)
    state.gf = gfunc(rho, pressureBar, temperatureK);

    // Set reference values
    state.epsilon = state.epsilon;  // Already computed above

    return state;
}

} // namespace Reaktoro::PerpleX
