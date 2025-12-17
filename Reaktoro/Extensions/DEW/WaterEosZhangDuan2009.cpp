#include "WaterEosZhangDuan2009.hpp"

// C++ includes
#include <cmath>
#include <algorithm>

// Reaktoro includes
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>

namespace Reaktoro {
namespace {

using std::pow;
using std::exp;
using std::abs;

// -----------------------------------------------------------------------------
// Constants from DEW Excel/VBA (equation = 2, Zhang & Duan 2009 branch)
// -----------------------------------------------------------------------------

// Gas constant [dm3·bar / (mol·K)]
constexpr double ZD09_R  = 0.083145;

// Scaling constant c1 = epsilon / (3.0626 * omega^3)
// (copied from DEW code)
constexpr double ZD09_c1 = 6.971118009;

// NOTE: The VBA snippet references ZD09_c4 in calculate_drhodP but never defines it.
// To stay structurally faithful while remaining numerically defined, we set:
constexpr double ZD09_c4 = 1.0;

// Molar mass of water [g/mol]
constexpr double H2O_M = 18.01528;

// -----------------------------------------------------------------------------
// Helpers & unit conversions
// -----------------------------------------------------------------------------

inline double kelvinToCelsius(double T) { return T - 273.15; }
inline double pascalToBar(double P)     { return P * 1.0e-5; }

// Density along Psat curve (same polynomial as in DEW Psat=True branch)
inline double densityPsatPoly_DEW(double T_C)
{
    // From your DEW code: Psat density polynomial
    // (originally used in calculateDensity when Psat=True).
    return
        -1.01023381581205e-104 * pow(T_C, 40.0) +
        -1.13685997859530e-27  * pow(T_C, 10.0) +
        -2.11689207168779e-11  * pow(T_C,  4.0) +
         1.26878850169523e-08  * pow(T_C,  3.0) +
        -4.92010672693621e-06  * pow(T_C,  2.0) +
        -3.26665986126920e-05  * T_C +
         1.00046144613017; // g/cm3
}

// -----------------------------------------------------------------------------
// Pressure as function of density and T for ZD09 (calculatePressure, Case 2)
// -----------------------------------------------------------------------------

inline double calculatePressure_ZD09(double rho_g_cm3, double T_C)
{
    // Direct translation of Case 2 (Zhang & Duan 2009) in calculatePressure.
    //
    // Input:
    //   rho_g_cm3 : density [g/cm3]
    //   T_C       : temperature [°C]
    // Output:
    //   P in [bar]

    const double m  = H2O_M;             // [g/mol]
    const double TK = T_C + 273.15;      // [K]

    // Prefactor-converted quantities (as in DEW code comments)
    const double dm = 475.05656886 * rho_g_cm3;
    const double Vm = 0.0021050125 * (m / rho_g_cm3);
    const double Tm = 0.3019607843 * TK;

    const double B =
          0.029517729893
        - 6337.56452413  / (Tm * Tm)
        - 275265.428882  / (Tm * Tm * Tm);

    const double C =
          0.00129128089283
        - 145.797416153   / (Tm * Tm)
        + 76593.8947237   / (Tm * Tm * Tm);

    const double D =
          2.58661493537e-06
        + 0.52126532146   / (Tm * Tm)
        - 139.839523753   / (Tm * Tm * Tm);

    const double E =
        -2.36335007175e-08
        + 0.00535026383543 / (Tm * Tm)
        - 0.27110649951    / (Tm * Tm * Tm);

    const double f =
        25038.7836486 / (Tm * Tm * Tm);

    const double Vm2 = Vm * Vm;
    const double Vm4 = Vm2 * Vm2;
    const double Vm5 = Vm4 * Vm;

    const double expterm = exp(-0.015483335997 / Vm2);

    const double delta =
          1.0
        + B / Vm
        + C / Vm2
        + D / Vm4
        + E / Vm5
        + f / Vm2
          * (0.73226726041 + 0.015483335997 / Vm2)
          * expterm;

    // Pm in bar (since ZD09_R is in dm3·bar/(mol·K))
    const double Pm_bar = ZD09_R * Tm * delta / Vm;

    // Final pressure in bar as in VBA:
    return Pm_bar * ZD09_c1;
}

// -----------------------------------------------------------------------------
// Density as function of P and T: bisection (calculateDensity, equation=2)
// -----------------------------------------------------------------------------

inline double calculateDensity_ZD09(double P_bar,
                                    double T_C,
                                    const WaterZhangDuan2009Options& opts)
{
    if (opts.usePsat)
    {
        // Use Psat polynomial fit (same as DEW Psat density).
        return densityPsatPoly_DEW(T_C); // g/cm3
    }

    // Bisection method as in DEW:
    //   minGuess = 1e-5
    //   maxGuess = 7.5 * equation - 5 = 10.0 for equation=2
    //   50 iterations max

    double minGuess = 1.0e-5;   // g/cm3
    double maxGuess = 7.5 * 2.0 - 5.0; // = 10.0 g/cm3
    double rho      = minGuess;

    for (int i = 0; i < opts.maxIterations; ++i)
    {
        const double P_calc = calculatePressure_ZD09(rho, T_C);
        const double diff   = P_calc - P_bar;

        if (std::abs(diff) <= opts.pressureToleranceBar)
            break;

        if (P_calc > P_bar)
        {
            maxGuess = rho;
            rho = 0.5 * (rho + minGuess);
        }
        else
        {
            minGuess = rho;
            rho = 0.5 * (rho + maxGuess);
        }
    }

    return rho; // g/cm3
}

// -----------------------------------------------------------------------------
// (d rho / dP)_T for ZD09 (calculate_drhodP, Case 2)
// -----------------------------------------------------------------------------

inline double calculate_drhodP_ZD09(double rho_g_cm3, double T_C)
{
    // Direct translation of the Case 2 (Zhang & Duan 2009) branch
    // in calculate_drhodP from your DEW Excel code.
    //
    // NOTE: ZD09_c4 is undefined in the snippet; here it is set to 1.0
    //       to preserve the original structure m/(R*T*(...)) with the
    //       additional ZD09_c1 factor.

    const double m  = H2O_M;           // g/mol
    const double TK = T_C + 273.15;    // K

    const double dm = 475.05656886 * rho_g_cm3;
    const double Vm = 0.0021050125 * (m / rho_g_cm3);
    const double Tm = 0.3019607843 * TK;

    const double B =
          0.029517729893
        - 6337.56452413  / (Tm * Tm)
        - 275265.428882  / (Tm * Tm * Tm);

    const double C =
          0.00129128089283
        - 145.797416153   / (Tm * Tm)
        + 76593.8947237   / (Tm * Tm * Tm);

    const double D =
          2.58661493537e-06
        + 0.52126532146   / (Tm * Tm)
        - 139.839523753   / (Tm * Tm * Tm);

    const double E =
        -2.36335007175e-08
        + 0.00535026383543 / (Tm * Tm)
        - 0.27110649951    / (Tm * Tm * Tm);

    const double f =
        25038.7836486 / (Tm * Tm * Tm);

    const double Vm2 = Vm * Vm;
    const double Vm4 = Vm2 * Vm2;
    const double Vm5 = Vm4 * Vm;

    const double expterm = exp(-0.015483335997 / Vm2);

    const double delta =
          1.0
        + B / Vm
        + C / Vm2
        + D / Vm4
        + E / Vm5
        + f / Vm2
          * (0.73226726041 + 0.015483335997 / Vm2)
          * expterm;

    const double kappa =
        B / m
        + 2.0 * C * dm / (m * m)
        + 4.0 * D * pow(dm, 3.0) / pow(m, 4.0)
        + 5.0 * E * pow(dm, 4.0) / pow(m, 5.0)
        + (
              2.0 * f * dm / (m * m)
                * (0.73226726041 + 0.015483335997 / Vm2)
            + f / Vm2
                * (1.0 - 0.73226726041 - 0.015483335997 / Vm2)
                * (2.0 * 0.015483335997 * dm / (m * m))
          ) * expterm;

    // VBA:
    //   calculate_drhodP = ZD09_c1 * m /
    //                      (ZD09_c4 * ZD09_R * Tm * (delta + dm * kappa))
    //
    // Units: (g/cm3)/bar (by consistency with ZD05 branch)
    const double drho_dP_bar =
        ZD09_c1 * m / (ZD09_c4 * ZD09_R * Tm * (delta + dm * kappa));

    return drho_dP_bar; // g/cm3 per bar
}

} // namespace (end anonymous namespace)

// -----------------------------------------------------------------------------
// Public API wrapper: waterThermoPropsZhangDuan2009
// -----------------------------------------------------------------------------

auto waterThermoPropsZhangDuan2009(real T, real P,
                                   const WaterZhangDuan2009Options& opts)
    -> WaterThermoProps
{
    WaterThermoProps wt;

    // Convert external to DEW units
    const double T_C   = kelvinToCelsius(T);
    const double P_bar = pascalToBar(P);

    // Density from DEW-style bisection
    const double rho_g_cm3 = calculateDensity_ZD09(P_bar, T_C, opts);

    // Convert to kg/m3
    const double rho_kg_m3 = rho_g_cm3 * 1000.0;
    wt.D = rho_kg_m3;

    // (d rho / dP)_T from analytic DEW expression
    const double drho_dP_bar_g_cm3 = calculate_drhodP_ZD09(rho_g_cm3, T_C);

    // Convert (g/cm3)/bar -> (kg/m3)/Pa:
    // factor = 1000 / 1e5 = 1e-2
    wt.DP = drho_dP_bar_g_cm3 * 1.0e-2;

    // Leave other derivatives zero here: no extra theory beyond the
    // provided Excel code is assumed.
    wt.DT  = 0.0;
    wt.DTT = 0.0;
    wt.DTP = 0.0;
    wt.DPP = 0.0;

    return wt;
}

} // namespace Reaktoro
