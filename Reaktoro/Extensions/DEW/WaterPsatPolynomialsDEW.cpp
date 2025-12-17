#include <Reaktoro/Extensions/DEW/WaterPsatPolynomialsDEW.hpp>

#include <cmath>

namespace Reaktoro {
namespace {

inline double toCelsius(double T_K)
{
    return T_K - 273.15;
}

/// Clamp helper for very low or unphysical temperatures if needed.
inline bool isValidForPsatPolys(double T_C)
{
    // Excel polynomials were fit over a finite range (~ -20 to 1000 °C region
    // depending on the function). We do not hard clamp here, but callers
    // should only use within their calibration ranges.
    return true;
}

} // namespace

//------------------------------------------------------------------------------
// rho_l(T) along Psat
// Excel (Psat=True in calculateDensity):
//   rho_g_cm3 =
//     -1.01023381581205E-104 * T^40
//     -1.1368599785953E-27  * T^10
//     -2.11689207168779E-11 * T^4
//      1.26878850169523E-08 * T^3
//     -4.92010672693621E-06 * T^2
//     -3.2666598612692E-05 * T
//      1.00046144613017
// with T in °C, rho in g/cm3.
//------------------------------------------------------------------------------
auto waterPsatDensityDEW(real T_K) -> real
{
    const double T_C = toCelsius(T_K);
    if (!isValidForPsatPolys(T_C))
        return 0.0;

    const double T2  = T_C * T_C;
    const double T3  = T2 * T_C;
    const double T4  = T2 * T2;
    const double T10 = std::pow(T_C, 10.0);
    const double T40 = std::pow(T_C, 40.0);

    const double rho_g_cm3 =
        -1.01023381581205e-104 * T40 +
        -1.1368599785953e-27   * T10 +
        -2.11689207168779e-11  * T4  +
         1.26878850169523e-08  * T3  +
        -4.92010672693621e-06  * T2  +
        -3.2666598612692e-05   * T_C +
         1.00046144613017;

    // Convert g/cm3 -> kg/m3
    return rho_g_cm3 * 1000.0;
}

//------------------------------------------------------------------------------
// epsilon(T) along Psat
// Excel (Psat=True in calculateEpsilon):
//   eps =
//     -1.66686763214295E-77 * T^30
//     -9.02887020379887E-07 * T^3
//      8.4590281449009E-04 * T^2
//     -0.396542037778945   * T
//      87.605024245432
// with T in °C.
//------------------------------------------------------------------------------
auto waterPsatEpsilonDEW(real T_K) -> real
{
    const double T_C = toCelsius(T_K);
    if (!isValidForPsatPolys(T_C))
        return 0.0;

    const double T2  = T_C * T_C;
    const double T3  = T2 * T_C;
    const double T30 = std::pow(T_C, 30.0);

    const double eps =
        -1.66686763214295e-77 * T30 +
        -9.02887020379887e-07 * T3  +
         8.4590281449009e-04  * T2  +
        -0.396542037778945    * T_C +
         87.605024245432;

    return eps;
}

//------------------------------------------------------------------------------
// G(T) along Psat
// Excel (Psat=True in calculateGibbsOfWater):
//   G [cal/mol] =
//     -2.72980941772081E-103 * T^40
//      2.88918186300446E-25  * T^10
//     -2.21891314234246E-08  * T^4
//      3.0912103873633E-05   * T^3
//     -3.20873264480928E-02  * T^2
//     -15.169458452209       * T
//     -56289.0379433809
// with T in °C.
// We convert to J/mol (×4.184).
//------------------------------------------------------------------------------
auto waterPsatGibbsDEW(real T_K) -> real
{
    const double T_C = toCelsius(T_K);
    if (!isValidForPsatPolys(T_C))
        return 0.0;

    const double T2  = T_C * T_C;
    const double T3  = T2 * T_C;
    const double T4  = T2 * T2;
    const double T10 = std::pow(T_C, 10.0);
    const double T40 = std::pow(T_C, 40.0);

    const double G_cal_mol =
        -2.72980941772081e-103 * T40 +
         2.88918186300446e-25  * T10 +
        -2.21891314234246e-08  * T4  +
         3.0912103873633e-05   * T3  +
        -3.20873264480928e-02  * T2  +
        -15.169458452209       * T_C +
        -56289.0379433809;

    // cal/mol -> J/mol
    return G_cal_mol * 4.184;
}

//------------------------------------------------------------------------------
// Q(T) along Psat
// Excel (Psat=True in calculateQ):
//   Q [1/bar] = (poly in T_C) * 1.0E-6
// We convert to 1/Pa: 1/bar = 1e-5 1/Pa => Q_Pa = Q_bar * 1e-5.
//------------------------------------------------------------------------------
auto waterPsatBornQDEW(real T_K) -> real
{
    const double T_C = toCelsius(T_K);
    if (!isValidForPsatPolys(T_C))
        return 0.0;

    const double T2  = T_C * T_C;
    const double T3  = T2 * T_C;
    const double T4  = T2 * T2;
    const double T5  = T4 * T_C;
    const double T6  = T3 * T3;
    const double T20 = std::pow(T_C, 20.0);

    const double poly =
         1.99258688758345e-49 * T20 +
        -4.43690270750774e-14 * T6  +
         4.29110215680165e-11 * T5  +
        -1.07146606081182e-08 * T4  +
         1.09982931856694e-06 * T3  +
         9.60705240954956e-06 * T2  +
         0.642579832259358;

    const double Q_bar = poly * 1.0e-6;  // [1/bar]

    // Convert to 1/Pa
    const double Q_Pa = Q_bar * 1.0e-5;
    return Q_Pa;
}

//------------------------------------------------------------------------------
// dgdP(T) along Psat
// Excel (Psat=True in calculate_dgdP):
//   if T < 0.01 then 0
//   else dgdP [Å/bar] = exp( A*ln(T)^15 + B*ln(T)^10 + C*ln(T) + D ) * 1.0E-6
// where T is in °C in that code.
// Returns dgdP in Å/bar (since g is in Ångströms, P is in bar).
//------------------------------------------------------------------------------
auto waterPsatDgdPDEW(real T_K) -> real
{
    const double T_C = toCelsius(T_K);

    if (T_C < 0.01)
        return 0.0;

    const double lnT = std::log(T_C);

    const double expo =
          1.37105493109451e-10 * std::pow(lnT, 15.0)
        - 1.43605469318795e-06 * std::pow(lnT, 10.0)
        + 26.2649453651117     * lnT
        - 125.108856715714;

    const double dgdP_A_per_bar = std::exp(expo) * 1.0e-6; // [Å/bar]

    return dgdP_A_per_bar; // Å/bar
}

} // namespace Reaktoro
