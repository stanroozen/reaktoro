// WaterSolventFunctionDEW.cpp

#include <Reaktoro/Extensions/DEW/WaterSolventFunctionDEW.hpp>
#include <Reaktoro/Extensions/DEW/WaterPsatPolynomialsDEW.hpp>
#include <Reaktoro/Extensions/DEW/WaterHelmholtzPropsWagnerPruss.hpp> // for Psat(T) if needed
#include <Reaktoro/Extensions/DEW/WaterUtils.hpp>

#include <cmath>

namespace Reaktoro {
namespace {

inline double toCelsius(double T_K) { return T_K - 273.15; }

// density: kg/m3 -> g/cm3
inline double rho_si_to_g_cm3(double rho_si) { return rho_si / 1000.0; }

// Convert DP from (kg/m3)/Pa to (g/cm3)/bar for DEW-style formula.
inline double drho_si_dP_to_g_cm3_per_bar(double DP_si)
{
    // ρ_g = ρ_SI / 1000
    // dρ_g/dP_bar = dρ_SI/dP_Pa * (1/1000) * (1e5) = DP_si * 100
    return DP_si * 100.0;
}

} // namespace

//----------------------------------------------------------------------------//
// g(P, T, ρ) - DEW / Shock et al. solvent function
//----------------------------------------------------------------------------//

auto waterSolventFunctionDEW(real T,
                             real P,
                             const WaterThermoProps& wt,
                             const WaterSolventFunctionOptions& opt)
    -> real
{
    const double T_K = static_cast<double>(T);
    const double T_C = toCelsius(T_K);

    // --- Psat branch --------------------------------------------------------
    //
    // Excel itself does not introduce a separate polynomial for g(Psat),
    // it just uses the same a_g, b_g, f(P,T) formula evaluated with
    // the Psat density. For completeness, we do exactly that here.
    //
    if (opt.Psat)
    {
        // Saturation density from DEW Psat polynomial (kg/m3 -> g/cm3)
        const double rho_psat_kg_m3 = waterPsatDensityDEW(T_K);
        const double rho_psat_g_cm3 = rho_psat_kg_m3 / 1000.0;

        if (rho_psat_g_cm3 >= 1.0)
            return 0.0;

        // Psat pressure from Wagner-Pruss (consistent with rest of framework)
        const double Psat_Pa  = waterSaturationPressureWagnerPruss(T_K);
        const double Psat_bar = Psat_Pa * 1.0e-5;

        const double a_g =
            -2.037662
            + 0.005747      * T_C
            - 6.557892e-6   * T_C * T_C;

        const double b_g =
            6.107361
            - 0.01074377    * T_C
            + 1.268348e-5   * T_C * T_C;

        double f = 0.0;
        if (Psat_bar <= 1000.0 && T_C >= 155.0 && T_C <= 355.0)
        {
            const double x     = (T_C - 155.0) / 300.0;
            const double x_4_8 = std::pow(x, 4.8);
            const double x_16  = std::pow(x, 16.0);

            f = (x_4_8 + 36.66666 * x_16) *
                (-1.504956e-10 * std::pow(1000.0 - Psat_bar, 3.0)
                 + 5.017997e-14 * std::pow(1000.0 - Psat_bar, 4.0));
        }

        const double one_minus_rho = 1.0 - rho_psat_g_cm3;
        const double g =
            a_g * std::pow(one_minus_rho, b_g) - f;

        return g;
    }

    // --- General branch (Psat = False) -------------------------------------

    const double P_bar = static_cast<double>(P) * 1.0e-5;
    const double rho_g = rho_si_to_g_cm3(wt.D);

    if (rho_g >= 1.0)
        return 0.0; // Excel: if density >= 1 g/cm3, g = 0

    const double a_g =
        -2.037662
        + 0.005747      * T_C
        - 6.557892e-6   * T_C * T_C;

    const double b_g =
        6.107361
        - 0.01074377    * T_C
        + 1.268348e-5   * T_C * T_C;

    double f = 0.0;
    if (P_bar <= 1000.0 && T_C >= 155.0 && T_C <= 355.0)
    {
        const double x     = (T_C - 155.0) / 300.0;
        const double x_4_8 = std::pow(x, 4.8);
        const double x_16  = std::pow(x, 16.0);

        f = (x_4_8 + 36.66666 * x_16) *
            (-1.504956e-10 * std::pow(1000.0 - P_bar, 3.0)
             + 5.017997e-14 * std::pow(1000.0 - P_bar, 4.0));
    }

    const double one_minus_rho = 1.0 - rho_g;
    const double g =
        a_g * std::pow(one_minus_rho, b_g) - f;

    return g;
}

//----------------------------------------------------------------------------//
// d g / dP - DEW / Shock et al.
//----------------------------------------------------------------------------//

auto waterSolventFunctionDgdP_DEW(real T,
                                  real P,
                                  const WaterThermoProps& wt,
                                  real g,
                                  const WaterSolventFunctionOptions& opt)
    -> real
{
    const double T_K = static_cast<double>(T);
    const double T_C = toCelsius(T_K);

    // Psat branch: use dedicated DEW polynomial (already SI-converted).
    if (opt.Psat)
    {
        return waterPsatDgdPDEW(T_K); // [1/Pa]
    }

    const double P_bar = static_cast<double>(P) * 1.0e-5;
    const double rho_g = rho_si_to_g_cm3(wt.D);

    if (rho_g >= 1.0)
        return 0.0; // Excel: if density >= 1 g/cm3, dgdP = 0

    const double b_g =
        6.107361
        - 0.01074377    * T_C
        + 1.268348e-5   * T_C * T_C;

    // d f / dP (only non-zero in same T/P window)
    double dfdP_bar = 0.0;
    if (P_bar <= 1000.0 && T_C >= 155.0 && T_C <= 355.0)
    {
        const double x     = (T_C - 155.0) / 300.0;
        const double x_4_8 = std::pow(x, 4.8);
        const double x_16  = std::pow(x, 16.0);

        dfdP_bar =
            -(x_4_8 + 36.66666 * x_16) *
            (3.0 * -1.504956e-10 * std::pow(1000.0 - P_bar, 2.0)
             + 4.0 * 5.017997e-14 * std::pow(1000.0 - P_bar, 3.0));
    }

    const double drho_g_dP_bar = drho_si_dP_to_g_cm3_per_bar(wt.DP);
    const double one_minus_rho = 1.0 - rho_g;

    // Excel:
    //   dgdP_bar = -b_g * drhodP_bar * g / (1 - rho) - dfdP_bar
    const double dgdP_bar =
        -b_g * drho_g_dP_bar * g / one_minus_rho - dfdP_bar;

    // Convert 1/bar -> 1/Pa
    return dgdP_bar * 1.0e-5;
}

} // namespace Reaktoro
