#include "WaterEosZhangDuan2005.hpp"

#include <cmath>

namespace Reaktoro {
namespace {

// -----------------------------------------------------------------------------
// Exact constants and formulas from DEW Excel (equation = 1, Zhang & Duan 2005)
// -----------------------------------------------------------------------------

// Molar mass of water [g/mol]
constexpr double M_H2O = 18.01528;

// Zhang & Duan (2005) constants (from Excel code)
constexpr double ZD05_R  = 83.144;        // [cm3·bar/(mol·K)]
constexpr double ZD05_Vc = 55.9480373;    // [cm3/mol]
constexpr double ZD05_Tc = 647.25;        // [K]

// Conversion helpers
inline double bar_from_P_Pa(double P_Pa)
{
    return P_Pa * 1.0e-5; // 1 bar = 1e5 Pa
}

inline double rho_kg_m3_from_g_cm3(double rho_g_cm3)
{
    return rho_g_cm3 * 1000.0; // 1 g/cm3 = 1000 kg/m3
}

inline double drho_kg_m3_per_Pa_from_g_cm3_per_bar(double drho_g_cm3_per_bar)
{
    // (g/cm3)/bar * (1000 kg/m3)/(1 g/cm3) * (1 bar / 1e5 Pa)
    // = drho * (1000 / 1e5) = drho * 0.01
    return drho_g_cm3_per_bar * 0.01;
}

// -----------------------------------------------------------------------------
// Pressure as function of density and temperature: calculatePressure (equation=1)
// Direct translation of the Case 1 block in Excel `calculatePressure`.
// Inputs:
//  - rho_g_cm3 : density [g/cm3]
//  - T_C       : temperature [°C]
// Output:
//  - P_bar     : pressure [bar]
// -----------------------------------------------------------------------------
double zd05_pressure_bar(double rho_g_cm3, double T_C)
{
    const double T_K = T_C + 273.15;

    const double Vr = M_H2O / (rho_g_cm3 * ZD05_Vc);
    const double Tr = T_K / ZD05_Tc;

    const double B = 0.349824207
        - 2.91046273 / (Tr*Tr)
        + 2.00914688 / (Tr*Tr*Tr);

    const double C = 0.112819964
        + 0.748997714 / (Tr*Tr)
        - 0.87320704 / (Tr*Tr*Tr);

    const double D = 0.0170609505
        - 0.0146355822 / (Tr*Tr)
        + 0.0579768283 / (Tr*Tr*Tr);

    const double E = -0.000841246372
        + 0.00495186474 / (Tr*Tr)
        - 0.00916248538 / (Tr*Tr*Tr);

    const double f = -0.100358152 / Tr;
    const double g = -0.00182674744 * Tr;

    const double delta =
          1.0
        + B / Vr
        + C / (Vr*Vr)
        + D / std::pow(Vr, 4)
        + E / std::pow(Vr, 5)
        + (f / (Vr*Vr) + g / std::pow(Vr, 4))
            * std::exp(-0.0105999998 / (Vr*Vr));

    const double P_bar = ZD05_R * T_K * rho_g_cm3 * delta / M_H2O;

    return P_bar;
}

// -----------------------------------------------------------------------------
// Density as function of pressure and temperature: calculateDensity (equation=1)
// Non-Psat branch only (Psat=False), as used by DEW for general P-T.
// This is the bisection method from Excel:
//  - minGuess = 1e-5 g/cm3
//  - maxGuess = 7.5*equation - 5 = 2.5 g/cm3 for equation=1
//  - up to 50 iterations
//  - tolerance `error` in bar, default 0.01
// -----------------------------------------------------------------------------
double zd05_density_g_cm3(double P_bar_target, double T_C,
                          double error_bar = 0.01)
{
    double rho_min = 1.0e-5;   // [g/cm3]
    double rho_max = 2.5;      // [g/cm3] as in Excel: 7.5*1 - 5
    double rho     = rho_min;

    for(int iter = 0; iter < 50; ++iter)
    {
        const double P_bar = zd05_pressure_bar(rho, T_C);
        const double diff  = P_bar - P_bar_target;

        if(std::fabs(diff) <= error_bar)
            return rho;

        if(diff > 0.0)
        {
            rho_max = rho;
            rho = 0.5 * (rho + rho_min);
        }
        else
        {
            rho_min = rho;
            rho = 0.5 * (rho + rho_max);
        }
    }

    // Return last iterate if not converged inside tolerance (Excel does same effectively)
    return rho;
}

// -----------------------------------------------------------------------------
// (∂ρ/∂P)_T : calculate_drhodP (equation=1)
// Direct translation of the Case 1 block in Excel `calculate_drhodP` for
// Zhang & Duan (2005).
//
// Inputs:
//  - rho_g_cm3 : density [g/cm3]
//  - T_C       : temperature [°C]
// Output:
//  - drho/dP in [g/cm3/bar]
// -----------------------------------------------------------------------------
double zd05_drhodP_g_cm3_per_bar(double rho_g_cm3, double T_C)
{
    const double T_K = T_C + 273.15;
    const double Tr  = T_K / ZD05_Tc;

    const double cc = ZD05_Vc / M_H2O;       // appears repeatedly
    const double Vr = M_H2O / (rho_g_cm3 * ZD05_Vc);

    const double B = 0.349824207
        - 2.91046273 / (Tr*Tr)
        + 2.00914688 / (Tr*Tr*Tr);

    const double C = 0.112819964
        + 0.748997714 / (Tr*Tr)
        - 0.87320704 / (Tr*Tr*Tr);

    const double D = 0.0170609505
        - 0.0146355822 / (Tr*Tr)
        + 0.0579768283 / (Tr*Tr*Tr);

    const double E = -0.000841246372
        + 0.00495186474 / (Tr*Tr)
        - 0.00916248538 / (Tr*Tr*Tr);

    const double f = -0.100358152 / Tr;
    const double g =  0.0105999998 * Tr; // NOTE: this g is different from calcPressure g

    const double Vr2 = Vr*Vr;
    const double Vr4 = Vr2*Vr2;
    const double Vr5 = Vr4*Vr;

    const double expterm = std::exp(-0.0105999998 / Vr2);

    const double delta =
          1.0
        + B / Vr
        + C / Vr2
        + D / Vr4
        + E / Vr5
        + (f / Vr2 + g / Vr4) * expterm;

    // kappa expression from Excel (copied exactly, just spaced)
    const double rho = rho_g_cm3;

    const double kappa =
          B * cc
        + 2.0 * C * (cc*cc) * rho
        + 4.0 * D * std::pow(cc, 4) * std::pow(rho, 3)
        + 5.0 * E * std::pow(cc, 5) * std::pow(rho, 4)
        + (
              2.0 * f * (cc*cc) * rho
            + 4.0 * g * std::pow(cc, 4) * std::pow(rho, 3)
            - (f / Vr2 + g / Vr4) * (2.0 * 0.0105999998 * (cc*cc) * rho)
          ) * expterm;

    // Excel:
    // calculate_drhodP = m / (ZD05_R * TK * (delta + density * kappa))
    const double drhodP_g_cm3_per_bar =
        M_H2O / (ZD05_R * T_K * (delta + rho * kappa));

    return drhodP_g_cm3_per_bar;
}

} // namespace

// -----------------------------------------------------------------------------
// Public interface: waterThermoPropsZhangDuan2005
// -----------------------------------------------------------------------------
auto waterThermoPropsZhangDuan2005(real T, real P, double densityTolerance) -> WaterThermoProps
{
    WaterThermoProps wt;

    const double T_K   = static_cast<double>(T);
    const double T_C   = T_K - 273.15;
    const double P_bar = bar_from_P_Pa(static_cast<double>(P));

    // Density from exact DEW bisection logic
    const double rho_g_cm3 = zd05_density_g_cm3(P_bar, T_C, densityTolerance);
    const double rho_kg_m3 = rho_kg_m3_from_g_cm3(rho_g_cm3);

    // (∂ρ/∂P)_T from exact DEW analytic expression
    const double drho_g_cm3_per_bar = zd05_drhodP_g_cm3_per_bar(rho_g_cm3, T_C);
    const double drho_kg_m3_per_Pa  = drho_kg_m3_per_Pa_from_g_cm3_per_bar(drho_g_cm3_per_bar);

    wt.D  = rho_kg_m3;
    wt.DP = drho_kg_m3_per_Pa;

    // The Excel snippet you provided does not define (∂ρ/∂T)_P or higher T-derivatives
    // for Zhang & Duan (2005). To keep this module an exact translation of that code,
    // we do NOT invent additional formulas here.
    // If Reaktoro needs DT, DTT, DTP, DPP, etc., those can be added separately using
    // consistent finite-difference approximations around this EOS.
    wt.DT  = 0.0;
    wt.DTT = 0.0;
    wt.DTP = 0.0;
    wt.DPP = 0.0;

    return wt;
}

} // namespace Reaktoro
