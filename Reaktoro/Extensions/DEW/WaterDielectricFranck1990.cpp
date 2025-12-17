#include "WaterDielectricFranck1990.hpp"

// C++ includes
#include <cmath>
using std::pow;

// Reaktoro includes
#include <Reaktoro/Extensions/DEW/WaterElectroProps.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>

namespace Reaktoro {
namespace {

// -----------------------------------------------------------------------------
// Constants (identical to DEW Excel/VBA Franck branch)
// -----------------------------------------------------------------------------

constexpr double pi = 3.14159265358979;

// Lennard-Jones distance [cm]
constexpr double omega = 0.0000000268; // 2.68e-8 cm

// Boltzmann constant [erg/K]
constexpr double kB = 1.380648e-16;

// Avogadro's number [1/mol]
constexpr double Na = 6.022e23;

// Dipole moment [statC·cm] (≈ Debye in cgs units)
constexpr double mu = 2.33e-18;

// Conversion: density kg/m3 -> g/cm3
inline double density_si_to_g_cm3(double rho_si)
{
    return rho_si / 1000.0;
}

// -----------------------------------------------------------------------------
// Epsilon (Franck 1990) as in calculateEpsilon, Case 2
// -----------------------------------------------------------------------------

inline double epsilon_franck(double T, double rho_g_cm3)
{
    // T in K (here)
    // rho_g_cm3 in g/cm3
    //
    // Excel code:
    //   rhostar   = (density * 0.055508) * omega^3 * Na
    //   mustarsq  = mu^2 / (k * (T+273.15) * omega^3)
    //   y         = (4*pi/9)*rhostar*mustarsq
    //   f1, f2, f3 as below
    //   epsilon   = ((3*y)/(1-f1*y))*(1 + (1-f1)*y + f2*y^2 + f3*y^3) + 1

    // In DEW, input temperature to Franck is Celsius; they use (T_C + 273.15).
    // Here T is K already, so we use TK = T directly.

    const double TK = T;

    // density [mol/cm3]
    const double rho_mol_cm3 = rho_g_cm3 * 0.055508; // 1/18.01528

    const double cc      = pow(omega, 3.0) * Na;
    const double rhostar = rho_mol_cm3 * cc;

    const double mustarsq =
        (mu * mu) / (kB * TK * pow(omega, 3.0));

    const double y  = (4.0 * pi / 9.0) * rhostar * mustarsq;
    const double f1 = 0.4341 * pow(rhostar, 2.0);
    const double f2 = -(0.05 + 0.75 * pow(rhostar, 3.0));
    const double f3 = -0.026 * pow(rhostar, 2.0)
                      + 0.173 * pow(rhostar, 4.0);

    const double one_minus_f1y = 1.0 - f1 * y;

    const double term =
        1.0
        + (1.0 - f1) * y
        + f2 * y * y
        + f3 * y * y * y;

    const double eps =
        ((3.0 * y) / one_minus_f1y) * term + 1.0;

    return eps;
}

// -----------------------------------------------------------------------------
// (d epsilon / d rho)_T for Franck (calculate_depsdrho, Case 2)
// -----------------------------------------------------------------------------

inline double depsdrho_franck(double T, double rho_g_cm3)
{
    // Direct translation of DEW's calculate_depsdrho, Case 2 (Franck),
    // including the 0.05508 factor as in the provided code.

    const double TK = T;

    // density in mol/cm3
    double density_mol_cm3 = rho_g_cm3 * 0.055508;

    const double cc      = pow(omega, 3.0) * Na;
    const double rhostar = density_mol_cm3 * cc;

    const double mustarsq =
        (mu * mu) / (kB * TK * pow(omega, 3.0));

    const double y  = (4.0 * pi / 9.0) * rhostar * mustarsq;
    const double f1 = 0.4341 * pow(rhostar, 2.0);
    const double f2 = -(0.05 + 0.75 * pow(rhostar, 3.0));
    const double f3 = -0.026 * pow(rhostar, 2.0)
                      + 0.173 * pow(rhostar, 4.0);

    const double dydrho =
        (4.0 * pi / 9.0) * mustarsq * cc;

    const double df1drho =
        2.0 * 0.4341 * pow(cc, 2.0) * density_mol_cm3;

    const double df2drho =
        -3.0 * 0.75 * pow(cc, 3.0)
        * pow(density_mol_cm3, 2.0);

    const double df3drho =
        -2.0 * 0.026 * pow(cc, 2.0) * density_mol_cm3
        + 4.0 * 0.173 * pow(cc, 4.0)
          * pow(density_mol_cm3, 3.0);

    const double eps = epsilon_franck(T, rho_g_cm3);

    const double one_minus_f1y = 1.0 - f1 * y;

    // This is exactly the DEW expression:
    const double term1 =
        ((dydrho + y * y * df1drho) / one_minus_f1y)
        * (eps - 1.0) / y;

    const double term2 =
        (3.0 * y / one_minus_f1y) *
        (
            -df1drho * y
            + df2drho * y * y
            + df3drho * y * y * y
            + (1.0 - f1 + 2.0 * f2 * y + 3.0 * f3 * y * y) * dydrho
        );

    // Note: DEW uses 0.05508 (slight truncation of 0.055508).
    const double depsdrho =
        0.05508 * (term1 + term2);

    // This is ∂ε/∂ρ where ρ is in g/cm3.
    return depsdrho;
}

} // namespace

// -----------------------------------------------------------------------------
// Public interface
// -----------------------------------------------------------------------------

auto waterElectroPropsFranck1990(real T,
                                 real P,
                                 const WaterThermoProps& wt)
    -> WaterElectroProps
{
    WaterElectroProps we;

    // Convert density to DEW convention (g/cm3)
    const double rho_si     = wt.D;                  // kg/m3
    const double rho_g_cm3  = density_si_to_g_cm3(rho_si);

    // Dielectric constant
    const double eps = epsilon_franck(T, rho_g_cm3);
    we.epsilon = eps;

    // Derivative of epsilon wrt density (rho in g/cm3)
    const double deps_drho_g = depsdrho_franck(T, rho_g_cm3);

    // Map drho/dP from SI to g/cm3 per Pa:
    //
    //   rho_g = rho_si / 1000
    //   => d rho_g / dP = (1/1000) * wt.DP
    //
    const double drho_g_dP_SI = wt.DP / 1000.0;

    // epsilon_P (per Pa), using chain rule:
    //   epsilon_P = (dε/dρ_g) * (dρ_g/dP_SI)
    we.epsilonP = deps_drho_g * drho_g_dP_SI;

    // We do NOT fabricate epsilonT, epsilonTT, etc.,
    // because Franck branch in DEW only defines dε/dρ and
    // uses it to form Q. Set them to zero for clarity.
    we.epsilonT  = 0.0;
    we.epsilonTT = 0.0;
    we.epsilonTP = 0.0;
    we.epsilonPP = 0.0; // Only epsilonP (first derivative) is defined above.

    // Born-style coefficients, matched to DEW logic:
    //
    // Z = -1/ε
    // Q = (1/ε²) * dε/dP
    // Other second-order terms left zero here; they are not
    // defined in the Franck DEW implementation.

    we.bornZ = -1.0 / eps;

    we.bornQ = (eps != 0.0)
               ? (we.epsilonP / (eps * eps))
               : real(0.0);

    we.bornY = 0.0;
    we.bornX = 0.0;
    we.bornN = 0.0;
    we.bornU = 0.0;

    return we;
}

} // namespace Reaktoro
