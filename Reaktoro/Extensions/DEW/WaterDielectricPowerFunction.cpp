#include "WaterDielectricPowerFunction.hpp"

// C++ includes
#include <cmath>
using std::exp;
using std::pow;
using std::sqrt;

// Reaktoro includes
#include <Reaktoro/Extensions/DEW/WaterElectroProps.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>

namespace Reaktoro {
namespace {

// -----------------------------------------------------------------------------
// Parameters (exactly as in your DEW VBA Power Function branch)
// -----------------------------------------------------------------------------

constexpr double a1 = -1.57637700752506e-03;
constexpr double a2 =  6.81028783422197e-02;
constexpr double a3 =  0.754875480393944;

constexpr double b1 = -8.01665106535394e-05;
constexpr double b2 = -6.87161761831994e-02;
constexpr double b3 =  4.74797272182151;

// Convert density from kg/m^3 (SI) to g/cm^3 (DEW convention)
inline double density_si_to_g_cm3(double rho_si)
{
    // 1 kg/m^3 = 0.001 g/cm^3
    return rho_si / 1000.0;
}

// -----------------------------------------------------------------------------
// Power Function epsilon as in calculateEpsilon, Case 4
// rho_g_cm3 : density in g/cm^3
// T_K       : temperature in K
// -----------------------------------------------------------------------------

inline double epsilon_power(double T_K, double rho_g_cm3)
{
    // DEW uses T in Celsius in this branch:
    //   A = a1*T + a2*sqrt(T) + a3
    //   B = b1*T + b2*sqrt(T) + b3
    //   epsilon = exp(B) * rho^A
    //
    // There T is °C. Here T_K is K, so:
    const double T_C = T_K - 273.15;

    // Note: DEW implicitly assumes T_C > 0 for sqrt(T_C); that is also
    // the physically relevant range for this model’s calibration.
    const double sqrtT = (T_C > 0.0) ? sqrt(T_C) : 0.0;

    const double A = a1 * T_C + a2 * sqrtT + a3;
    const double B = b1 * T_C + b2 * sqrtT + b3;

    // epsilon(rho) = exp(B) * rho^A
    // rho is in g/cm^3
    if(rho_g_cm3 <= 0.0)
        return 1.0; // graceful fallback; unphysical density

    const double eps = exp(B) * pow(rho_g_cm3, A);
    return eps;
}

// -----------------------------------------------------------------------------
// (d epsilon / d rho)_T as in calculate_depsdrho, Case 4
// rho_g_cm3 : density in g/cm^3
// T_K       : temperature in K
// -----------------------------------------------------------------------------

inline double depsdrho_power(double T_K, double rho_g_cm3)
{
    // From DEW:
    //   A = a1*T + a2*sqrt(T) + a3
    //   B = b1*T + b2*sqrt(T) + b3
    //   dε/dρ = A * exp(B) * ρ^(A-1)
    //
    // with T in °C, ρ in g/cm^3.

    const double T_C = T_K - 273.15;
    const double sqrtT = (T_C > 0.0) ? sqrt(T_C) : 0.0;

    const double A = a1 * T_C + a2 * sqrtT + a3;
    const double B = b1 * T_C + b2 * sqrtT + b3;

    if(rho_g_cm3 <= 0.0)
        return 0.0;

    const double depsdrho =
        A * exp(B) * pow(rho_g_cm3, A - 1.0);

    return depsdrho;
}

} // namespace

// -----------------------------------------------------------------------------
// Public interface
// -----------------------------------------------------------------------------

auto waterElectroPropsPowerFunction(real T,
                                    real P,
                                    const WaterThermoProps& wt)
    -> WaterElectroProps
{
    WaterElectroProps we;

    // Convert density from EOS to DEW units
    const double rho_si    = wt.D;                 // kg/m^3
    const double rho_g_cm3 = density_si_to_g_cm3(rho_si);

    // Dielectric constant epsilon (Power Function)
    const double eps = epsilon_power(T, rho_g_cm3);
    we.epsilon = eps;

    // dε/dρ (ρ in g/cm^3)
    const double deps_drho_g = depsdrho_power(T, rho_g_cm3);

    // dρ_g/dP with P in Pa:
    //
    //   ρ_g = ρ_SI / 1000
    //   dρ_g/dP = (1/1000) * dρ_SI/dP
    //
    const double drho_g_dP = wt.DP / 1000.0;

    // ε_P via chain rule:
    //   ε_P = (dε/dρ_g) * (dρ_g/dP)
    we.epsilonP = deps_drho_g * drho_g_dP;

    // As with Franck/Fernandez modules, we only define what DEW defines.
    // Higher-order T/P derivatives are set to zero to avoid adding
    // non-DEW behavior.
    we.epsilonT  = 0.0;
    we.epsilonTT = 0.0;
    we.epsilonTP = 0.0;
    we.epsilonPP = 0.0;

    // Born-style coefficients in DEW spirit:
    //
    // Z = -1/ε
    // Q = (1/ε²) * ε_P
    if(eps != 0.0)
    {
        const double eps2 = eps * eps;
        we.bornZ = -1.0 / eps;
        we.bornQ = we.epsilonP / eps2;
    }
    else
    {
        we.bornZ = 0.0;
        we.bornQ = 0.0;
    }

    we.bornY = 0.0;
    we.bornX = 0.0;
    we.bornN = 0.0;
    we.bornU = 0.0;

    return we;
}

} // namespace Reaktoro
