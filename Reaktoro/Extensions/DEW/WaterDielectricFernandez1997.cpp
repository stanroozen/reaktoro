#include "WaterDielectricFernandez1997.hpp"

// C++ includes
#include <cmath>
using std::pow;
using std::sqrt;

// Reaktoro includes
#include <Reaktoro/Extensions/DEW/WaterElectroProps.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>

namespace Reaktoro {
namespace {

// -----------------------------------------------------------------------------
// Constants exactly as in DEW / Fernandez et al. (1997) branch
// -----------------------------------------------------------------------------

// Avogadro's number [1/mol]
constexpr double avogadro = 6.0221367e23;

// Dipole moment of water [C·m] (1.84 Debye)
constexpr double dipole = 6.138e-30;

// Permittivity of free space [C^2 J^-1 m^-1]
constexpr double epsilon0 = 8.8541878176204e-12;

// Boltzmann constant [J K^-1]
constexpr double boltzmann = 1.380658e-23;

// Mean molecular polarizability of water [C^2 J^-1 m^-2]
constexpr double alpha = 1.636e-40;

// Critical density [mol m^-3]
constexpr double density_c = 17873.728; // 322 / Mw, from Fernandez, via DEW

// Critical temperature [K]
constexpr double T_c = 647.096;

// Density unit conversion:
//   - Input from EOS: kg/m^3
//   - DEW/Fernandez expects: density in g/cm^3 for the outer function,
//     then converts to mol/m^3 as:
//       density_molm3 = density_g_cm3 * 0.055508 * 1e6
//   (0.055508 ≈ 1/18.01528)
inline double density_si_to_g_cm3(double rho_si)
{
    return rho_si / 1000.0;
}

// -----------------------------------------------------------------------------
// N_k, i_k, j_k arrays (from DEW/Excel Fernandez implementation)
// -----------------------------------------------------------------------------

// N_k(0..11)
const double N_k[12] =
{
    0.978224486826,
   -0.957771379375,
    0.237511794148,
    0.714692224396,
   -0.298217036956,
   -0.108863472196,
    0.0949327488264,
   -0.00980469816509,
    0.000016516763497,
    0.0000937359795772,
   -1.2317921872e-10,
    0.00196096504426
};

// i_k(0..10)
const double i_k[11] =
{
    1, 1, 1, 2, 3, 3, 4, 5, 6, 7, 10
};

// j_k(0..10)
const double j_k[11] =
{
    0.25, 1, 2.5, 1.5, 1.5, 2.5, 2, 2, 5, 0.5, 10
};

// -----------------------------------------------------------------------------
// epsilon(T, rho) - Fernandez et al. (1997), DEW Case 3
// rho_g_cm3: density in g/cm^3
// T: temperature in K
// -----------------------------------------------------------------------------

inline double epsilon_fernandez(double T, double rho_g_cm3)
{
    // Convert density to mol/m^3 (as in DEW):
    const double density_molm3 = rho_g_cm3 * 0.055508 * 1.0e6;
    const double T_K = T;

    // Compute g(T, rho)
    double g = 1.0;

    const double x = density_molm3 / density_c;
    const double Tratio = T_c / T_K;

    for(int ii = 0; ii <= 10; ++ii)
    {
        g += N_k[ii]
             * pow(x, i_k[ii])
             * pow(Tratio, j_k[ii]);
    }

    g += N_k[11]
         * x
         * pow((T_K / 228.0) - 1.0, -1.2);

    // A, B, C as in DEW:

    const double A =
        (avogadro * dipole * dipole * density_molm3 * g)
        / (epsilon0 * boltzmann * T_K);

    const double B =
        (avogadro * alpha * density_molm3)
        / (3.0 * epsilon0);

    const double C =
          9.0
        + 2.0 * A
        + 18.0 * B
        + A * A
        + 10.0 * A * B
        + 9.0 * B * B;

    const double eps =
        (1.0 + A + 5.0 * B + std::sqrt(C))
        / (4.0 - 4.0 * B);

    return eps;
}

// -----------------------------------------------------------------------------
// d epsilon / d rho (rho in g/cm^3), DEW Case 3 (Fernandez)
// -----------------------------------------------------------------------------

inline double depsdrho_fernandez(double T, double rho_g_cm3)
{
    // This is a direct transcription of DEW's calculate_depsdrho, Case 3.

    const double density_molm3 = rho_g_cm3 * 0.055508 * 1.0e6;
    const double T_K = T;

    const double x = density_molm3 / density_c;
    const double Tratio = T_c / T_K;

    // g:
    double g = 1.0;
    for(int ii = 0; ii <= 10; ++ii)
    {
        g += N_k[ii]
             * pow(x, i_k[ii])
             * pow(Tratio, j_k[ii]);
    }
    g += N_k[11]
         * x
         * pow((T_K / 228.0) - 1.0, -1.2);

    // dg/drho  (rho in g/cm^3 => via density_molm3):
    //
    // In DEW they differentiate w.r.t density_molm3 first, then convert to
    // per (g/cm^3) using 55508 = 0.055508 * 1e6.
    //
    // Here we compute dg/drho_molm3, but formula below already expresses in
    // terms of density_molm3 and density_c.

    double dgdrho_molm3 = 0.0;

    for(int ii = 0; ii <= 10; ++ii)
    {
        const double ik = i_k[ii];
        dgdrho_molm3 += ik * N_k[ii]
                        * pow(density_molm3, ik - 1.0)
                        / pow(density_c, ik)
                        * pow(Tratio, j_k[ii]);
    }

    dgdrho_molm3 += (N_k[11] / density_c)
                    * pow((T_K / 228.0) - 1.0, -1.2);

    // Now A, B, C as functions of density_molm3:

    const double A =
        (avogadro * dipole * dipole * density_molm3 * g)
        / (epsilon0 * boltzmann * T_K);

    const double B =
        (avogadro * alpha * density_molm3)
        / (3.0 * epsilon0);

    const double C =
          9.0
        + 2.0 * A
        + 18.0 * B
        + A * A
        + 10.0 * A * B
        + 9.0 * B * B;

    const double eps =
        (1.0 + A + 5.0 * B + std::sqrt(C))
        / (4.0 - 4.0 * B);

    // dA/drho_molm3 and dB/drho_molm3 and dC/drho_molm3, as in DEW:

    const double dAdrho_molm3 =
        A / density_molm3 + (A / g) * dgdrho_molm3;

    const double dBdrho_molm3 =
        B / density_molm3;

    const double dCdrho_molm3 =
          2.0 * dAdrho_molm3
        + 18.0 * dBdrho_molm3
        + 2.0 * A * dAdrho_molm3
        + 10.0 * (dAdrho_molm3 * B + A * dBdrho_molm3)
        + 18.0 * B * dBdrho_molm3;

    // Convert from derivative wrt density_molm3 to wrt rho_g_cm3:
    //
    // density_molm3 = rho_g_cm3 * 0.055508 * 1e6
    // => d/d(rho_g_cm3) = (0.055508 * 1e6) * d/d(density_molm3) = 55508 * d/d(density_molm3)
    //
    // DEW packs this factor into the final formula as 55508.

    const double factor = 55508.0;

    // epsilon = (1 + A + 5B + sqrt(C)) / (4 - 4B)
    // dε/dρ = factor * 1/(4-4B) * (4 dBdrho ε + dAdrho + 5 dBdrho + 0.5 C^(-1/2) dCdrho)
    //
    // where dAdrho, dBdrho, dCdrho here are w.r.t density_molm3.

    const double denom = (4.0 - 4.0 * B);
    const double sqrtC = std::sqrt(C);

    const double depsdrho =
        factor * (1.0 / denom) *
        (
            4.0 * dBdrho_molm3 * eps
          + dAdrho_molm3
          + 5.0 * dBdrho_molm3
          + 0.5 * dCdrho_molm3 / sqrtC
        );

    // This is dε/dρ with ρ in g/cm^3, exactly as DEW.
    return depsdrho;
}

} // namespace

// -----------------------------------------------------------------------------
// Public interface
// -----------------------------------------------------------------------------

auto waterElectroPropsFernandez1997(real T,
                                    real P,
                                    const WaterThermoProps& wt)
    -> WaterElectroProps
{
    WaterElectroProps we;

    // Map density from EOS: kg/m^3 -> g/cm^3
    const double rho_si     = wt.D;                 // kg/m^3
    const double rho_g_cm3  = density_si_to_g_cm3(rho_si);

    // Dielectric constant ε (Fernandez 1997)
    const double eps = epsilon_fernandez(T, rho_g_cm3);
    we.epsilon = eps;

    // dε/dρ (ρ in g/cm^3)
    const double deps_drho_g = depsdrho_fernandez(T, rho_g_cm3);

    // drho_g/dP (P in Pa):
    //   rho_g = rho_si / 1000
    //   => d rho_g / dP = (1/1000) * wt.DP
    const double drho_g_dP = wt.DP / 1000.0;

    // ε_P via chain rule: ε_P = (dε/dρ_g) * (dρ_g/dP)
    we.epsilonP = deps_drho_g * drho_g_dP;

    // As in the DEW usage pattern for Fernandez:
    // we do not fabricate ε_T, higher T/P derivatives, or mixed terms here.
    we.epsilonT  = 0.0;
    we.epsilonTT = 0.0;
    we.epsilonTP = 0.0;
    we.epsilonPP = 0.0; // we only defined first derivative ε_P explicitly.

    // Born-style coefficients, consistent with DEW-style usage:
    //
    // Z = -1/ε
    // Q = (1/ε²) * ε_P
    // (others zeroed unless you later extend consistently)

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
