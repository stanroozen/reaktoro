#include <Reaktoro/Extensions/DEW/WaterBornOmegaDEW.hpp>

#include <cmath>

namespace Reaktoro {
namespace {

using std::abs;

// DEW constant eta in units of (Å · cal / mol)
constexpr double eta_cal_per_A = 166027.0;

// Wrapper to compute g using the DEW solvent function module.
inline double compute_g(real T,
                        real P,
                        const WaterThermoProps& wt,
                        const WaterBornOmegaOptions& opt)
{
    return waterSolventFunctionDEW(T, P, wt, opt.solvent);
}

// Wrapper to compute dgdP [1/Pa] using the DEW solvent function module.
inline double compute_dgdP(real T,
                           real P,
                           const WaterThermoProps& wt,
                           double g,
                           const WaterBornOmegaOptions& opt)
{
    return waterSolventFunctionDgdP_DEW(T, P, wt, g, opt.solvent);
}

} // namespace

//----------------------------------------------------------------------------//
// omega(P, T)  [J/mol]
//----------------------------------------------------------------------------//

auto waterBornOmegaDEW(real T,
                       real P,
                       const WaterThermoProps& wt,
                       real wref_Jmol,
                       real Z,
                       const WaterBornOmegaOptions& opt)
    -> real
{
    // Trivial / Excel-like bypass cases
    if (Z == 0.0 || opt.isHydrogenLike || P > opt.maxPressureForVariation)
        return wref_Jmol;

    // Convert wref from J/mol to cal/mol for DEW formula.
    const double wref_cal = wref_Jmol / 4.184;

    // Excel:
    //   reref = Z^2 / (wref/eta + Z/3.082)
    const double denom = (wref_cal / eta_cal_per_A) + Z / 3.082;
    if (denom == 0.0)
        return wref_Jmol; // safeguard

    const double reref_A = (Z * Z) / denom; // [Å]

    // Solvent function g(T,P,ρ)
    const double g = compute_g(T, P, wt, opt);

    // Electrostatic radius at (P,T)
    const double re_A = reref_A + abs(Z) * g;
    if (re_A <= 0.0)
        return wref_Jmol; // safeguard

    // Excel DEW:
    //   omega(cal/mol) = eta * (Z^2 / re - Z / (3.082 + g))
    const double omega_cal =
        eta_cal_per_A * ((Z * Z) / re_A - Z / (3.082 + g));

    // Convert back to J/mol
    const double omega_J = omega_cal * 4.184;

    return omega_J;
}

//----------------------------------------------------------------------------//
// dω/dP  [J/mol/Pa]
//----------------------------------------------------------------------------//

auto waterBornDOmegaDP_DEW(real T,
                           real P,
                           const WaterThermoProps& wt,
                           real wref_Jmol,
                           real Z,
                           const WaterBornOmegaOptions& opt)
    -> real
{
    // Same trivial / cutoff logic as omega
    if (Z == 0.0 || opt.isHydrogenLike || P > opt.maxPressureForVariation)
        return 0.0;

    const double wref_cal = wref_Jmol / 4.184;
    const double denom = (wref_cal / eta_cal_per_A) + Z / 3.082;
    if (denom == 0.0)
        return 0.0;

    const double reref_A = (Z * Z) / denom;

    // g and dgdP from solvent function module
    const double g    = compute_g(T, P, wt, opt);
    const double dgdP = compute_dgdP(T, P, wt, g, opt); // [1/Pa]

    const double re_A = reref_A + abs(Z) * g;
    if (re_A <= 0.0)
        return 0.0;

    // Excel (conceptually, in cal/mol/bar):
    //   dω/dP = -eta * ( |Z|^3 / re^2 - Z / (3.082 + g)^2 ) * dgdP
    //
    // Our dgdP is already in 1/Pa, so we stay in SI and convert only
    // the energy units (cal -> J).

    const double term =
          std::pow(autodiff::detail::abs(Z), 3.0) / (re_A * re_A)
        - Z / ((3.082 + g) * (3.082 + g));

    // First in cal/mol/Pa
    const double domega_dP_cal_per_Pa =
        -eta_cal_per_A * term * dgdP;

    // Convert to J/mol/Pa
    const double domega_dP_J_per_Pa =
        domega_dP_cal_per_Pa * 4.184;

    return domega_dP_J_per_Pa;
}

} // namespace Reaktoro
