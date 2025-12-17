#include <Reaktoro/Extensions/DEW/WaterGibbsModel.hpp>

#include <cmath>

#include <Reaktoro/Extensions/DEW/WaterPsatPolynomialsDEW.hpp>
#include <Reaktoro/Extensions/DEW/WaterHelmholtzPropsWagnerPruss.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoModel.hpp>

namespace Reaktoro {
namespace {

using std::abs;
using std::pow;

//----------------------------------------------------------------------------//
// Helpers
//----------------------------------------------------------------------------//

inline double toCelsius(double T_K)
{
    return T_K - 273.15;
}

// Use Wagner–Pruß saturation pressure as the reference Psat(T).
inline double saturationPressure(double T_K)
{
    // Adjust to your actual function name if different.
    return waterSaturationPressureWagnerPruss(T_K); // [Pa]
}

inline bool isNearPsat(double T_K,
                       double P_Pa,
                       double relTol)
{
    if (!std::isfinite(T_K) || !std::isfinite(P_Pa) || relTol <= 0.0)
        return false;

    const double Psat = saturationPressure(T_K);
    if (!(Psat > 0.0))
        return false;

    const double diff = std::abs(P_Pa - Psat);
    return diff <= relTol * Psat;
}

//----------------------------------------------------------------------------//
// 1) Delaney & Helgeson (1978) polynomial (Excel equation = 1)
//----------------------------------------------------------------------------//
//
// Excel Psat=False, equation=1:
//
//   coeff(0..14) given,
//   G[cal/mol] = sum_{j=0..4} sum_{k=0..4-j} coeff[count] * T^j * P^k
//   with T in °C, P in bar.
//
// We port it exactly and convert to SI: J/mol.
//
double gibbsDelaneyHelgeson1978_J_per_mol(double T_K, double P_Pa)
{
    const double T_C   = toCelsius(T_K);
    const double P_bar = P_Pa * 1.0e-5;

    const double c[15] = {
        -56130.073,
         0.38101798,
        -2.1167697e-6,
         2.0266445e-11,
        -8.3225572e-17,
        -15.285559,
         1.375239e-4,
        -1.5586868e-9,
         6.6329577e-15,
        -0.026092451,
         3.5988857e-8,
        -2.7916588e-14,
         1.7140501e-5,
        -1.6860893e-11,
        -6.0126987e-9
    };

    double G_cal = 0.0;
    int idx = 0;

    for (int j = 0; j <= 4; ++j)
    {
        const double Tj = pow(T_C, j);
        for (int k = 0; k <= 4 - j; ++k)
        {
            const double Pk = pow(P_bar, k);
            G_cal += c[idx++] * Tj * Pk;
        }
    }

    // cal/mol -> J/mol
    return G_cal * 4.184;
}

//----------------------------------------------------------------------------//
// 2) DEW integral model (Excel equation = 2)
//----------------------------------------------------------------------------//
//
// Excel idea (simplified):
//   - Fit G(T) at 1 kbar as polynomial in T (100–1000 °C).
//   - For P > 1 kbar, integrate V dP using density from selected EOS.
//   - Work in bar, cm3, cal; we do the same physics in SI:
//
//   dG/dP = V_m
//   G(P)  = G(1 kbar) + ∫_{1kbar}^{P} V_m(T, P') dP'
//
// We mirror:
//   - Threshold at 1000 bar (1e8 Pa)
//   - Step size as in Excel: spacing_bar = max(20, (Pbar - 1000)/500)
//   - Use waterThermoPropsModel with provided thermo options.
//
double GAtOneKb_cal_per_mol(double T_C)
{
    // From Excel GAtOneKb polynomial in equation=2 branch.
    // Valid roughly 100–1000 °C.
    return
        2.6880734e-9   * pow(T_C, 4.0) +
        6.3163061e-7   * pow(T_C, 3.0) -
        1.9372355e-2   * (T_C * T_C)   -
        16.945093      * T_C          -
        55769.287;
}

double gibbsDewIntegral_J_per_mol(double T_K,
                                  double P_Pa,
                                  const WaterThermoModelOptions& thermo)
{
    const double T_C = toCelsius(T_K);
    const double P_bar = P_Pa * 1.0e-5;

    // Excel: if P < 1000 bar, return 0 (not defined in that formulation).
    if (P_bar < 1000.0)
        return 0.0;

    // Base G at 1000 bar from polynomial (cal/mol).
    const double G1k_cal = GAtOneKb_cal_per_mol(T_C);

    if (P_bar == 1000.0)
        return G1k_cal * 4.184; // J/mol

    // For P > 1000 bar: integrate V_m dP from 1000 bar to P.
    //
    // Excel:
    //   spacing = max(20 bar, (P - 1000)/500)
    //   integral += (18.01528 / rho / 41.84) * spacing
    //
    // which is exactly V_m (cm3/mol) * spacing(bar) * (0.1 J/cm3/bar)/(4.184 J/cal)
    // -> cal/mol. We'll do equivalent in SI directly:
    //
    //   V_m = M / rho   [m3/mol], M=18.01528e-3 kg/mol
    //   dG = V_m * dP   [J/mol]

    const double M = 18.01528e-3; // kg/mol

    // Step size in bar as in Excel:
    const double spacing_bar = std::max(20.0, (P_bar - 1000.0) / 500.0);
    const double spacing_Pa  = spacing_bar * 1.0e5;

    double G_int_J = 0.0;

    // Integrate from 1000 bar up to P_bar.
    // Excel VBA uses: For i = 1000 To pressure Step spacing
    // VBA For...To loops are INCLUSIVE of both start and end points.
    // For P=6000 with spacing=20: evaluates at 1000, 1020, 1040, ..., 5980, 6000 (251 iterations)
    // This is NOT a proper LEFT Riemann sum (which would exclude the endpoint),
    // but we must match Excel exactly. Use <= with tolerance for floating-point safety.
    for (double Pstep_bar = 1000.0; Pstep_bar <= P_bar + 1e-9; Pstep_bar += spacing_bar)
    {
        const double Pstep_Pa = Pstep_bar * 1.0e5;

        // Use chosen EOS to get density at (T, Pstep).
        const auto wt = waterThermoPropsModel(T_K, Pstep_Pa, thermo);

        if (wt.D <= 0.0)
            continue; // skip unphysical; Excel code would just accumulate less

        const double Vm = M / wt.D; // [m3/mol]
        G_int_J += Vm * spacing_Pa;
    }

    // Total Gibbs:
    const double G_total_J = G1k_cal * 4.184 + G_int_J;
    return G_total_J;
}

//----------------------------------------------------------------------------//
// 3) Optional Psat(T) override (Excel Psat=True)
//----------------------------------------------------------------------------//
//
// Excel has a Psat branch for calculateGibbsOfWater:
//   G(T) along Psat(T) as a dedicated polynomial in T_C.
//
// We already implemented that in WaterPsatPolynomialsDEW.
// Here we optionally override the chosen model near Psat.
//
double maybeOverrideWithPsatGibbs(double T_K,
                                  double P_Pa,
                                  const WaterGibbsModelOptions& opt,
                                  double currentG_J_per_mol)
{
    if (!opt.usePsatPolynomials)
        return currentG_J_per_mol;

    if (!isNearPsat(T_K, P_Pa, opt.psatRelativeTolerance))
        return currentG_J_per_mol;

    const double G_psat = waterPsatGibbsDEW(T_K); // [J/mol] from Module 2
    return G_psat;
}

} // namespace

//----------------------------------------------------------------------------//
// Public interface
//----------------------------------------------------------------------------//

auto waterGibbsModel(real T,
                     real P,
                     const WaterGibbsModelOptions& opt)
    -> real
{
    const double T_K = static_cast<double>(T);
    const double P_Pa = static_cast<double>(P);

    double G_J_per_mol = 0.0;

    switch (opt.model)
    {
        case WaterGibbsModel::DelaneyHelgeson1978:
        {
            G_J_per_mol = gibbsDelaneyHelgeson1978_J_per_mol(T_K, P_Pa);
            break;
        }

        case WaterGibbsModel::DewIntegral:
        {
            G_J_per_mol = gibbsDewIntegral_J_per_mol(T_K, P_Pa, opt.thermo);
            break;
        }
    }

    // Optional DEW Psat polynomial override
    G_J_per_mol = maybeOverrideWithPsatGibbs(T_K, P_Pa, opt, G_J_per_mol);

    return G_J_per_mol;
}

} // namespace Reaktoro
