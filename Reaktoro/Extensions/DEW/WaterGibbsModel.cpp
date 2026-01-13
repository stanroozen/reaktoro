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
// 2) DEW integral model with multiple integration methods
//----------------------------------------------------------------------------//

// Helper: 16-point Gauss-Legendre quadrature nodes and weights
// For interval [-1, 1]. Transform to [a, b] using: x = (b-a)/2 * xi + (a+b)/2
struct GaussLegendre16
{
    static constexpr int N = 16;
    
    // Nodes (abscissae) for [-1, 1]
    static constexpr double nodes[16] = {
        -0.9894009349916500,  -0.9445750230732326,  -0.8656312023878921,  -0.7554044083550030,
        -0.6178762444026437,  -0.4545454545454545,  -0.2736629541353671,  -0.0914297953593871,
         0.0914297953593871,   0.2736629541353671,   0.4545454545454545,   0.6178762444026437,
         0.7554044083550030,   0.8656312023878921,   0.9445750230732326,   0.9894009349916500
    };
    
    // Weights for [-1, 1]
    static constexpr double weights[16] = {
        0.0271524594117540,  0.0622535239386479,  0.0951585116824927,  0.1246289712555339,
        0.1495959888165767,  0.1691423829631435,  0.1826034150449235,  0.1894506104550585,
        0.1894506104550585,  0.1826034150449235,  0.1691423829631435,  0.1495959888165767,
        0.1246289712555339,  0.0951585116824927,  0.0622535239386479,  0.0271524594117540
    };
};

//----------------------------------------------------------------------------//
// Simpson's Rule (O(h⁴))
//----------------------------------------------------------------------------//
// Simpson's 1/3 rule: ∫f(x)dx ≈ (h/3) * (f₀ + 4f₁ + 2f₂ + 4f₃ + ... + fₙ)
// Requires: even number of intervals
double simpsonRule(double T_K,
                   double P_start_Pa,
                   double P_end_Pa,
                   int nsteps,
                   const WaterThermoModelOptions& thermoWithTol,
                   const double M)
{
    // Ensure even number of steps
    if (nsteps % 2 != 0) nsteps++;

    const double h = (P_end_Pa - P_start_Pa) / nsteps;
    
    auto wt = waterThermoPropsModel(T_K, P_start_Pa, thermoWithTol);
    double Vm_left = (wt.D > 0.0) ? (M / wt.D) : 0.0;
    double sum = Vm_left;

    // Odd indices (multiplied by 4)
    for (int i = 1; i < nsteps; i += 2)
    {
        wt = waterThermoPropsModel(T_K, P_start_Pa + i * h, thermoWithTol);
        double Vm = (wt.D > 0.0) ? (M / wt.D) : 0.0;
        sum += 4.0 * Vm;
    }

    // Even indices (multiplied by 2)
    for (int i = 2; i < nsteps; i += 2)
    {
        wt = waterThermoPropsModel(T_K, P_start_Pa + i * h, thermoWithTol);
        double Vm = (wt.D > 0.0) ? (M / wt.D) : 0.0;
        sum += 2.0 * Vm;
    }

    // Right endpoint
    wt = waterThermoPropsModel(T_K, P_end_Pa, thermoWithTol);
    double Vm_right = (wt.D > 0.0) ? (M / wt.D) : 0.0;
    sum += Vm_right;

    return (h / 3.0) * sum;
}

//----------------------------------------------------------------------------//
// 16-Point Gauss-Legendre Quadrature (O(1/n³²))
//----------------------------------------------------------------------------//
// High-order quadrature: divide [P_start, P_end] into nsegments 16-node segments
double gaussLegendre16(double T_K,
                       double P_start_Pa,
                       double P_end_Pa,
                       int nsegments,
                       const WaterThermoModelOptions& thermoWithTol,
                       const double M)
{
    const double segment_width = (P_end_Pa - P_start_Pa) / nsegments;
    double integral = 0.0;

    for (int seg = 0; seg < nsegments; ++seg)
    {
        double P_seg_start = P_start_Pa + seg * segment_width;
        double P_seg_end = P_seg_start + segment_width;

        // Transform from [-1, 1] to [P_seg_start, P_seg_end]
        double center = (P_seg_start + P_seg_end) / 2.0;
        double half_width = segment_width / 2.0;

        double seg_integral = 0.0;
        for (int i = 0; i < GaussLegendre16::N; ++i)
        {
            double P_node = center + half_width * GaussLegendre16::nodes[i];
            auto wt = waterThermoPropsModel(T_K, P_node, thermoWithTol);
            double Vm = (wt.D > 0.0) ? (M / wt.D) : 0.0;
            seg_integral += GaussLegendre16::weights[i] * Vm;
        }

        integral += half_width * seg_integral;
    }

    return integral;
}

//----------------------------------------------------------------------------//
// Adaptive Simpson's Rule with Automatic Subdivision
//----------------------------------------------------------------------------//
// Recursively subdivides intervals where |V(P_mid) - (V_L + V_R)/2| > tolerance
double adaptiveSimpsonsHelper(double T_K,
                             double P_L,
                             double P_R,
                             double tol,
                             int max_depth,
                             int current_depth,
                             const WaterThermoModelOptions& thermoWithTol,
                             const double M,
                             double Vm_L,
                             double Vm_R)
{
    if (current_depth > max_depth)
        return 0.0;  // Safety: stop recursing

    double P_mid = (P_L + P_R) / 2.0;
    auto wt_mid = waterThermoPropsModel(T_K, P_mid, thermoWithTol);
    double Vm_mid = (wt_mid.D > 0.0) ? (M / wt_mid.D) : 0.0;

    // Simpson's rule on [P_L, P_R]
    double h = (P_R - P_L) / 6.0;
    double S = h * (Vm_L + 4.0 * Vm_mid + Vm_R);

    // Check convergence criterion: |Vm_mid - (Vm_L + Vm_R)/2| vs tolerance
    double error_estimator = std::abs(Vm_mid - 0.5 * (Vm_L + Vm_R));
    
    if (error_estimator <= tol)
        return S;  // Converged, use this Simpson's estimate

    // Subdivide into [P_L, P_mid] and [P_mid, P_R]
    double P_left_mid = (P_L + P_mid) / 2.0;
    auto wt_left_mid = waterThermoPropsModel(T_K, P_left_mid, thermoWithTol);
    double Vm_left_mid = (wt_left_mid.D > 0.0) ? (M / wt_left_mid.D) : 0.0;

    double P_right_mid = (P_mid + P_R) / 2.0;
    auto wt_right_mid = waterThermoPropsModel(T_K, P_right_mid, thermoWithTol);
    double Vm_right_mid = (wt_right_mid.D > 0.0) ? (M / wt_right_mid.D) : 0.0;

    double left_integral = adaptiveSimpsonsHelper(T_K, P_L, P_mid, tol / 2.0, max_depth, current_depth + 1,
                                                  thermoWithTol, M, Vm_L, Vm_mid);
    double right_integral = adaptiveSimpsonsHelper(T_K, P_mid, P_R, tol / 2.0, max_depth, current_depth + 1,
                                                   thermoWithTol, M, Vm_mid, Vm_R);

    return left_integral + right_integral;
}

double adaptiveSimpson(double T_K,
                      double P_start_Pa,
                      double P_end_Pa,
                      double tol,
                      int max_depth,
                      const WaterThermoModelOptions& thermoWithTol,
                      const double M)
{
    auto wt_start = waterThermoPropsModel(T_K, P_start_Pa, thermoWithTol);
    double Vm_start = (wt_start.D > 0.0) ? (M / wt_start.D) : 0.0;

    auto wt_end = waterThermoPropsModel(T_K, P_end_Pa, thermoWithTol);
    double Vm_end = (wt_end.D > 0.0) ? (M / wt_end.D) : 0.0;

    return adaptiveSimpsonsHelper(T_K, P_start_Pa, P_end_Pa, tol, max_depth, 0,
                                 thermoWithTol, M, Vm_start, Vm_end);
}

//
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
                                  const WaterThermoModelOptions& thermo,
                                  const WaterGibbsModelOptions& opt)
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

    // Create a copy of thermo options with the density tolerance from Gibbs options
    WaterThermoModelOptions thermoWithTol = thermo;
    thermoWithTol.densityTolerance = opt.densityTolerance;

    double G_int_J = 0.0;

    if (opt.useExcelIntegration)
    {
        // Excel VBA compatibility mode: adaptive step size with 5000 steps max
        // Step size in bar as in Excel:
        const double spacing_bar = std::max(1.0, (P_bar - 1000.0) / 5000.0);
        const double spacing_Pa  = spacing_bar * 1.0e5;

        // Integrate from 1000 bar up to P_bar.
        // Excel VBA uses: For i = 1000 To pressure Step spacing
        // VBA For...To loops are INCLUSIVE of both start and end points.
        // For P=6000 with spacing=1: evaluates at 1000, 1001, 1002, ..., 5999, 6000 (5001 iterations)
        // This is NOT a proper LEFT Riemann sum (which would exclude the endpoint),
        // but we must match Excel exactly. Use <= with tolerance for floating-point safety.
        for (double Pstep_bar = 1000.0; Pstep_bar <= P_bar + 1e-9; Pstep_bar += spacing_bar)
        {
            const double Pstep_Pa = Pstep_bar * 1.0e5;

            // Use chosen EOS to get density at (T, Pstep).
            const auto wt = waterThermoPropsModel(T_K, Pstep_Pa, thermoWithTol);

            if (wt.D <= 0.0)
                continue; // skip unphysical; Excel code would just accumulate less

            const double Vm = M / wt.D; // [m3/mol]
            G_int_J += Vm * spacing_Pa;
        }
    }
    else
    {
        // High-precision mode: choose integration method
        const double P_start_Pa = 1000.0 * 1.0e5;

        switch (opt.integrationMethod)
        {
            case WaterIntegrationMethod::Trapezoidal:
            {
                // Fixed step trapezoidal rule: O(h²)
                const int nsteps = opt.integrationSteps;
                const double dP = (P_Pa - P_start_Pa) / nsteps;

                auto wt_prev = waterThermoPropsModel(T_K, P_start_Pa, thermoWithTol);
                real Vm_prev = (wt_prev.D > 0.0) ? (M / wt_prev.D) : real(0.0);

                for (int i = 1; i <= nsteps; ++i)
                {
                    const double Pstep_Pa = P_start_Pa + i * dP;
                    const auto wt = waterThermoPropsModel(T_K, Pstep_Pa, thermoWithTol);

                    if (wt.D <= 0.0)
                        continue;

                    const double Vm = M / wt.D;
                    G_int_J += 0.5 * (Vm_prev + Vm) * dP;
                    Vm_prev = Vm;
                }
                break;
            }

            case WaterIntegrationMethod::Simpson:
            {
                // Simpson's 1/3 rule: O(h⁴)
                G_int_J = simpsonRule(T_K, P_start_Pa, P_Pa, opt.integrationSteps, thermoWithTol, M);
                break;
            }

            case WaterIntegrationMethod::GaussLegendre16:
            {
                // 16-point Gauss-Legendre quadrature: O(1/n³²)
                // integrationSteps = number of 16-node segments
                int nsegments = std::max(1, opt.integrationSteps / 16);
                G_int_J = gaussLegendre16(T_K, P_start_Pa, P_Pa, nsegments, thermoWithTol, M);
                break;
            }

            case WaterIntegrationMethod::AdaptiveSimpson:
            {
                // Adaptive Simpson's with automatic subdivision
                // integrationSteps used as initial subdivision count
                G_int_J = adaptiveSimpson(T_K, P_start_Pa, P_Pa, 
                                         opt.adaptiveIntegrationTolerance,
                                         opt.maxAdaptiveSubdivisions,
                                         thermoWithTol, M);
                break;
            }
        }
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
            G_J_per_mol = gibbsDewIntegral_J_per_mol(T_K, P_Pa, opt.thermo, opt);
            break;
        }
    }

    // Optional DEW Psat polynomial override
    G_J_per_mol = maybeOverrideWithPsatGibbs(T_K, P_Pa, opt, G_J_per_mol);

    return G_J_per_mol;
}

} // namespace Reaktoro
