#pragma once

// Unified Gibbs free energy model for pure water.
//
// Mirrors the DEW Excel/VBA logic via a clean modular interface:
//
//   - Delaney & Helgeson (1978) polynomial:
//       Excel: calculateGibbsOfWater(equation = 1)
//   - DEW-style integral over volume from 1 kbar (1000 bar):
//       Excel: calculateGibbsOfWater(equation = 2)
//       using a chosen density EOS.
//
// Additionally, can optionally:
//
//   - Use DEW Psat(T) Gibbs polynomial along saturation
//     (Excel: Psat = True branch).
//
// Inputs:
//   - T [K], P [Pa]
//   - WaterGibbsModelOptions: choice of formulation and EOS
//
// Output:
//   - Gibbs free energy of water G [J/mol]

#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoModel.hpp>
#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>

namespace Reaktoro {

/// Enumeration of available numerical integration methods for Gibbs volume integral.
enum class WaterIntegrationMethod
{
    Trapezoidal     = 0,  ///< Trapezoidal rule: O(h²) error
    Simpson         = 1,  ///< Simpson's 1/3 rule: O(h⁴) error
    GaussLegendre16 = 2,  ///< 16-point Gauss-Legendre quadrature: O(1/n³²) error
    AdaptiveSimpson = 3,  ///< Adaptive Simpson's rule with auto-subdivision
};

/// Options to control Gibbs calculation.
struct WaterGibbsModelOptions
{
    /// Which Gibbs formulation to use.
    WaterGibbsModel model = WaterGibbsModel::DelaneyHelgeson1978;

    /// EOS + Psat options used when Gibbs depends on density (DewIntegral).
    WaterThermoModelOptions thermo;

    /// If true, use DEW Psat(T) Gibbs polynomial when (T,P) is
    /// sufficiently close to saturation (Psat=True behavior).
    ///
    /// This is applied on top of the chosen model when:
    ///   |P - Psat(T)| / Psat(T) <= psatRelativeTolerance.
    bool usePsatPolynomials = false;

    /// Relative tolerance for Psat proximity.
    real psatRelativeTolerance = 1e-3;

    /// Numerical integration method for volume integral.
    /// Default: Trapezoidal (O(h²), good balance of speed/accuracy)
    /// Options:
    ///   - Trapezoidal:     Fast, O(h²) error, good for 5000+ steps
    ///   - Simpson:         O(h⁴), ~1.5x slower than trapezoidal, better accuracy
    ///   - GaussLegendre16: Very high accuracy O(1/n³²), fewer function evals
    ///   - AdaptiveSimpson: Auto-subdivides to reach target tolerance
    WaterIntegrationMethod integrationMethod = WaterIntegrationMethod::Trapezoidal;

    /// Integration steps for DewIntegral model (when integrating V dP).
    /// Meaning depends on integrationMethod:
    ///   - Trapezoidal/Simpson: Total number of intervals
    ///   - GaussLegendre16: Number of 16-node segments (16*N evals total)
    ///   - AdaptiveSimpson: Initial number of intervals before subdivision
    /// Default 5000 gives high precision (~0.01% error) for trapezoidal.
    int integrationSteps = 5000;

    /// If true, use Excel VBA's exact integration logic (adaptive step size).
    /// If false, use fixed step size based on integrationSteps.
    /// Note: Only applies when integrationMethod = Trapezoidal.
    bool useExcelIntegration = false;

    /// Adaptive integration tolerance [J/mol] for AdaptiveSimpson method.
    /// Controls how much error is acceptable when subdividing intervals.
    /// Smaller = more accurate but slower (more function evaluations).
    /// Default 0.1 J/mol gives ~0.001% relative accuracy.
    double adaptiveIntegrationTolerance = 0.1;

    /// Maximum subdivisions for adaptive integration (safety limit).
    /// Prevents infinite recursion if tolerance is too tight.
    /// Default 20 allows up to 2^20 = 1,048,576 subdivisions.
    int maxAdaptiveSubdivisions = 20;

    /// Density calculation tolerance [bar] for integration.
    /// Default 0.001 bar gives high accuracy.
    /// Only affects Zhang & Duan EOS during Gibbs integration.
    double densityTolerance = 0.001;
};

/// Compute the Gibbs free energy of pure water at (T, P).
///
/// @param T    Temperature [K]
/// @param P    Pressure [Pa]
/// @param opt  Gibbs model and EOS/Psat options
/// @return     G [J/mol]
auto waterGibbsModel(real T,
                     real P,
                     const WaterGibbsModelOptions& opt = {})
    -> real;

} // namespace Reaktoro
