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
    real psatRelativeTolerance = 1e-3; // 0.1
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
