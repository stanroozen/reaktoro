#pragma once

// High-level selector for dielectric models of water.
//
// Wraps modular implementations:
//   - waterElectroPropsJohnsonNorton
//   - waterElectroPropsFranck1990
//   - waterElectroPropsFernandez1997
//   - waterElectroPropsPowerFunction
//
// and optionally applies DEW-style Psat polynomials (epsilon, Q, dgdP)
// near saturation, mimicking the Excel `Psat = True` behavior.
//
// Inputs:
//   - T  [K]
//   - P  [Pa]
//   - wt : WaterThermoProps (from your chosen EOS; see WaterThermoModel)
//   - options specifying which dielectric model + Psat handling.
//
// Output:
//   - WaterElectroProps consistent with the chosen model and (optionally)
//     DEW Psat polynomials, without polluting the EOS layer.

#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>
#include <Reaktoro/Extensions/DEW/WaterElectroProps.hpp>

namespace Reaktoro {

/// Dielectric models supported (aligned with DEW Excel `calculateEpsilon`).
enum class WaterDielectricModel
{
    JohnsonNorton1991,  ///< equation = 1 (SUPCRT / DEW default)
    Franck1990,         ///< equation = 2
    Fernandez1997,      ///< equation = 3
    PowerFunction       ///< equation = 4 (Sverjensky–Harrison)
};

/// Options for the unified dielectric selector.
struct WaterElectroModelOptions
{
    /// Which dielectric model to use.
    WaterDielectricModel model = WaterDielectricModel::JohnsonNorton1991;

    /// If true, use DEW Psat polynomials along saturation:
    ///  - epsilon(T) Psat polynomial
    ///  - Born Q(T) Psat polynomial
    ///  - (optionally) derive epsilonP from Q and epsilon.
    ///
    /// Semantics:
    ///  - Only applied when (T, P) is close to saturation according to
    ///    psatRelativeTolerance and using Wagner–Pruß Psat(T).
    bool usePsatPolynomials = false;

    /// Relative tolerance for |P - Psat(T)| / Psat(T) to trigger Psat override.
    real psatRelativeTolerance = 1e-3; // 0.1%
};

/// Unified dielectric properties of water at (T, P).
///
/// @param T     Temperature [K]
/// @param P     Pressure [Pa]
/// @param wt    WaterThermoProps from chosen EOS (density & derivatives)
/// @param opt   Dielectric & Psat configuration
/// @return      WaterElectroProps consistent with `opt`.
auto waterElectroPropsModel(real T,
                            real P,
                            const WaterThermoProps& wt,
                            const WaterElectroModelOptions& opt = {})
    -> WaterElectroProps;

} // namespace Reaktoro
