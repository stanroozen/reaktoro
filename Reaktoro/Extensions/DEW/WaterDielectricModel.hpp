#pragma once

// Unified dielectric model selector for water (DEW-style + modular).
//
// This wraps the individual dielectric implementations:
//
//   - Johnson & Norton (1991)      -> waterElectroPropsJohnsonNorton
//   - Franck et al. (1990)        -> waterElectroPropsFranck1990
//   - Fernandez et al. (1997)     -> waterElectroPropsFernandez1997
//   - Power Function (DEW fit)    -> waterElectroPropsPowerFunction
//
// and optionally applies DEW Psat polynomials near saturation:
//
//   - epsilon(T) along Psat
//   - Q(T) along Psat
//
// External units:
//   - T in K, P in Pa
//   - WaterThermoProps in SI
//
// This module is purely a router/combiner: it does NOT change the
// underlying model formulas, it just selects and (optionally) blends
// them in a DEW-consistent way.

#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Extensions/DEW/WaterElectroProps.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>

namespace Reaktoro {

/// Primary dielectric model choices (non-Psat branches).
enum class WaterDielectricPrimaryModel
{
    JohnsonNorton1991,
    Franck1990,
    Fernandez1997,
    PowerFunction
};

/// How to treat Psat-specific polynomials.
enum class WaterDielectricPsatMode
{
    /// Never use Psat polynomials. Always call the primary model.
    None,

    /// If (T, P) is near saturation, override epsilon (and Q if available)
    /// using DEW Psat(T) polynomials.
    UsePsatWhenNear,

    /// Always use Psat(T) polynomials (caller guarantees they are appropriate).
    ForcePsat
};

struct WaterDielectricModelOptions
{
    /// Which base model to use away from Psat.
    WaterDielectricPrimaryModel primary = WaterDielectricPrimaryModel::JohnsonNorton1991;

    /// Psat handling strategy.
    WaterDielectricPsatMode psatMode = WaterDielectricPsatMode::None;

    /// Relative tolerance for |P - Psat(T)| / Psat(T) when psatMode == UsePsatWhenNear.
    real psatRelativeTolerance = 1e-3; // 0.1%

    /// If true, when Psat polynomials are used, overwrite both epsilon and Born Q
    /// from DEW fits (if available). If false, only epsilon is overridden and
    /// other props are left as from the primary model.
    bool overrideQWithPsatFit = true;
};

/// High-level entry point:
///
/// Compute WaterElectroProps using the selected dielectric model and
/// optional DEW Psat behavior.
///
/// @param T     Temperature [K]
/// @param P     Pressure [Pa]
/// @param wt    WaterThermoProps at (T, P)
/// @param opts  Dielectric model + Psat handling options
/// @return      WaterElectroProps populated consistently
auto waterElectroPropsModel(real T,
                            real P,
                            const WaterThermoProps& wt,
                            const WaterDielectricModelOptions& opts)
    -> WaterElectroProps;

} // namespace Reaktoro
