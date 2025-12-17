// WaterState.hpp
//
// High-level unified water state interface.
//
// Orchestrates the modular components so you can reproduce the full
// DEW/Helgeson-style behavior in a clean, extensible way:
//
//   - WaterThermoModel       : EOS selector (Wagner–Pruss, HGK, ZD05, ZD09, Psat poly)
//   - WaterDielectricModel   : epsilon models (Johnson–Norton, Franck, Fernandez, Power, Psat poly)
//   - WaterGibbsModel        : G(T,P) (Delaney–Helgeson poly, DEW integral, Psat poly)
//   - WaterSolventFunctionDEW: Shock/DEW g(P,T,ρ), dgdP (including Psat poly branch)
//   - WaterBornOmegaDEW      : Born omega(P,T), domega/dP using g, dgdP
//
// This header adds no new physics; it just wires the pieces together.

#pragma once

#include <Reaktoro/Common/Real.hpp>

#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>
#include <Reaktoro/Extensions/DEW/WaterElectroProps.hpp>

// Modular models
#include <Reaktoro/Extensions/DEW/WaterThermoModel.hpp>
#include <Reaktoro/Extensions/DEW/WaterDielectricModel.hpp>
#include <Reaktoro/Extensions/DEW/WaterGibbsModel.hpp>
#include <Reaktoro/Extensions/DEW/WaterSolventFunctionDEW.hpp>
#include <Reaktoro/Extensions/DEW/WaterBornOmegaDEW.hpp>

namespace Reaktoro {

/// Aggregated state for water at (T, P).
struct WaterState
{
    // Always computed
    WaterThermoProps  thermo;     ///< ρ, derivatives, etc. (SI)
    WaterElectroProps electro;    ///< ε, Born Z/Q/... (per chosen dielectric model)

    // Optional subsystems (enabled via WaterStateOptions)
    bool hasGibbs    = false;
    bool hasSolventG = false;
    bool hasOmega    = false;

    real gibbs       = 0.0;  ///< G(T,P) [J/mol]

    real g_solv      = 0.0;  ///< DEW solvent function g(P,T,ρ) [-]
    real dgdP        = 0.0;  ///< ∂g/∂P [1/Pa]

    real omega       = 0.0;  ///< Born omega(P,T) [J/mol]
    real domega_dP   = 0.0;  ///< ∂omega/∂P [J/mol/Pa]
};

/// Options controlling which models and subsystems are used.
struct WaterStateOptions
{
    // --- Core models (always used) ---

    /// EOS selection (Wagner–Pruss, HGK, ZD05, ZD09, Psat override, etc.)
    WaterThermoModelOptions     thermo;

    /// Dielectric selection (Johnson–Norton, Franck, Fernandez, Power, Psat override).
    WaterDielectricModelOptions dielectric;

    // --- Optional subsystems ---

    /// If true, compute Gibbs free energy with WaterGibbsModel.
    bool computeGibbs = false;

    /// Gibbs model options (Delaney–Helgeson, DEW integral, Psat poly).
    WaterGibbsModelOptions gibbs;

    /// If true, compute DEW solvent function g and dgdP.
    bool computeSolventG = false;

    /// Solvent function options:
    ///   - EOS choice (via density in WaterThermoProps)
    ///   - Psat flag for using Psat(T) polynomials when desired.
    WaterSolventFunctionOptions solvent;

    /// If true, compute Born omega and dω/dP.
    bool computeOmega = false;

    /// Born omega options:
    ///   - Uses g,dgdP from `solvent` options inside.
    ///   - Handles H+-like species, P cutoffs, etc.
    WaterBornOmegaOptions omega;
};

/// Compute a complete WaterState at (T, P) using the selected models.
///
/// @param T     Temperature [K]
/// @param P     Pressure [Pa]
/// @param opts  Model selections and feature flags
/// @return      Filled WaterState (unused parts flagged via has* booleans)
auto waterState(real T,
                real P,
                const WaterStateOptions& opts) -> WaterState;

} // namespace Reaktoro
