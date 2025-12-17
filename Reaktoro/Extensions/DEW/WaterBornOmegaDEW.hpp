#pragma once

// DEW/Shock-style Born coefficient omega(P, T) and its pressure derivative.
//
// Modular C++ translation of the Excel/VBA functions:
//   - calculateOmega(P, T, density, name, wref, Z)
//   - calculate_domegadP(P, T, density, name, wref, Z, densityEquation, Psat)
//
// This interface:
//   - Uses g and dgdP from WaterSolventFunctionDEW (DEW solvent function).
//   - Takes wref in J/mol (the actual omega, not 1e-5-scaled).
//   - Returns omega [J/mol] and dω/dP [J/mol/Pa].
//   - Lets the caller control hydrogen-like / neutral behavior and
//     pressure cutoff via options.
//   - Assumes WaterThermoProps is already computed by your chosen EOS.

#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>
#include <Reaktoro/Extensions/DEW/WaterSolventFunctionDEW.hpp>

namespace Reaktoro {

struct WaterBornOmegaOptions
{
    /// Solvent function configuration (EOS choice, Psat behavior).
    WaterSolventFunctionOptions solvent;

    /// If true, treat as H+-like or neutral in the DEW sense:
    ///   - For Z = 0 or hydrogen-like species, omega ≈ wref and dω/dP ≈ 0.
    bool isHydrogenLike = false;

    /// Maximum pressure [Pa] up to which DEW omega(P,T) variation is applied.
    /// Above this, we fall back to wref, as in Excel (P > 6000 bar).
    real maxPressureForVariation = 6000.0e5; // 6000 bar in Pa
};

/// Born coefficient omega(P, T) in J/mol.
///
/// Inputs:
///   - T [K], P [Pa]
///   - wt   : WaterThermoProps at (T,P)
///   - wref : reference omega at STP [J/mol]
///   - Z    : ionic charge
///   - opt  : options (solvent g, hydrogen-like flag, cutoff)
///
/// Behavior:
///   - If Z == 0, isHydrogenLike, or P > maxPressureForVariation:
///       returns wref.
///   - Else:
///       uses DEW/Shock formula with g(T,P) from WaterSolventFunctionDEW.
auto waterBornOmegaDEW(real T,
                       real P,
                       const WaterThermoProps& wt,
                       real wref,
                       real Z,
                       const WaterBornOmegaOptions& opt = {})
    -> real;

/// Pressure derivative dω/dP in J/mol/Pa.
///
/// Same inputs and logic as waterBornOmegaDEW. In the trivial cases
/// (Z == 0, hydrogen-like, or P above cutoff) returns 0.
auto waterBornDOmegaDP_DEW(real T,
                           real P,
                           const WaterThermoProps& wt,
                           real wref,
                           real Z,
                           const WaterBornOmegaOptions& opt = {})
    -> real;

} // namespace Reaktoro
