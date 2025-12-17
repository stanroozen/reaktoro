#pragma once

#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {

/// Saturated liquid water density along Psat(T) using the DEW Excel polynomial.
/// Input:
///   T [K]
/// Output:
///   rho_l [kg/m3]
auto waterPsatDensityDEW(real T) -> real;

/// Dielectric constant of water along Psat(T) using the DEW Excel polynomial.
/// Input:
///   T [K]
/// Output:
///   epsilon [-]
auto waterPsatEpsilonDEW(real T) -> real;

/// Gibbs free energy of water along Psat(T) using the DEW Excel polynomial
/// from `calculateGibbsOfWater` Psat=True branch.
/// Input:
///   T [K]
/// Output:
///   G [J/mol]
auto waterPsatGibbsDEW(real T) -> real;

/// Born coefficient Q along Psat(T) from the DEW Excel Psat branch of calculateQ.
/// Input:
///   T [K]
/// Output:
///   Q [1/Pa]
auto waterPsatBornQDEW(real T) -> real;

/// d(g)/dP along Psat(T) from DEW Excel Psat branch of calculate_dgdP.
/// Input:
///   T [K]
/// Output:
///   dgdP [1/Pa]
auto waterPsatDgdPDEW(real T) -> real;

} // namespace Reaktoro
