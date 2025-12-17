#pragma once

#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Water/WaterThermoProps.hpp>

namespace Reaktoro {

/// Zhang & Duan (2005) pure water equation of state as used in DEW.
/// This is a direct C++ translation of the DEW Excel/VBA implementation
/// for equation = 1 (Zhang & Duan 2005), including:
///  - P(ρ, T)
///  - ρ(P, T) via bisection
///  - (∂ρ/∂P)_T analytic, exactly as in the Excel code
///
/// Inputs:
///  - T: Temperature in K
///  - P: Pressure in Pa
///
/// Outputs (WaterThermoProps):
///  - D   : density in kg/m3
///  - DP  : (∂ρ/∂P)_T in kg/m3/Pa (from DEW drhodP, unit-converted)
///  - other fields are set to 0.0 here (no extra theory invented).
auto waterThermoPropsZhangDuan2005(real T, real P) -> WaterThermoProps;

} // namespace Reaktoro
