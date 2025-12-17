#pragma once

// Power Function dielectric model for water
// -----------------------------------------
// Exact translation of `equation = 4` branch in your DEW Excel/VBA:
//
//   epsilon = exp(B) * rho^A
//   A = a1*T_C + a2*sqrt(T_C) + a3
//   B = b1*T_C + b2*sqrt(T_C) + b3
//
//   (d epsilon / d rho)_T = A * exp(B) * rho^(A-1)
//
// where:
//   - T_C is temperature in °C
//   - rho is density in g/cm^3
//
// Here:
//   - Input T is in K, P in Pa (Reaktoro convention).
//   - WaterThermoProps provides density D in kg/m^3 and DP = ∂D/∂P in SI.
//   - We convert to g/cm^3 to preserve DEW’s implementation exactly.
//
// Output:
//   - WaterElectroProps:
//       epsilon   : Power Function epsilon
//       epsilonP  : via chain rule from dε/dρ and dρ/dP
//       bornZ     : -1/epsilon
//       bornQ     : epsilonP / epsilon²
//       others    : set to 0.0 to avoid introducing non-DEW behavior.
//

#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Water/WaterElectroProps.hpp>
#include <Reaktoro/Water/WaterThermoProps.hpp>

namespace Reaktoro {

auto waterElectroPropsPowerFunction(real T,
                                    real P,
                                    const WaterThermoProps& wt)
    -> WaterElectroProps;

} // namespace Reaktoro
