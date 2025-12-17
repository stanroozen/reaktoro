#pragma once

// Fernandez et al. (1997) dielectric model for water
// --------------------------------------------------
// Exact translation of the `equation = 3` branch in the DEW Excel/VBA:
//
//   - calculateEpsilon (Case 3, Fernandez 1997)
//   - calculate_depsdrho (Case 3, Fernandez 1997)
//
// Excel/DEW conventions:
//   - Input density: g/cm^3
//   - Internal: density_molm3 in mol/m^3
//   - T in °C, then TK = T + 273.15
//
// Here we:
//   - Take T in K, P in Pa.
//   - Use WaterThermoProps (wt) with density in kg/m^3 & derivatives.
//   - Convert to g/cm^3 and mol/m^3 exactly as in DEW.
//   - Return epsilon, epsilonP, bornZ, bornQ consistent with DEW’s use
//     (via deps/drho and drho/dP). Other derivatives are set to 0.0
//     to avoid introducing non-DEW behavior.
//

#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Extensions/DEW/WaterElectroProps.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>

namespace Reaktoro {

auto waterElectroPropsFernandez1997(real T,
                                    real P,
                                    const WaterThermoProps& wt)
    -> WaterElectroProps;

} // namespace Reaktoro
