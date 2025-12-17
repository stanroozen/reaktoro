#pragma once

// Franck et al. (1990) dielectric model for water
// ------------------------------------------------
// Exact translation of the `equation = 2` (Franck) branch in your DEW Excel/VBA:
//
//   Uses Lennard-Jones-based reduced density rhostar and reduced dipole mustarsq
//   to compute the dielectric constant epsilon(ρ, T).
//
//   Also uses the DEW expression for (d epsilon / d rho)_T to obtain
//   Q = (1/epsilon^2) * (d epsilon / dP)_T via chain rule with
//   drho/dP from WaterThermoProps.
//
// Interface:
//   - T  [K]
//   - P  [Pa]
//   - wt : WaterThermoProps (density & derivatives from chosen water EOS)
//   - returns WaterElectroProps with epsilon, epsilonP, bornZ, bornQ.
//     All other fields are set to 0.0 to stay consistent with what the
//     DEW Franck implementation actually defines.
//

#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Extensions/DEW/WaterElectroProps.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>

namespace Reaktoro {

auto waterElectroPropsFranck1990(real T,
                                 real P,
                                 const WaterThermoProps& wt)
    -> WaterElectroProps;

} // namespace Reaktoro
