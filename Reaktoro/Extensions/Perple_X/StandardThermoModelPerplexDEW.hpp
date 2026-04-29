#pragma once

// Reaktoro includes
#include <Reaktoro/Core/StandardThermoModel.hpp>

namespace Reaktoro {

/// DEW-compatible HKF parameter surface backed by Perple_X water/electrolyte internals.
///
/// The fields intentionally match StandardThermoModelParamsDEW so Python scripts
/// can switch backend with minimal/no parameter-object churn.
struct StandardThermoModelParamsPerplexDEW
{
    real Gf;     ///< The apparent standard Gibbs energy of formation at reference T and P (J/mol)
    real Hf;     ///< The apparent standard enthalpy of formation at reference T and P (J/mol)
    real Sr;     ///< The standard entropy at reference T and P (J/(mol·K))
    real a1;     ///< The HKF EoS coefficient a1 (J/(mol·Pa)) — same convention as StandardThermoModelParamsDEW
    real a2;     ///< The HKF EoS coefficient a2 (J/mol) — same convention as StandardThermoModelParamsDEW
    real a3;     ///< The HKF EoS coefficient a3 (J·K/(mol·Pa)) — same convention as StandardThermoModelParamsDEW
    real a4;     ///< The HKF EoS coefficient a4 (J·K/mol) — same convention as StandardThermoModelParamsDEW
    real c1;     ///< The HKF EoS coefficient c1 (J/(mol·K)) — same convention as StandardThermoModelParamsDEW
    real c2;     ///< The HKF EoS coefficient c2 (J·K/mol) — same convention as StandardThermoModelParamsDEW
    real wref;   ///< The Born coefficient omega at reference conditions (J/mol) — same convention as StandardThermoModelParamsDEW
    real charge; ///< The electrostatic charge of the aqueous species
    real Tmax;   ///< The maximum temperature for which the model is valid (K)
};

/// Return a Perple_X-backed HKF standard thermo model with DEW-compatible API.
auto StandardThermoModelPerplexDEW(
    const StandardThermoModelParamsPerplexDEW& params) -> StandardThermoModel;

} // namespace Reaktoro
