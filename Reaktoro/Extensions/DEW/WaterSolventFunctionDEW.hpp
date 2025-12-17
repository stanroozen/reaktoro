// WaterSolventFunctionDEW.hpp

#pragma once

#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>

namespace Reaktoro {

struct WaterSolventFunctionOptions
{
    /// If true, treat evaluation as along Psat(T).
    /// We use DEW Psat polynomials for density and dgdP.
    bool Psat = false;

    /// (Optional) reserved for future: which density equation (ZD05/ZD09) etc.
    /// Included only if you want to mirror the Excel interface more literally.
    int densityEquation = 0;
};

auto waterSolventFunctionDEW(real T,
                             real P,
                             const WaterThermoProps& wt,
                             const WaterSolventFunctionOptions& opt = {})
    -> real;

auto waterSolventFunctionDgdP_DEW(real T,
                                  real P,
                                  const WaterThermoProps& wt,
                                  real g,
                                  const WaterSolventFunctionOptions& opt = {})
    -> real;

} // namespace Reaktoro
