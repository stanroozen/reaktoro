#pragma once

// Zhang & Duan (2009) pure water EOS in DEW-style implementation.
//
// This module is a direct C++ translation of the `equation = 2` branch
// in the DEW Excel/VBA code you provided, including:
//   - P(rho, T)
//   - rho(P, T) via bisection
//   - (d rho / dP)_T from the DEW analytic expression
//
// External interface (Reaktoro-style):
//   - Input:  T in K, P in Pa
//   - Output: WaterThermoProps
//       D   = density [kg/m3]
//       DP  = (d rho / dP)_T [kg/(m3·Pa)]
//     Other derivatives left 0.0 here to keep this file a literal translation
//     of the provided DEW formulas without inventing missing theory.
//
// Use this side-by-side with the ZD05 module and select via options/enum
// in your higher-level water model selector.

#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>

namespace Reaktoro {

struct WaterZhangDuan2009Options
{
    /// Use Psat polynomial branch (as in DEW Psat=True).
    /// DEW uses the same fitted Psat density polynomial; we expose the flag.
    bool usePsat = false;

    /// Bisection tolerance in pressure (|Pcalc - Ptarget|) [bar].
    double pressureToleranceBar = 0.01;

    /// Maximum number of bisection iterations.
    int maxIterations = 50;
};

auto waterThermoPropsZhangDuan2009(real T, real P,
                                   const WaterZhangDuan2009Options& opts = {}) -> WaterThermoProps;

} // namespace Reaktoro
