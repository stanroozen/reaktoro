#pragma once

#include <array>

namespace Reaktoro::PerpleX {

/// MRK a/b parameters (cm3/mol; R = 83.1441) for Perple_X hybrid fluids.
struct MrkParameters
{
    /// a(i) coefficients, 1-based indexing (index 0 unused).
    std::array<double, 19> a{};

    /// b(i) coefficients, 1-based indexing (index 0 unused).
    std::array<double, 19> b{};
};

/// Compute MRK parameters at temperature T (K) using Perple_X rkparm logic.
MrkParameters mrkParameters(double temperatureK);

} // namespace Reaktoro::PerpleX
