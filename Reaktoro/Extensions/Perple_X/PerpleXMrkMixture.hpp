#pragma once

#include <array>
#include <vector>

#include "PerpleXMrkParameters.hpp"

namespace Reaktoro::PerpleX {

/// ============================================================================
/// MRK Mixture: Mixing-Rule EOS Framework
/// ============================================================================
///
/// MRK (Modified Redlich-Kwong) Mixture provides the mixing-rule EOS used as
/// the baseline for GFSM calculations.
///
/// MRK Mixture Properties:
/// - EOS Type: Modified Redlich-Kwong with temperature-dependent coefficients
/// - All 12 Species: Can be included (CO2, CH4, H2S, SO2, H2, CO, N2, NH3, HF, C2H6, HCl, H2O)
/// - Uses: Base for GFSM (ifug=39)
/// - Species interactions: Defined by MRK mixing rules
/// - No pure EOS switching: Fixed MRK parameters for all species
///
/// Relationship to GFSM:
/// - MRK Mixture: Provides the baseline 12-species system
/// - GFSM framework: Takes MRK mixture results, then replaces H2O/CO2/CH4
///   pure EOS with optional hybrid choices (HSMRK, CORK, etc.) in speciation space
///

struct MrkRootState
{
    /// Previous mixture volume (for Newton-Raphson continuation)
    double vrt = 0.0;

    /// Previous number of roots (cubic equation)
    int irt = 0;

    /// If true, select root consistent with previous iteration
    bool sroot = false;

    /// Root selection flag (unused in some implementations)
    bool max = false;
};

struct MrkMixOptions
{
    /// Averaging scheme for interaction parameters
    /// 1 = geometric (default, recommended)
    /// 2 = arithmetic
    /// else = harmonic
    int iavg = 1;

    /// Floor for mole fractions in log expressions (avoid log(0))
    double minY = 1e-12;
};

struct Roots3Result
{
    std::array<double, 3> roots{};
    double vmin = 0.0;
    double vmax = 0.0;
    int iroots = 0;
    int ineg = 0;
    int ipos = 0;
};

struct MrkMixResult
{
    std::array<double, 19> ln_f{};
    std::array<double, 19> g{};
    std::array<double, 19> v{};
    double vol = 0.0;
    Roots3Result roots{};
};

/// Solve x^3 + a1 x^2 + a2 x + a3 = 0 (real roots), Perple_X roots3.
Roots3Result roots3(double a1, double a2, double a3);

/// Compute MRK mixture fugacity coefficients and pmvs (Perple_X mrkmix).
MrkMixResult mrkMix(const std::vector<int>& species,
                    const std::array<double, 19>& y,
                    double pressureBar,
                    double temperatureK,
                    const MrkMixOptions& options,
                    MrkRootState* rootState = nullptr);

/// Apply hybrid pure-species fugacity coefficients to selected species.
void applyHybridFugacity(std::array<double, 19>& g,
                         const std::vector<int>& hybridSpecies,
                         const std::array<double, 19>& gh);

} // namespace Reaktoro::PerpleX
