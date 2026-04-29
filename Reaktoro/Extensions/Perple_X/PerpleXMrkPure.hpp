#pragma once

#include <array>
#include <vector>

#include "PerpleXMrkMixture.hpp"

namespace Reaktoro::PerpleX {

/// ============================================================================
/// MRK Pure-Species Properties
/// ============================================================================
///
/// These functions compute pure-species MRK properties, which are used as:
/// 1. Baseline for GFSM hybrid framework
/// 2. Component of MRK mixture models
///
/// The pure-species MRK properties are then optionally replaced by
/// alternative pure EOS (HSMRK, CORK, etc.) in the GFSM framework.
///

struct MrkPureResult
{
    /// Log fugacity coefficients for all 19 potential species
    std::array<double, 19> ln_f{};

    /// Fugacity coefficients (Perple_X g array)
    std::array<double, 19> g{};

    /// Partial molar volumes cm³/mol (Perple_X v array)
    std::array<double, 19> v{};

    /// Total molar volume (for single species)
    double vol = 0.0;
};

/// ============================================================================
/// mrkPure - Compute Pure-Species MRK Properties
/// ============================================================================
///
/// Evaluates MRK equation of state for pure species.
///
/// This computes the baseline that GFSM uses before hybrid substitution.
///
/// Corresponds to Perple_X mrkpur() subroutine (flib.f).
///
MrkPureResult mrkPure(const std::vector<int>& species,
                      double pressureBar,
                      double temperatureK);

/// ============================================================================
/// loMrkMix - Low-Temperature MRK Mixture (near H2O critical point)
/// ============================================================================
///
/// Special low-temperature variant of MRK mixture for near-critical water.
///
/// Used when standard MRK produces numerical issues near critical point.
/// Rarely used; included for completeness.
///
/// Corresponds to Perple_X lomrk() subroutine (flib.f).
///
MrkMixResult loMrkMix(const std::vector<int>& species,
                      const std::array<double, 19>& y,
                      double pressureBar,
                      double temperatureK,
                      const MrkMixOptions& options);

} // namespace Reaktoro::PerpleX
