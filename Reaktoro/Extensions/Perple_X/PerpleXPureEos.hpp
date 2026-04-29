#pragma once

#include <array>

#include "PerpleXHybridEos.hpp"

namespace Reaktoro::PerpleX {

/// ============================================================================
/// Pure Fluid EOS Implementations for GFSM (Type 39)
/// ============================================================================
///
/// These functions provide alternative pure-species equations of state
/// used by the GFSM framework (Generic Fluid Solution Model, Perple_X Type 39).
///
/// IMPORTANT ARCHITECTURAL NOTE:
/// =============================
/// GFSM uses EXPLICIT SPECIATION SPACE:
/// - User provides all 12 mole fractions (Xn_CO2, Xn_H2O, ... Xn_HCl)
/// - Each species gets a pure-species EOS evaluation (independent)
/// - H2O, CO2, CH4 can use alternatives to MRK
/// - Other 9 species fixed to MRK
/// - No mixture law or coupling between species
///
/// GFSM Pure EOS Options:
/// - H2O (7 available): HSMRK, CORK, PSEOS, Haar, ZhangDuan05, ZhangDuan09
/// - CO2 (6 available): HSMRK, CORK, BRMRK, PSEOS, ZhangDuan09
/// - CH4 (3 available): HSMRK, ZhangDuan09
/// - Others (9 fixed): MRK (no switching)
///
/// Each function evaluates ln(fugacity) and partial molar volume for a
/// pure species at given P-T conditions.
///
/// How GFSM Uses These:
/// 1. User specifies: Xn_CO2, Xn_H2O, Xn_CH4, ... (all 12 species mole fractions)
/// 2. For each species i: Evaluate pure EOS at (P, T)
/// 3. Combine all 12 pure EOS evaluations directly (no mixing law)
/// 4. Result: Properties in explicit speciation space
///

struct PerpleXPureEosOptions
{
    /// Maximum iterations for root-finding in pure EOS (some use Newton-Raphson)
    int maxIter = 100;

    /// Maximum warnings to print during computation
    int maxWarn = 10;

    /// Convergence tolerance for nonlinear equations
    double tol = 1e-10;
};

/// ============================================================================
/// HSMRK - Hard-Sphere Modified Redlich-Kwong (Kerrick & Jacobs 1981)
/// ============================================================================
///
/// Pure-species ln(fugacity) with hard-sphere volume correction.
/// Callable for H2O, CO2, CH4 in GFSM framework.
///
/// Arguments:
///   - vol: [in/out] volume (cm3/mol), updated with solution
///   - species: species index (1=H2O, 2=CO2, 3=CH4, etc.)
///   - pressureBar: pressure in bar
///   - temperatureK: temperature in Kelvin
///   - options: convergence/iteration options
///
/// Returns: ln(fugacity coefficient)
///
double hsmrkf(double& vol, int species, double pressureBar, double temperatureK,
              const PerpleXPureEosOptions& options = {});

/// ============================================================================
/// CORK - Holland & Powell 1990 High-Pressure Water EOS
/// ============================================================================
///
/// Water-specific EOS optimized for high pressure and temperature.
/// Used for H2O in GFSM when high-P/T accuracy is needed.
///
void crkH2O(double pressureBar, double temperatureK, double& vol, double& lnfug);

/// ============================================================================
/// CORK - Holland & Powell 1990 High-Pressure CO2 EOS
/// ============================================================================
///
/// CO2-specific EOS optimized for high pressure and temperature.
/// Used for CO2 in GFSM when high-P/T accuracy is needed.
///
void crkCO2(double pressureBar, double temperatureK, double& vol, double& lnfug);

/// ============================================================================
/// PSEOS - Pitzer & Sterner Pseudo-Virial EOS
/// ============================================================================
///
/// Pseudo-virial expansion for H2O and CO2.
/// Callable for H2O (species=1) and CO2 (species=2) in GFSM.
///
double pseos(double& vol, int species, double pressureBar, double temperatureK,
             const PerpleXPureEosOptions& options = {});

/// ============================================================================
/// BRMRK - Bottinga & Richet 1981 CO2 EOS
/// ============================================================================
///
/// CO2-specific EOS from Bottinga & Richet 1981.
/// Used for CO2 in GFSM as alternative to MRK or CORK.
///
void brmrk(double pressureBar, double temperatureK, double& vol, double& lnfug,
           const PerpleXPureEosOptions& options = {});

/// ============================================================================
/// Haar - Haar et al. 1979 Water EOS
/// ============================================================================
///
/// High-accuracy water EOS based on Haar et al. (1979).
/// Most accurate for water in wide P-T range.
/// Used for H2O in GFSM when maximum accuracy needed.
///
void haar(double pressureBar, double temperatureK, double& vol, double& lnfug,
          const PerpleXPureEosOptions& options = {});

/// ============================================================================
/// ZhangDuan05 - Zhang & Duan 2005 Water Correlation
/// ============================================================================
///
/// Water fugacity correlation valid to high pressures.
/// Used for H2O in GFSM as alternative to MRK or CORK.
///
double zhdh2o(double& vol, double pressureBar, double temperatureK,
              const PerpleXPureEosOptions& options = {});

/// ============================================================================
/// ZhangDuan09 - Zhang & Duan 2009 Multi-Species Correlation
/// ============================================================================
///
/// Unified correlation for H2O, CO2, CH4, and other species.
/// Used in GFSM for H2O, CO2, CH4 as alternative to MRK.
///
double zd09pr(double& vol, int species, double pressureBar, double temperatureK,
              const PerpleXPureEosOptions& options = {});

/// ============================================================================
/// Build HybridEosOptions from Perple_X Pure EOS Implementations
/// ============================================================================
///
/// Factory function that creates a HybridEosOptions struct with all GFSM
/// pure EOS callbacks initialized to the Perple_X implementations above.
///
/// Use this to set up GFSM with all available pure EOS options.
///
HybridEosOptions makePerpleXHybridEosOptions(const PerpleXPureEosOptions& options = {});

/// Create HybridEosOptions with ZhangDuan09 for H2O, CO2, and CH4, matching
/// the recommended configuration for Perple_X COH-Fluid+ (iopt(25)=7,
/// iopt(26)=7, iopt(27)=7). All other callbacks are set to the standard
/// Perple_X pure EOS implementations.
HybridEosOptions makePerplexCOHFluidPlusEosOptions(const PerpleXPureEosOptions& options = {});

} // namespace Reaktoro::PerpleX
