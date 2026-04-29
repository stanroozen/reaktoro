#pragma once

#include "PerpleXElectrolyte.hpp"
#include <array>

namespace Reaktoro::PerpleX {

// HKF (Helgeson-Kirkham-Flowers) aqueous species parameters
struct HKFParams
{
    // Standard state properties
    double G0 = 0.0;           // thermo(1) - Standard Gibbs energy
    double S0 = 0.0;           // thermo(2) - Standard entropy
    double omega0 = 0.0;       // thermo(5) - Born coefficient at ref conditions
    double charge = 0.0;       // thermo(6) - Ionic charge

    // HKF equation-of-state parameters
    double a1 = 0.0;           // thermo(7)
    double a2 = 0.0;           // thermo(8)
    double a3 = 0.0;           // thermo(9)
    double a4 = 0.0;           // thermo(10)
    double c1 = 0.0;           // thermo(11)
    double c2 = 0.0;           // thermo(12)

    // Pre-computed coefficients (from Perple_X preprocessing)
    double b8 = 0.0;           // thermo(13)
    double b9 = 0.0;           // thermo(14)
    double b10 = 0.0;          // thermo(15)
    double b11 = 0.0;          // thermo(16)
    double b12 = 0.0;          // thermo(17)
    double b13 = 0.0;          // thermo(18)
    double bornRadius = 0.0;   // thermo(19) - Reference born radius
    double chargeSquared = 0.0; // thermo(20)
};

// HKF calculation state
struct HKFState
{
    double G = 0.0;            // Apparent Gibbs energy (J/mol)
    double omega = 0.0;        // Born coefficient at P-T
    double bornTerm = 0.0;     // omega*(1/ε - 1) - omega0/ε0
};

/// HKF constants (from Perple_X rlib.f)
namespace HKFConstants {
    constexpr double PSI = 2600.0;           // Pressure reference parameter (bar)
    constexpr double THETA = 228.0;          // Temperature reference parameter (K)
    constexpr double ETA = 694656.968;       // Born function parameter
    constexpr double EPSILON0 = 78.244;      // Reference dielectric constant at 25°C, 1 bar (SUPCRT92/DEW, J&N1991 @ 997 kg/m³)
    constexpr double TR = 298.15;            // Reference temperature (K)
    constexpr double PR = 1.0;               // Reference pressure (bar)
    constexpr double NEUTRAL_RADIUS = 3.082; // Neutral species effective radius (Å)
}

/// Calculate Born omega coefficient for aqueous species
/// @param charge Ionic charge
/// @param bornRadius Reference Born radius (Å)
/// @param omega0 Born coefficient at reference conditions
/// @param gf Shock g-function value (Å)
/// @return Born coefficient omega at P-T
double calculateBornOmega(double charge,
                          double bornRadius,
                          double omega0,
                          double gf);

/// Calculate apparent Gibbs energy for HKF aqueous species
/// @param params HKF parameters for the species
/// @param pressureBar Pressure in bar
/// @param temperatureK Temperature in Kelvin
/// @param epsilon Dielectric constant
/// @param gf Shock g-function (Å)
/// @return HKF state with G, omega, and Born term
HKFState computeHKFGibbs(const HKFParams& params,
                         double pressureBar,
                         double temperatureK,
                         double epsilon,
                         double gf);

/// Pre-compute HKF coefficients b8-b13 from raw parameters
/// This matches Perple_X preprocessing done when loading thermodynamic data
/// @param params HKF parameters (a1-a4, c1, c2, G0, S0, omega0)
/// @param Tr Reference temperature (K)
/// @param Pr Reference pressure (bar)
/// @return Updated HKF parameters with b8-b13 filled
HKFParams preprocessHKFParams(const HKFParams& params,
                              double Tr = HKFConstants::TR,
                              double Pr = HKFConstants::PR);

/// Calculate water density at given P-T conditions
/// Required for g-function calculation
/// @param pressureBar Pressure in bar
/// @param temperatureK Temperature in Kelvin
/// @return Density in g/cm³
double waterDensity(double pressureBar, double temperatureK);

/// Get HKF solvent properties for pure water
/// This is the slvnt0 routine from Perple_X
/// @param pressureBar Pressure in bar
/// @param temperatureK Temperature in Kelvin
/// @param waterVolume Molar volume of water (cm³/mol) - output
/// @return DielectricState with all solvent properties
DielectricState getWaterSolventState(double pressureBar,
                                     double temperatureK,
                                     double& waterVolume);

} // namespace Reaktoro::PerpleX
