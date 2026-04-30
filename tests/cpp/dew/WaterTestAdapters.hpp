#pragma once

#include <string>

// Adapter API: these functions must be implemented in WaterTestAdapters.cpp
// by calling your existing DEW classes / functions and doing any necessary
// unit conversions so they match Excel-based truth tables.
//
// All temperatures here are in °C and pressures in bar (to match the CSVs).

// --------------------------- Density & drho/dP ---------------------------

// ρ(T, P) from Zhang & Duan 2005 and 2009, in g/cm³
double dew_density_ZD2005(double T_C, double P_bar);
double dew_density_ZD2009(double T_C, double P_bar);

// Saturated liquid density ρ_l(T) along Psat(T) using DEW Psat polynomial, in g/cm³
double dew_density_psat(double T_C);

// dρ/dP(T, P) for each EOS, in (g/cm³)/bar
double dew_drhodP_ZD2005(double T_C, double P_bar);
double dew_drhodP_ZD2009(double T_C, double P_bar);

// ------------------------- Dielectric constant ε ------------------------

// Relative dielectric constant ε_r(T, P) for each DEW dielectric model.
// These return the dimensionless dielectric constant eps_r columns.
double dew_epsilon_JN1991(double T_C, double P_bar);
double dew_epsilon_Franck1990(double T_C, double P_bar);
double dew_epsilon_Fernandez1997(double T_C, double P_bar);
double dew_epsilon_Power(double T_C, double P_bar);

// ε_r(T) along Psat(T) using DEW Psat polynomial
double dew_epsilon_psat(double T_C);

// ------------------------ dε/dρ (depsdrho) tests ------------------------

// (dε/dρ)_T as in your Excel calculate_depsdrho for each dielectric model.
// Input density is in g/cm³, output units are those of your CSVs (per g/cm³).
double dew_depsdrho_JN1991(double T_C, double rho_g_cm3);
double dew_depsdrho_Franck1990(double T_C, double rho_g_cm3);
double dew_depsdrho_Fernandez1997(double T_C, double rho_g_cm3);
double dew_depsdrho_Power(double T_C, double rho_g_cm3);

// ------------------------- Solvent function g(T,P) ----------------------

/// DEW solvent function g(T,P) used in Ω and Gibbs calculations.
/// Returns g in the same units as the "g" column in truth_g.csv.
double dew_g_eq2(double T_C, double P_bar);

/// d(g)/dP at constant T using DEW "equation 2" branch.
/// Returns in 1/bar, matching truth_dgdP_eq2.csv.
double dew_dgdP_eq2(double T_C, double P_bar);

/// d(g)/dP along Psat(T). Returns in 1/bar.
double dew_dgdP_psat(double T_C);

// ---------------------------- Gibbs free energy -------------------------

/// Delaney & Helgeson (1978) Gibbs polynomial (equation=1).
/// Returns G in cal/mol, matching truth_G_DH1978.csv.
double dew_G_DH1978(double T_C, double P_bar);

/// DEW integral formulation for G (equation=2), in cal/mol.
double dew_G_integral(double T_C, double P_bar);

/// DEW integral formulation for G with high precision (5000 steps), in cal/mol.
double dew_G_integral_highprec(double T_C, double P_bar);

/// Psat(T) Gibbs polynomial, in cal/mol.
double dew_G_psat(double T_C);

// ---------------------- Species Born Omega & dOmega/dP ------------------

// omega(P,T) for a given aqueous species in cal/mol.
// Inputs:
//   - speciesName      : name string as in truth_Omega_AllSpecies.csv
//   - T_C, P_bar       : state
//   - rho_g_cm3        : density (for consistency with Excel truth)
double dew_omega_species(const std::string& speciesName,
                         double T_C, double P_bar,
                         double rho_g_cm3);

// d(omega)/dP for a given aqueous species in cal/mol/bar
double dew_domegadP_species(const std::string& speciesName,
                            double T_C, double P_bar,
                            double rho_g_cm3);

// ------------------------------ Born Q(T,P) -----------------------------

// Final Q regression: density from densEq1 (ZD2005) and epsilon from epsEq4
// (Power-law dielectric). Returns Q in 1/bar to match truth_Q_densEq1_epsEq4.csv.
double dew_Q_densEq1_epsEq4(double T_C, double P_bar);

