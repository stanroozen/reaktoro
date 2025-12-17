// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// High-level selector for pure water equations of state used in DEW-style
// and Reaktoro-style calculations.
//
// This module unifies:
//   - Wagner & Pruß (IAPWS-95) via WaterHelmholtzPropsWagnerPruss
//   - Haar-Gallagher-Kell (HGK) via WaterHelmholtzPropsHGK
//   - Zhang & Duan (2005)  (DEW equation = 1)
//   - Zhang & Duan (2009)  (DEW equation = 2)
// with an optional DEW-style Psat-density polynomial override along saturation.
//
// Inputs:
//   - T [K], P [Pa]
//   - WaterThermoModelOptions to choose EOS and Psat handling
//
// Output:
//   - WaterThermoProps (SI units), consistent with the chosen EOS and,
//     if requested, DEW's Psat density polynomial near saturation.

#pragma once

#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>
#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>

// Existing EOS / utilities
#include <Reaktoro/Water/WaterHelmholtzProps.hpp>
#include <Reaktoro/Water/WaterHelmholtzPropsWagnerPruss.hpp>
#include <Reaktoro/Water/WaterHelmholtzPropsHGK.hpp>
#include <Reaktoro/Water/WaterThermoPropsUtils.hpp>
#include <Reaktoro/Water/WaterUtils.hpp>

// DEW-style Zhang & Duan EOS modules (you already have / just added these)
#include <Reaktoro/Extensions/DEW/WaterEosZhangDuan2005.hpp>
#include <Reaktoro/Extensions/DEW/WaterEosZhangDuan2009.hpp>

namespace Reaktoro {

/// Options controlling how waterThermoPropsModel selects and augments EOS.
struct WaterThermoModelOptions
{
    /// Which EOS to use.
    WaterEosModel eosModel = WaterEosModel::WagnerPruss;

    /// Use the DEW Psat density polynomial along the saturation curve
    /// for the Zhang & Duan EOS (equation = 1 or 2).
    ///
    /// When enabled:
    ///   - If (T,P) is within psatRelativeTolerance of the saturation
    ///     pressure from Wagner-Pruss, the density D is overridden by
    ///     the DEW Psat polynomial (liquid branch).
    ///   - This reproduces DEW’s behavior for Psat=True.
    bool usePsatPolynomials = false;

    /// Relative tolerance for |P - Psat(T)| / Psat(T) to trigger the
    /// Psat polynomial override. Ignored if usePsatPolynomials == false.
    double psatRelativeTolerance = 1e-3; // 0.1%

    /// Options forwarded to the Zhang & Duan (2009) DEW-style module.
    WaterZhangDuan2009Options zhangDuan2009Options;
};

/// High-level water EOS wrapper:
///   - Selects the requested EOS.
///   - Computes WaterThermoProps in SI units.
///   - Optionally applies DEW Psat density polynomial near saturation
///     for Zhang & Duan models when requested.
///
/// @param T    Temperature [K]
/// @param P    Pressure [Pa]
/// @param opt  EOS and Psat options
/// @return     WaterThermoProps, fully consistent with the chosen model.
auto waterThermoPropsModel(real T,
                           real P,
                           const WaterThermoModelOptions& opt = {})
    -> WaterThermoProps;

} // namespace Reaktoro
