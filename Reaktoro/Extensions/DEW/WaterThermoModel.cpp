#include "WaterThermoModel.hpp"

#include <cmath>

namespace Reaktoro {
namespace {

using std::abs;
using std::pow;

//------------------------------------------------------------------------------
// 1) DEW Psat density polynomial (liquid) from Excel, Psat=True branch
//------------------------------------------------------------------------------
//
// Excel DEW:
//   calculateDensity(Psat=True) uses:
//     rho(T_C) [g/cm3] =
//       -1.01023381581205E-104 * T_C^40
//       -1.1368599785953E-27   * T_C^10
//       -2.11689207168779E-11  * T_C^4
//       +1.26878850169523E-08  * T_C^3
//       -4.92010672693621E-06  * T_C^2
//       -3.2666598612692E-05   * T_C
//       +1.00046144613017
//
// with T_C in °C. This is a *fit to saturated liquid density* used by DEW.
// We keep it bitwise identical here.

inline double dewPsatDensityLiquid_g_cm3_from_T(double T_K)
{
    const double T_C = T_K - 273.15;

    const double t2  = T_C * T_C;
    const double t3  = t2 * T_C;
    const double t4  = t2 * t2;
    const double t10 = pow(T_C, 10.0);
    const double t40 = pow(T_C, 40.0);

    const double rho =
        -1.01023381581205e-104 * t40 +
        -1.1368599785953e-27   * t10 +
        -2.11689207168779e-11  * t4  +
         1.26878850169523e-08  * t3  +
        -4.92010672693621e-06  * t2  +
        -3.2666598612692e-05   * T_C +
         1.00046144613017;

    return rho; // g/cm3
}

//------------------------------------------------------------------------------
// 2) Helper: check if (T,P) is near saturation (using Wagner–Pruß Psat(T))
//------------------------------------------------------------------------------

inline bool isNearPsatWagnerPruss(double T_K,
                                  double P_Pa,
                                  double relTol)
{
    if (!std::isfinite(T_K) || !std::isfinite(P_Pa) || relTol <= 0.0)
        return false;

    const double Psat = waterSaturationPressureWagnerPruss(T_K); // [Pa]

    if (!(Psat > 0.0))
        return false;

    const double diff = std::abs(P_Pa - Psat);
    return diff <= relTol * Psat;
}

//------------------------------------------------------------------------------
// 3) Apply DEW Psat polynomial override for Zhang & Duan models
//------------------------------------------------------------------------------

inline void maybeApplyDewPsatOverride(double T_K,
                                      double P_Pa,
                                      const WaterThermoModelOptions& opt,
                                      WaterThermoProps& wt)
{
    if (!opt.usePsatPolynomials)
        return;

    // Only meaningful for DEW-style Zhang & Duan branches
    if (opt.eosModel != WaterEosModel::ZhangDuan2005 &&
        opt.eosModel != WaterEosModel::ZhangDuan2009)
        return;

    if (!isNearPsatWagnerPruss(T_K, P_Pa, opt.psatRelativeTolerance))
        return;

    // Use DEW Psat polynomial for saturated liquid density:
    const double rho_g_cm3 = dewPsatDensityLiquid_g_cm3_from_T(T_K);

    // Map to SI:
    const double rho_kg_m3 = rho_g_cm3 * 1000.0;

    wt.D   = rho_kg_m3;

    // The original DEW Psat polynomial branch is a standalone fit; it does
    // not provide consistent analytic derivatives with respect to P or T.
    // In DEW usage, Psat=True is a special-path evaluation.
    //
    // To remain faithful and avoid inventing theory, we zero derivatives
    // here instead of mixing mismatched slopes from the bulk EOS.
    wt.DP  = 0.0;
    wt.DT  = 0.0;
    wt.DTT = 0.0;
    wt.DTP = 0.0;
    wt.DPP = 0.0;
}

} // namespace

//------------------------------------------------------------------------------
// Public interface
//------------------------------------------------------------------------------

auto waterThermoPropsModel(real T,
                           real P,
                           const WaterThermoModelOptions& opt)
    -> WaterThermoProps
{
    const double T_K = static_cast<double>(T);
    const double P_Pa = static_cast<double>(P);

    WaterThermoProps wt;

    switch (opt.eosModel)
    {
        case WaterEosModel::WagnerPruss:
        {
            // Use existing Helmholtz-based implementation
            const auto whp = waterHelmholtzPropsWagnerPruss(T_K, P_Pa);
            wt = waterThermoProps(T_K, P_Pa, whp);
            break;
        }

        case WaterEosModel::HGK:
        {
            const auto whp = waterHelmholtzPropsHGK(T_K, P_Pa);
            wt = waterThermoProps(T_K, P_Pa, whp);
            break;
        }

        case WaterEosModel::ZhangDuan2005:
        {
            // Exact DEW-style Zhang & Duan (2005) translation
            wt = waterThermoPropsZhangDuan2005(T_K, P_Pa);
            break;
        }

        case WaterEosModel::ZhangDuan2009:
        {
            // Exact DEW-style Zhang & Duan (2009) translation
            wt = waterThermoPropsZhangDuan2009(T_K, P_Pa,
                                               opt.zhangDuan2009Options);
            break;
        }
    }

    // Optional DEW Psat polynomial behavior (Psat=True analogue)
    maybeApplyDewPsatOverride(T_K, P_Pa, opt, wt);

    return wt;
}

} // namespace Reaktoro
