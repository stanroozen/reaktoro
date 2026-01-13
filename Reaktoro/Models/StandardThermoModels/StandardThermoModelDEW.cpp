// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "StandardThermoModelDEW.hpp"

// C++ includes
#include <cmath>
using std::log;

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Extensions/DEW/WaterState.hpp>
#include <Reaktoro/Extensions/DEW/WaterBornOmegaDEW.hpp>
#include <Reaktoro/Serialization/Models/StandardThermoModels.hpp>

namespace Reaktoro {
namespace {

/// The reference temperature assumed in the HKF/DEW equations of state (in units of K)
const auto Tr = 298.15;

/// The reference pressure assumed in the HKF/DEW equations of state (in units of Pa)
const auto Pr = 1.0e+05;

/// The reference Born coefficient Z at Tr, Pr (dimensionless)
/// NOTE: DEW uses different water models, so Zr/Yr might differ from HKF
/// For now, using HKF reference values - may need adjustment for DEW
const auto Zr = -1.278055636e-02;

/// The reference Born coefficient Y at Tr, Pr (dimensionless)
const auto Yr = -5.795424563e-05;

/// The constant characteristic θ of the solvent (in units of K)
const auto theta = 228.0;

/// The constant characteristic ψ of the solvent (in units of Pa)
const auto psi = 2600.0e+05;

/// Configure WaterStateOptions for DEW thermo calculations
auto configureWaterStateOptions(const WaterModelOptions& waterOpts) -> WaterStateOptions
{
    WaterStateOptions opts;

    // Configure thermo model (EOS)
    opts.thermo.eosModel = waterOpts.eosModel;
    opts.thermo.usePsatPolynomials = waterOpts.usePsatPolynomials;
    opts.thermo.psatRelativeTolerance = waterOpts.psatRelTol;
    opts.thermo.densityTolerance = waterOpts.densityTolerance;

    // Configure dielectric model
    opts.dielectric.primary = static_cast<WaterDielectricPrimaryModel>(waterOpts.dielectricModel);
    if (waterOpts.usePsatPolynomials)
        opts.dielectric.psatMode = WaterDielectricPsatMode::UsePsatWhenNear;
    else
        opts.dielectric.psatMode = WaterDielectricPsatMode::None;
    opts.dielectric.psatRelativeTolerance = waterOpts.psatRelTol;

    // Configure Gibbs calculation (always required for species thermodynamics)
    opts.computeGibbs = true;
    opts.gibbs.model = waterOpts.gibbsModel;
    opts.gibbs.thermo = opts.thermo;  // Use same EOS for Gibbs integral
    // Use high-precision integration (5000 steps) by default
    opts.gibbs.integrationSteps = 5000;
    opts.gibbs.useExcelIntegration = false;
    opts.gibbs.densityTolerance = waterOpts.densityTolerance;

    // Enable solvent function (for omega calculation)
    opts.computeSolventG = true;
    opts.solvent.Psat = waterOpts.usePsatPolynomials;

    // Enable Born omega calculation if requested
    if (waterOpts.bornModel != WaterBornModel::None)
    {
        opts.computeOmega = true;
        opts.omega.solvent = opts.solvent;
    }

    return opts;
}} // namespace

auto StandardThermoModelDEW(const StandardThermoModelParamsDEW& params) -> StandardThermoModel
{
    auto evalfn = [=](StandardThermoProps& props, real T, real P)
    {
        auto& [G0, H0, V0, Cp0, VT0, VP0] = props;
        const auto& [Gf, Hf, Sr, a1, a2, a3, a4, c1, c2, wr, charge, Tmax, waterOpts] = params;

        // Configure and compute DEW water state
        const auto wsOpts = configureWaterStateOptions(waterOpts);
        const auto ws = waterState(T, P, wsOpts);

        // Extract water thermo properties
        const auto& wt = ws.thermo;
        const auto& we = ws.electro;

        // Born omega values (using DEW models if enabled)
        real w = 0.0;
        real wP = 0.0;
        real wT = 0.0;
        real wTP = 0.0;
        real wTT = 0.0;
        real wPP = 0.0;

        if (waterOpts.bornModel != WaterBornModel::None)
        {
            // Compute DEW Born omega and derivatives for ALL species (charged and neutral)
            // Neutral species have constant omega = wref (polarization/quadrupole)
            // Charged species have pressure-dependent omega from Born theory
            WaterBornOmegaOptions omegaOpts;  // Use default options
            w = waterBornOmegaDEW(T, P, wt, wr, charge, omegaOpts);
            wP = waterBornDOmegaDP_DEW(T, P, wt, wr, charge, omegaOpts);

            // For temperature derivatives, we'd need to compute at T±ε
            // Simplified approach: use Born function derivatives from WaterElectroProps
            // These are computed via the chosen dielectric model
            const auto& Z = we.bornZ;
            const auto& Y = we.bornY;
            const auto& Q = we.bornQ;
            const auto& U = we.bornU;
            const auto& N = we.bornN;
            const auto& X = we.bornX;

            // Approximate temperature derivatives using Born functions
            // (Full DEW would compute omega(T±ε, P) numerically)
            wT = -w * Y / Z;  // Simplified approximation
            wTP = 0.0;  // Would need numerical differentiation
            wTT = -w * X / Z;  // Simplified approximation
            wPP = 0.0;  // Would need numerical differentiation
        }

        // HKF/DEW non-solvation contributions
        const auto Tth = T - theta;
        const auto Tth2 = Tth * Tth;
        const auto Tth3 = Tth * Tth2;
        const auto psiP = psi + P;
        const auto psiPr = psi + Pr;

        // Standard molar volume V0
        V0 = a1 + a2/psiP + (a3 + a4/psiP)/Tth - w*we.bornQ - (we.bornZ + 1)*wP;

        // Temperature derivative of V0
        VT0 = -(a3 + a4/psiP)/(Tth*Tth) - wT*we.bornQ - w*we.bornU - we.bornY*wP - (we.bornZ + 1)*wTP;

        // Pressure derivative of V0
        VP0 = -a2/(psiP*psiP) - a4/(psiP*psiP*Tth) - wP*we.bornQ - w*we.bornN - we.bornQ*wP - (we.bornZ + 1)*wPP;

        // Standard molar Gibbs energy G0
        G0 = Gf - Sr*(T - Tr) - c1*(T*log(T/Tr) - T + Tr)
            + a1*(P - Pr) + a2*log(psiP/psiPr)
            - c2*((1.0/(T - theta) - 1.0/(Tr - theta))*(theta - T)/theta
            - T/(theta*theta)*log(Tr/T * (T - theta)/(Tr - theta)))
            + 1.0/Tth*(a3*(P - Pr) + a4*log(psiP/psiPr))
            - w*(we.bornZ + 1) + wr*(Zr + 1) + wr*Yr*(T - Tr);

        // Standard molar enthalpy H0
        H0 = Hf + c1*(T - Tr) - c2*(1.0/Tth - 1.0/(Tr - theta))
            + a1*(P - Pr) + a2*log(psiP/psiPr)
            + (2.0*T - theta)/Tth2*(a3*(P - Pr) + a4*log(psiP/psiPr))
            - w*(we.bornZ + 1) + w*T*we.bornY + T*(we.bornZ + 1)*wT + wr*(Zr + 1) - wr*Tr*Yr;

        // Standard molar heat capacity Cp0
        Cp0 = c1 + c2/Tth2 - 2.0*T/Tth3*(a3*(P - Pr) + a4*log(psiP/psiPr))
            + w*T*we.bornX + 2.0*T*we.bornY*wT + T*(we.bornZ + 1.0)*wTT;

        // Note: Entropy S0 can be derived from S0 = (H0 - G0)/T if needed
    };

    Data paramsdata;
    paramsdata["DEW"] = params;

    return StandardThermoModel(evalfn, paramsdata);
}

} // namespace Reaktoro
