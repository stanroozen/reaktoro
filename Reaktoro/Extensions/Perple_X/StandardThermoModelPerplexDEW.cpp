#include "StandardThermoModelPerplexDEW.hpp"

// C++ includes
#include <algorithm>
#include <cmath>

// Reaktoro includes
#include <Reaktoro/Extensions/Perple_X/PerpleXHKF.hpp>

namespace Reaktoro {

namespace {

auto makePerpleXHKFParams(const StandardThermoModelParamsPerplexDEW& params) -> PerpleX::HKFParams
{
    // The PerplexHKF engine (ported from Perple_X rlib.f) uses pure SI-J units throughout,
    // with the exception that pressure is in bar rather than Pa.  This matches the Perple_X
    // datafile convention (e.g. DEW24HP633ver_elements.dat) where G0/S0/a2/a4/c1/c2/wref
    // are already in J and a1/a3 are in J/(mol·bar) and J·K/(mol·bar) respectively.
    //
    // The DEW database in Reaktoro (dew2024-aqueous.yaml) stores the same quantities but with
    // a1 in J/(mol·Pa) and a3 in J·K/(mol·Pa) — a factor of 1e5 smaller.  All other
    // parameters are identical to the Perple_X datafile (already in J).
    //
    // Conversion: multiply a1 and a3 by 1e5 (Pa→bar); everything else is unchanged.
    constexpr double bar = 1.0e5;   // Pa per bar

    PerpleX::HKFParams hkf;
    hkf.G0     = params.Gf;             // J/mol            (no conversion)
    hkf.S0     = params.Sr;             // J/(mol·K)         (no conversion)
    hkf.omega0 = params.wref;           // J/mol            (no conversion)
    hkf.charge = params.charge;
    hkf.a1     = params.a1 * bar;       // J/(mol·Pa) → J/(mol·bar)
    hkf.a2     = params.a2;             // J/mol            (no conversion)
    hkf.a3     = params.a3 * bar;       // J·K/(mol·Pa) → J·K/(mol·bar)
    hkf.a4     = params.a4;             // J·K/mol          (no conversion)
    hkf.c1     = params.c1;             // J/(mol·K)        (no conversion)
    hkf.c2     = params.c2;             // J·K/mol          (no conversion)
    return PerpleX::preprocessHKFParams(hkf);
}

auto evalGibbs(const PerpleX::HKFParams& hkf, real T, real Pbar) -> real
{
    double waterVolume = 0.0;
    const auto solvent = PerpleX::getWaterSolventState(Pbar, T, waterVolume);
    const auto hkfState = PerpleX::computeHKFGibbs(hkf, Pbar, T, solvent.epsilon, solvent.gf);
    return hkfState.G;
}

} // namespace

auto StandardThermoModelPerplexDEW(const StandardThermoModelParamsPerplexDEW& params) -> StandardThermoModel
{
    const auto hkf = makePerpleXHKFParams(params);

    auto evalfn = [=](StandardThermoProps& props, real T, real P)
    {
        auto& [G0, H0, V0, Cp0, VT0, VP0] = props;

        const auto Tmax = params.Tmax > 0.0 ? params.Tmax : real(1200.0);
        const auto Tclip = std::min(T, Tmax);
        const auto Pbar = P / 1.0e5;

        const auto dT = 0.1;
        const auto dPbar = 0.05;
        const auto Tm = std::max(real(273.16), Tclip - dT);
        const auto Tp = Tclip + dT;
        const auto Pmbar = std::max(real(1.0e-3), Pbar - dPbar);
        const auto Ppbar = Pbar + dPbar;

        const auto G = evalGibbs(hkf, Tclip, Pbar);
        const auto GTm = evalGibbs(hkf, Tm, Pbar);
        const auto GTp = evalGibbs(hkf, Tp, Pbar);
        const auto GPm = evalGibbs(hkf, Tclip, Pmbar);
        const auto GPp = evalGibbs(hkf, Tclip, Ppbar);
        const auto GTmPm = evalGibbs(hkf, Tm, Pmbar);
        const auto GTmPp = evalGibbs(hkf, Tm, Ppbar);
        const auto GTpPm = evalGibbs(hkf, Tp, Pmbar);
        const auto GTpPp = evalGibbs(hkf, Tp, Ppbar);

        const auto dGdT = (GTp - GTm) / (Tp - Tm);
        const auto d2GdT2 = (GTp - 2.0 * G + GTm) / ((Tp - Tclip) * (Tp - Tclip));

        const auto dPbarEff = Ppbar - Pmbar;
        const auto dGdPbar = (GPp - GPm) / dPbarEff;
        const auto d2GdPbar2 = (GPp - 2.0 * G + GPm) / ((dPbarEff * 0.5) * (dPbarEff * 0.5));
        const auto d2GdTdPbar = (GTpPp - GTpPm - GTmPp + GTmPm) / ((Tp - Tm) * dPbarEff);

        const auto S0 = -dGdT;
        G0 = G;
        H0 = G0 + Tclip * S0;
        Cp0 = -Tclip * d2GdT2;

        // Convert bar derivatives to Pa derivatives for SI consistency.
        V0 = dGdPbar / 1.0e5;
        VT0 = d2GdTdPbar / 1.0e5;
        VP0 = d2GdPbar2 / 1.0e10;
    };

    Data paramsdata;
    paramsdata["PerplexDEW"]  = true;
    paramsdata["omega0"]      = hkf.omega0;     // J/mol — Born coefficient ω₀ at ref conditions
    paramsdata["bornRadius"]  = hkf.bornRadius; // Å — effective Born radius (from preprocessHKFParams)
    paramsdata["charge"]      = hkf.charge;     // unitless ionic charge

    return StandardThermoModel(evalfn, paramsdata);
}

} // namespace Reaktoro
