#include "StandardThermoModelPerplexGFSM.hpp"

// C++ includes
#include <array>
#include <cmath>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

auto StandardThermoModelPerplexGFSM(
    const StandardThermoModelParamsPerplexGFSM& params) -> StandardThermoModel
{
    constexpr double R = 8.31446261815324; // J/(mol*K)
    constexpr double Tref = 298.15;        // K
    constexpr double Pref = 1.0e5;         // Pa

    errorif(params.speciesIndex < 1 || params.speciesIndex > 18,
        "Invalid Perple_X species index ", params.speciesIndex,
        " in StandardThermoModelPerplexGFSM. Expected range is [1, 18].");

    auto buildGfsmOptions = [&]() -> PerpleX::GFSMFluidOptions
    {
        PerpleX::GFSMFluidOptions opts;
        opts.hybridEosOptions = params.hybridEosOptions;
        opts.mrkMixOptions = params.mrkMixOptions;
        opts.useLowTMrk = params.useLowTMrk;
        opts.enableElectrolyte = false;
        opts.hybridSpeciesIndices = {1, 2, 4}; // Perple_X hybrid-capable species
        return opts;
    };

    auto evalLnF = [&](double T, double P) -> double
    {
        PerpleX::GFSMFluidModel model;
        const auto opts = buildGfsmOptions();
        std::vector<int> species = {params.speciesIndex};
        std::array<double, 19> y{};
        y[params.speciesIndex] = 1.0;

        const auto state = model.compute(species, y, P / 1.0e5, T, opts); // Pa->bar
        const double lnf = state.ln_f[params.speciesIndex];
        errorif(!std::isfinite(lnf),
            "GFSM ln(f) became non-finite for species ", params.speciesIndex,
            " at T=", T, " K, P=", P, " Pa.");
        return lnf;
    };

    const double lnf_ref = evalLnF(Tref, Pref);

    return [=](real T, real P) -> StandardThermoProps
    {
        // Temperature range check
        errorif(T > params.Tmax,
            "Temperature ", T, " K exceeds maximum ", params.Tmax, " K for Perple_X GFSM model");

        const double Tval = static_cast<double>(T);
        const double Pval = static_cast<double>(P);

        // Finite-difference spacing (stable across wide P-T range).
        const double dT = std::max(1.0e-2, 1.0e-3 * Tval);
        const double dP = std::max(1.0e2, 1.0e-5 * std::max(Pval, 1.0e5));

        const double Tm = std::max(273.16, Tval - dT);
        const double Tp = Tval + dT;
        const double Pm = std::max(1.0, Pval - dP);
        const double Pp = Pval + dP;

        // Anchor to user-provided reference Gibbs at (Tref, Pref), and use GFSM ln(f)
        // as the P-T correction term for the pure reference species.
        auto evalG = [&](double Te, double Pe) -> double
        {
            const double lnf = evalLnF(Te, Pe);
            return static_cast<double>(params.G0) + R * Te * (lnf - lnf_ref);
        };

        const double G = evalG(Tval, Pval);
        const double GTm = evalG(Tm, Pval);
        const double GTp = evalG(Tp, Pval);
        const double GPm = evalG(Tval, Pm);
        const double GPp = evalG(Tval, Pp);
        const double GTmPm = evalG(Tm, Pm);
        const double GTmPp = evalG(Tm, Pp);
        const double GTpPm = evalG(Tp, Pm);
        const double GTpPp = evalG(Tp, Pp);

        const double dTd = Tp - Tm;
        const double dPd = Pp - Pm;
        const double dGdT = (GTp - GTm) / dTd;
        const double d2GdT2 = (GTp - 2.0 * G + GTm) / std::pow((Tp - Tm) * 0.5, 2);
        const double dGdP = (GPp - GPm) / dPd;
        const double d2GdP2 = (GPp - 2.0 * G + GPm) / std::pow((Pp - Pm) * 0.5, 2);
        const double d2GdTdP = (GTpPp - GTpPm - GTmPp + GTmPm) / (dTd * dPd);

        const double S = -dGdT;

        StandardThermoProps props;
        props.G0 = G;
        props.H0 = G + Tval * S;
        props.V0 = dGdP;
        props.Cp0 = -Tval * d2GdT2;
        props.VT0 = d2GdTdP;
        props.VP0 = d2GdP2;

        return props;
    };
}

} // namespace Reaktoro
