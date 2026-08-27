#include "StandardThermoModelPerplexGFSM.hpp"

// C++ includes
#include <array>
#include <cmath>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Extensions/DEW/WaterEosZhangDuan2005.hpp>

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

    // -----------------------------------------------------------------------
    // Special case: H2O (speciesIndex == 1) in DEW/aqueous context.
    //
    // GFSMFluidModel::compute() routes H2O through the MRK gas-mixture
    // machinery before substituting the ZD05 correction.  At low pressures
    // (e.g. the 1-bar standard-state reference) or in an aqueous equilibration
    // context, MRK can be numerically unstable, causing an access violation.
    //
    // ActivityModelPerplexDEW already computes water ACTIVITY internally via
    // the same ZD05 EoS (through getWaterSolventState / computeFluidSolventState).
    // What is needed here is only the standard Gibbs energy G°(T,P) for water.
    //
    // We compute this directly from waterThermoPropsZhangDuan2005, anchored to
    // params.G0 at (Tref, Pref), using analytical derivatives for V, Cp, etc.
    // This is thermodynamically consistent with ActivityModelPerplexDEW and
    // avoids the MRK machinery entirely.
    // -----------------------------------------------------------------------
    if(params.speciesIndex == 1)
    {
        constexpr double Mwater = 0.01801528; // kg/mol (molar mass of H2O)

        // Capture standard-state properties at (Tref, Pref = 1 bar).
        const auto wtp_ref = waterThermoPropsZhangDuan2005(Tref, Pref);
        const double G_ref = static_cast<double>(wtp_ref.G); // J/kg at (Tref,Pref)
        const double H_ref = static_cast<double>(wtp_ref.H); // J/kg at (Tref,Pref)

        // Anchor to user-provided reference G0 and H0 from the DEW database.
        const double G0_anchor = static_cast<double>(params.G0); // J/mol
        const double H0_anchor = static_cast<double>(params.H0); // J/mol

        return [=](real T, real P) -> StandardThermoProps
        {
            // Evaluate ZD05 water EoS at (T, P).
            const auto wtp = waterThermoPropsZhangDuan2005(T, P);

            const double G_val  = static_cast<double>(wtp.G);  // J/kg
            const double H_val  = static_cast<double>(wtp.H);  // J/kg
            const double D_val  = static_cast<double>(wtp.D);  // kg/m³
            const double DT_val = static_cast<double>(wtp.DT); // (kg/m³)/K
            const double DP_val = static_cast<double>(wtp.DP); // (kg/m³)/Pa
            const double Cp_val = static_cast<double>(wtp.Cp); // J/(kg·K)

            // Standard molar Gibbs energy anchored to database reference:
            //   G°(T,P) = G0_anchor + Mwater × [G_ZD05(T,P) − G_ZD05(Tref,Pref)]
            const double G0 = G0_anchor + Mwater * (G_val - G_ref);

            // Standard molar enthalpy anchored to database reference:
            //   H°(T,P) = H0_anchor + Mwater × [H_ZD05(T,P) − H_ZD05(Tref,Pref)]
            const double H0 = H0_anchor + Mwater * (H_val - H_ref);

            // Molar volume: V = Mwater / D  [m³/mol]
            const double V0 = (D_val > 0.0) ? Mwater / D_val : 0.0;

            // Molar heat capacity: Cp = Mwater × Cp_specific  [J/(mol·K)]
            const double Cp0 = Mwater * Cp_val;

            // ∂V/∂T at const P = −Mwater × DT / D²  [m³/(mol·K)]
            const double VT0 = (D_val > 0.0) ? -Mwater * DT_val / (D_val * D_val) : 0.0;

            // ∂V/∂P at const T = −Mwater × DP / D²  [m³/(mol·Pa)]
            const double VP0 = (D_val > 0.0) ? -Mwater * DP_val / (D_val * D_val) : 0.0;

            StandardThermoProps props;
            props.G0  = G0;
            props.H0  = H0;
            props.V0  = V0;
            props.Cp0 = Cp0;
            props.VT0 = VT0;
            props.VP0 = VP0;
            return props;
        };
    }

    // Non-hybrid gases (all species except H2O/CO2/CH4) can trigger unstable
    // behavior in pure-species GFSM/MRK evaluation at the standard-state
    // reference conditions used by fugacity constraints. For these species,
    // use the database-provided constant reference model.
    if(params.speciesIndex != 2 && params.speciesIndex != 4)
    {
        return [=](real /*T*/, real /*P*/) -> StandardThermoProps
        {
            StandardThermoProps props;
            props.G0 = params.G0;
            props.H0 = params.H0;
            props.V0 = params.V0;
            props.Cp0 = 0.0;
            props.VT0 = 0.0;
            props.VP0 = 0.0;
            return props;
        };
    }

    // -----------------------------------------------------------------------
    // General case: non-water GFSM species (CO2, CO, CH4, H2, H2S, …)
    // Use GFSMFluidModel::compute() via the hybrid EoS framework.
    // -----------------------------------------------------------------------
    auto buildGfsmOptions = [&]() -> PerpleX::GFSMFluidOptions
    {
        PerpleX::GFSMFluidOptions opts;
        opts.hybridEosOptions = params.hybridEosOptions;
        opts.mrkMixOptions = params.mrkMixOptions;
        opts.useLowTMrk = params.useLowTMrk;
        opts.enableElectrolyte = false;

        // Pure-species standard-state evaluation calls GFSM with exactly one
        // species in the mixture. Only include that species in hybrid handling
        // when it is one of the hybrid-capable species.
        if(params.speciesIndex == 1 || params.speciesIndex == 2 || params.speciesIndex == 4)
            opts.hybridSpeciesIndices = {params.speciesIndex};

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
