// WaterState.cpp

#include <Reaktoro/Extensions/DEW/WaterState.hpp>

namespace Reaktoro {

auto waterState(real T,
                real P,
                const WaterStateOptions& opts) -> WaterState
{
    WaterState ws;

    // 1) Core EOS: WaterThermoProps (always)
    ws.thermo = waterThermoPropsModel(T, P, opts.thermo);

    // 2) Core dielectric: WaterElectroProps (always)
    //    Uses ws.thermo.D and ws.thermo.DP internally (and optional Psat overrides)
    ws.electro = waterElectroPropsModel(T, P, ws.thermo, opts.dielectric);

    // 3) Optional: Gibbs free energy of water
    if (opts.computeGibbs)
    {
        ws.gibbs    = waterGibbsModel(T, P, opts.gibbs);
        ws.hasGibbs = true;
    }

    // 4) Optional: DEW solvent function g(P,T,ρ) and its pressure derivative
    if (opts.computeSolventG)
    {
        ws.g_solv = waterSolventFunctionDEW(T, P, ws.thermo, opts.solvent);
        ws.dgdP   = waterSolventFunctionDgdP_DEW(T, P, ws.thermo,
                                                ws.g_solv, opts.solvent);
        ws.hasSolventG = true;
    }

    // 5) Optional: DEW Born omega(P,T) and dω/dP
    //
    // The Born formulas depend on:
    //   - thermo (ρ, dρ/dP)
    //   - solvent function g, dgdP (via opts.omega.solvent)
    //
    // Callers supply wref & Z when they use these helpers for species.
    if (opts.computeOmega)
    {
        // Note:
        //   waterBornOmegaDEW / waterBornDOmegaDP_DEW require (wref, Z).
        //   WaterState itself doesn’t know species, so here we compute
        //   omega/domega_dP only if the user has encoded an effective
        //   "generic" system via opts.omega (e.g. by setting isHydrogenLike
        //   or using a representative wref/Z externally).
        //
        // Typical pattern for species:
        //   auto ws = waterState(T,P, base_opts);
        //   auto omega_i = waterBornOmegaDEW(T,P, ws.thermo, wref_i, Zi, optsOmega);
        //   ...
        //
        // To keep WaterState generic, we only evaluate omega using
        // the options struct provided (which may encode a specific Z/wref
        // in your higher-level API). If you prefer purely species-level
        // calls, you can ignore ws.omega/ws.domega_dP entirely.

        // If you don't bind species here, you can leave these zero and
        // just use the standalone waterBornOmegaDEW(...) APIs.
        ws.omega     = 0.0;
        ws.domega_dP = 0.0;
        ws.hasOmega  = false;
    }

    return ws;
}

} // namespace Reaktoro
