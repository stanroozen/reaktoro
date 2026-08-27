#include "ActivityModelPerplexGFSM.hpp"

// C++ includes
#include <cmath>
#include <map>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

auto ActivityModelPerplexGFSM(
    const ActivityModelParamsPerplexGFSM& params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& species) -> ActivityModel
    {
        // Map Reaktoro species names to Perple_X indices (1-13)
        const std::map<std::string, int> perplex_indices = {
            {"H2O", 1},
            {"CO2", 2},
            {"CO", 3},
            {"CH4", 4},
            {"H2", 5},
            {"H2S", 6},
            {"O2", 7},
            {"SO2", 8},
            {"N2", 10},
            {"NH3", 11},
            {"HF", 17},
            {"C2H6", 16},
            {"HCl", 18}
        };

        // Build species index mapping and validate
        std::vector<int> perplex_species;
        std::map<int, size_t> perplex_to_reaktoro_idx;

        for(size_t i = 0; i < species.size(); ++i)
        {
            const auto& sp = species[i];
            std::string base_name = sp.name();

            // Try to match base formula
            auto it = perplex_indices.find(sp.formula().str());
            if(it == perplex_indices.end())
            {
                it = perplex_indices.find(base_name);
            }
            if(it == perplex_indices.end())
            {
                // Try without (g) suffix
                if(base_name.size() > 2 && base_name.substr(base_name.size()-3) == "(g)")
                {
                    base_name = base_name.substr(0, base_name.size()-3);
                    it = perplex_indices.find(base_name);
                }
            }

            errorifnot(it != perplex_indices.end(),
                "Species '", sp.name(), "' is not a valid Perple_X GFSM fluid species. ",
                "Valid species: H2O, CO2, CH4, H2, CO, H2S, O2, SO2, N2, NH3, HF, C2H6, HCl");

            const int px_idx = it->second;
            perplex_species.push_back(px_idx);
            perplex_to_reaktoro_idx[px_idx] = i;
        }

        errorifnot(perplex_to_reaktoro_idx.count(1),
            "ActivityModelPerplexGFSM requires H2O to be present in the gaseous phase "
            "(Perple_X species index 1). The provided phase has no H2O species, which "
            "is unsupported for GFSM and can lead to invalid EOS evaluations.");

        // Create Perple_X GFSM model instance
        PerpleX::GFSMFluidModel model;

        // Configure GFSM options
        PerpleX::GFSMFluidOptions gfsm_options;
        gfsm_options.mrkMixOptions = params.mrkMixOptions;
        gfsm_options.useLowTMrk = params.useLowTMrk;
        gfsm_options.enableElectrolyte = params.enableElectrolyte;
        gfsm_options.hybridEosOptions = params.hybridEosOptions;

        // Only keep hybrid species that are actually present in this phase.
        // Passing absent species indices downstream can trigger invalid paths in
        // low-level hybrid routines during fugacity-constrained solves.
        const std::array<int, 3> default_hybrid_species = {1, 2, 4}; // H2O, CO2, CH4
        for(const int idx : default_hybrid_species)
        {
            if(perplex_to_reaktoro_idx.count(idx))
                gfsm_options.hybridSpeciesIndices.push_back(idx);
        }

        // Define the activity model function
        ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args)
        {
            // Extract T, P, composition
            const auto& [T, P, x] = args;

            // Convert to SI units for Perple_X
            const double T_K = T;
            const double P_bar = P / 1.0e5;  // Pa to bar

            // Prepare species composition array for Perple_X (19-element)
            std::array<double, 19> y{};
            for(size_t i = 0; i < species.size(); ++i)
            {
                const int px_idx = perplex_species[i];
                y[px_idx] = x[i];
            }

            // Compute GFSM fluid state
            const auto gfsm_state = model.compute(
                perplex_species,
                y,
                P_bar,
                T_K,
                gfsm_options);

            // Set state of matter
            props.som = StateOfMatter::Fluid;

            // Transfer molar volume and derivatives
            props.Vx = gfsm_state.molarVolume * 1.0e-6;  // cm³/mol to m³/mol

            // Initialize Gibbs energy, enthalpy, heat capacity (not computed by GFSM)
            props.Gx = 0.0;
            props.Hx = 0.0;
            props.Cpx = 0.0;
            props.VxT = 0.0;
            props.VxP = 0.0;

            // Transfer activity coefficient and volume arrays
            props.ln_g = ArrayXr::Zero(species.size());
            props.ln_a = ArrayXr::Zero(species.size());
            props.Vxi = ArrayXr::Zero(species.size());

            for(size_t i = 0; i < species.size(); ++i)
            {
                const int px_idx = perplex_species[i];
                props.ln_g[i] = std::log(gfsm_state.g_mrk[px_idx]);
                props.ln_a[i] = gfsm_state.ln_f[px_idx];
                props.Vxi[i] = gfsm_state.v_mrk[px_idx] * 1.0e-6;  // cm³/mol to m³/mol
            }

            // Publish a coupled-fluid handoff payload for aqueous activity models.
            // This enables gas->aqueous water-activity coupling even when CO2 exists
            // only in the gas phase and is absent from aqueous neutral species.
            auto water_it = perplex_to_reaktoro_idx.find(1); // Perple_X index 1 = H2O
            if(water_it != perplex_to_reaktoro_idx.end())
            {
                std::array<double, 19> y_pure{};
                y_pure[1] = 1.0;

                auto pure_opts = gfsm_options;
                pure_opts.enableElectrolyte = true;

                const auto pure_state = model.compute(
                    perplex_species,
                    y_pure,
                    P_bar,
                    T_K,
                    pure_opts);

                const double ln_f_mix_h2o = gfsm_state.ln_f[1];
                const double ln_f_pure_h2o = pure_state.ln_f[1];
                const double ln_f_ratio_h2o = ln_f_mix_h2o - ln_f_pure_h2o;

                const auto stid_it = props.extra.find("Reaktoro::ChemicalProps::StateId");
                if(stid_it != props.extra.end())
                    props.extra["PerplexGFSM::WaterActivity::StateId"] = stid_it->second;

                props.extra["PerplexGFSM::WaterActivity::Source"] = String("coupled-fluid");
                props.extra["PerplexGFSM::WaterActivity::ln_f_h2o_mix"] = ln_f_mix_h2o;
                props.extra["PerplexGFSM::WaterActivity::ln_f_h2o_pure"] = ln_f_pure_h2o;
                props.extra["PerplexGFSM::WaterActivity::ln_f_ratio_h2o"] = ln_f_ratio_h2o;
                props.extra["PerplexGFSM::WaterActivity::x_h2o_fluid"] = y[1];

                // Optional solvent state payload for PerplexDEW mixed-fluid dielectric/Born updates.
                props.extra["PerplexGFSM::WaterActivity::epsilon_mix"] = gfsm_state.dielectric.epsilon;
                props.extra["PerplexGFSM::WaterActivity::epsilon_pure"] = pure_state.dielectric.epsilon;
                props.extra["PerplexGFSM::WaterActivity::gf_mix"] = gfsm_state.dielectric.gf;
                props.extra["PerplexGFSM::WaterActivity::gf_pure"] = pure_state.dielectric.gf;
                props.extra["PerplexGFSM::WaterActivity::msol_mix"] = gfsm_state.dielectric.msol;
                props.extra["PerplexGFSM::WaterActivity::vsolv_mix"] = gfsm_state.dielectric.vsolv;
            }
        };

        return fn;
    };
}

} // namespace Reaktoro
