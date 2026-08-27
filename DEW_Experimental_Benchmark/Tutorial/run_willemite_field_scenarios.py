import importlib.util
import os
import re

import numpy as np


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TUTORIAL_PATH = os.path.join(
    SCRIPT_DIR,
    "willemite_solubility_tutorial_dew17hp622_zn.py",
)


def slugify(name):
    slug = name.lower().strip()
    slug = re.sub(r"[^a-z0-9]+", "_", slug)
    slug = re.sub(r"_+", "_", slug).strip("_")
    return slug


def load_tutorial_module(path):
    spec = importlib.util.spec_from_file_location("willemite_tutorial", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def scenario_output_path(base_path, slug):
    root, ext = os.path.splitext(base_path)
    return f"{root}_{slug}{ext}"


def run_scenario(module, name, h2o_mol, wlm_mol, ph_min, ph_max, nacl_m):
    slug = slugify(name)

    module.INITIAL_SPECIES_AMOUNTS_MOL = dict(module._orig_initial)
    module.INITIAL_SPECIES_AMOUNTS_MOL["H2O"] = float(h2o_mol)
    module.INITIAL_SPECIES_AMOUNTS_MOL[module.MINERAL_NAME] = float(wlm_mol)
    module.PH_RANGE = np.linspace(float(ph_min), float(ph_max), int(module.SENS_POINTS))
    module.USE_NACL_BRINE_BACKGROUND = True
    module.NACL_BRINE_MOL_PER_KG_H2O = float(nacl_m)

    module.OUTPUT_T_CURVE = scenario_output_path(module._orig_output_t_curve, slug)
    module.OUTPUT_SENS = scenario_output_path(module._orig_output_sens, slug)
    module.OUTPUT_TRUE_FO2 = scenario_output_path(module._orig_output_true_fo2, slug)

    solvent_species_name = module.infer_solvent_species_name(
        module.AQUEOUS_SPECIES,
        module.SOLVENT_SPECIES_NAME,
        module.INITIAL_SPECIES_AMOUNTS_MOL,
    )

    system_no_gas = module.build_tutorial_system()
    temperatures_c, curves = module.compute_temperature_curves(
        system_no_gas, solvent_species_name
    )
    module.save_temperature_plot(temperatures_c, curves)

    pH_x, pH_y = module.compute_pH_sensitivity(system_no_gas, solvent_species_name)

    mu_sio2_x, mu_sio2_y = module.compute_silica_potential_sensitivity(
        system_no_gas,
        solvent_species_name,
        module.MU_SIO2_RANGE,
    )

    system_with_h2s_gas = module.build_system_with_h2s_gas_phase()
    fH2S_x, fH2S_y = module.compute_true_fh2s_sensitivity(
        system_with_h2s_gas,
        solvent_species_name,
    )

    system_with_o2_gas = module.build_system_with_o2_gas_phase()
    fO2_x, fO2_y = module.compute_true_fo2_sensitivity(
        system_with_o2_gas,
        solvent_species_name,
    )

    module.save_sensitivity_plot(
        pH_x,
        pH_y,
        mu_sio2_x,
        mu_sio2_y,
        fH2S_x,
        fH2S_y,
        fO2_x,
        fO2_y,
    )

    module.save_true_fo2_plot(fO2_x, fO2_y)

    return {
        "temperature": module.OUTPUT_T_CURVE,
        "sensitivity": module.OUTPUT_SENS,
        "fo2": module.OUTPUT_TRUE_FO2,
    }


def main():
    module = load_tutorial_module(TUTORIAL_PATH)

    try:
        module.Warnings.disable(906)
    except Exception:
        pass

    module._orig_initial = dict(module.INITIAL_SPECIES_AMOUNTS_MOL)
    module._orig_ph = np.array(module.PH_RANGE, dtype=float)
    module._orig_brine = float(module.NACL_BRINE_MOL_PER_KG_H2O)
    module._orig_use_brine = bool(module.USE_NACL_BRINE_BACKGROUND)
    module._orig_output_t_curve = module.OUTPUT_T_CURVE
    module._orig_output_sens = module.OUTPUT_SENS
    module._orig_output_true_fo2 = module.OUTPUT_TRUE_FO2

    scenarios = [
        ("field_lite_a", 555.0, 0.1, 4.5, 8.5, 0.2),
        ("field_lite_b", 555.0, 0.05, 4.5, 8.5, 0.2),
        ("field_mid_c", 1110.0, 0.05, 5.0, 8.0, 0.1),
        ("field_dilute_d", 5550.0, 0.01, 5.0, 8.0, 0.1),
    ]

    created = {}

    try:
        for scenario in scenarios:
            name = scenario[0]
            print("=" * 78)
            print(f"Running scenario: {name}")
            created[name] = run_scenario(module, *scenario)
    finally:
        module.INITIAL_SPECIES_AMOUNTS_MOL = dict(module._orig_initial)
        module.PH_RANGE = np.array(module._orig_ph, dtype=float)
        module.NACL_BRINE_MOL_PER_KG_H2O = module._orig_brine
        module.USE_NACL_BRINE_BACKGROUND = module._orig_use_brine
        module.OUTPUT_T_CURVE = module._orig_output_t_curve
        module.OUTPUT_SENS = module._orig_output_sens
        module.OUTPUT_TRUE_FO2 = module._orig_output_true_fo2

    print("\nGenerated files:")
    for name, outputs in created.items():
        print(f"- {name}")
        for key, path in outputs.items():
            print(f"  {key}: {path}")


if __name__ == "__main__":
    main()
