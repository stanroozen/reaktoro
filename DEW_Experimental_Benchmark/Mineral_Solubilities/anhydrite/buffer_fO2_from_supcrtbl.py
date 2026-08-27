"""
Compute buffer fugacities from SUPCRTBL thermodynamics.

This script does NOT add buffer phases to a chemical system. It only computes
fugacities from thermodynamic data in SUPCRTBL, following the fixed-fugacity
specification approach (EquilibriumSpecs/EquilibriumConditions).
"""

import os
import sys
import math

# Try to import Reaktoro; fall back to local extension module if not installed
try:
    from reaktoro import *  # noqa: F401,F403
except ModuleNotFoundError:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
    ROOT_DIR = os.path.dirname(BENCHMARK_DIR)
    PYD_DIR = os.path.join(ROOT_DIR, "build", "Reaktoro", "Release")
    if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
        sys.path.insert(0, PYD_DIR)
    try:
        from reaktoro4py import *  # noqa: F401,F403

        print("Using local reaktoro4py extension from build.")
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "Could not import Reaktoro. Install the 'reaktoro' package or ensure reaktoro4py is on PYTHONPATH."
        ) from e

# Path to SUPCRTBL database file embedded in this repository.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_SUPCRTBL_DB_FILE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(_THIS_DIR))),
    "embedded",
    "databases",
    "reaktoro",
    "supcrtbl.yaml",
)


def _load_supcrtbl_db():
    """Load SUPCRTBL database from the embedded file."""
    return Database.fromFile(_SUPCRTBL_DB_FILE)


BUFFER_REACTIONS = {
    "fO2": {
        # IW: Fe + 1/2 O2 = FeO (SUPCRTBL uses Ferropericlase for FeO)
        "IW": [
            ("Ferropericlase", 1.0),
            ("Iron", -1.0),
            ("O2(g)", 0.5),
        ],
        # FMQ: 3 Fe2SiO4 + O2 = 2 Fe3O4 + 3 SiO2
        "FMQ": [
            ("Magnetite", 2.0),
            ("Quartz", 3.0),
            ("Fayalite", -3.0),
            ("O2(g)", 1.0),
        ],
        # NNO: Ni + 1/2 O2 = NiO
        "NNO": [
            ("Nickel-Oxide", 1.0),
            ("Nickel", -1.0),
            ("O2(g)", 0.5),
        ],
        # MH: 2 Fe3O4 + 1/2 O2 = 3 Fe2O3
        "MH": [
            ("Hematite", 3.0),
            ("Magnetite", -2.0),
            ("O2(g)", 0.5),
        ],
        # HM alias
        "HM": [
            ("Hematite", 3.0),
            ("Magnetite", -2.0),
            ("O2(g)", 0.5),
        ],
        # Mn2O3â€“MnO: 2 MnO + 1/2 O2 = Mn2O3
        "Mn2O3-MnO": [
            ("Bixbyite", 1.0),
            ("Manganosite", -2.0),
            ("O2(g)", 0.5),
        ],
        # CCO: C + O2 = CO2
        "CCO": [
            ("CO2(g)", 1.0),
            ("Graphite", -1.0),
            ("O2(g)", 1.0),
        ],
    },
    "fCO2": {
        # Câ€“Qâ€“Wo: CaCO3 + SiO2 = CaSiO3 + CO2
        "C-Q-Wo": [
            ("Wollastonite", 1.0),
            ("CO2(g)", 1.0),
            ("Calcite", -1.0),
            ("Quartz", -1.0),
        ],
        # Dolâ€“Calâ€“Qtz: CaMg(CO3)2 + SiO2 = CaCO3 + MgSiO3 + CO2
        "Dol-Cal-Qtz": [
            ("Calcite", 1.0),
            ("Enstatite", 1.0),
            ("CO2(g)", 1.0),
            ("Dolomite", -1.0),
            ("Quartz", -1.0),
        ],
        # Magâ€“Enâ€“Qtz: MgCO3 + SiO2 = MgSiO3 + CO2
        "Mag-En-Qtz": [
            ("Enstatite", 1.0),
            ("CO2(g)", 1.0),
            ("Magnesite", -1.0),
            ("Quartz", -1.0),
        ],
        # Grâ€“COâ€“CO2: C + O2 = CO2 (uses gas O2; included for completeness)
        "Gr-CO-CO2": [
            ("CO2(g)", 1.0),
            ("Graphite", -1.0),
            ("O2(g)", 1.0),
        ],
    },
    "fS2": {
        # Pyâ€“Po: FeS2 = FeS + 1/2 S2 (use pyrrhotite,trot for FeS)
        "Py-Po": [
            ("Pyrrhotite,trot", 1.0),
            ("S2(g)", 0.5),
            ("Pyrite", -1.0),
        ],
        # Poâ€“Mt: FeS + O2 = Fe3O4 + S2 (stoichiometry generalized)
        "Po-Mt": [
            ("Magnetite", 1.0),
            ("S2(g)", 1.0),
            ("Pyrrhotite,trot", -1.0),
            ("O2(g)", -1.0),
        ],
        # Pyâ€“Mt: FeS2 + O2 = Fe3O4 + S2 (stoichiometry generalized)
        "Py-Mt": [
            ("Magnetite", 1.0),
            ("S2(g)", 1.0),
            ("Pyrite", -1.0),
            ("O2(g)", -1.0),
        ],
    },
    "fH2": {
        # IWâ€“H2: Fe + H2O = FeO + H2 (FeO -> Ferropericlase)
        "IW-H2": [
            ("Ferropericlase", 1.0),
            ("H2(g)", 1.0),
            ("Iron", -1.0),
            ("H2O(g)", -1.0),
        ],
        # FMQâ€“H2: Fe2SiO4 + H2O = Fe3O4 + H2 + SiO2 (generalized)
        "FMQ-H2": [
            ("Magnetite", 1.0),
            ("Quartz", 1.0),
            ("H2(g)", 1.0),
            ("Fayalite", -1.0),
            ("H2O(g)", -1.0),
        ],
        # Câ€“H2O: C + H2O = CO + H2
        "C-H2O": [
            ("CO(g)", 1.0),
            ("H2(g)", 1.0),
            ("Graphite", -1.0),
            ("H2O(g)", -1.0),
        ],
    },
}

BUFFER_ALIASES = {
    "hm": ("fO2", "HM"),
    "m-h": ("fO2", "MH"),
    "hematite-magnetite": ("fO2", "HM"),
    "mn2o3-mno": ("fO2", "Mn2O3-MnO"),
    "ni-nio": ("fO2", "NNO"),
    "iw": ("fO2", "IW"),
    "fmq": ("fO2", "FMQ"),
    "nno": ("fO2", "NNO"),
    "mh": ("fO2", "MH"),
    "cco": ("fO2", "CCO"),
    "c-q-wo": ("fCO2", "C-Q-Wo"),
    "dol-cal-qtz": ("fCO2", "Dol-Cal-Qtz"),
    "mag-en-qtz": ("fCO2", "Mag-En-Qtz"),
    "gr-co-co2": ("fCO2", "Gr-CO-CO2"),
    "py-po": ("fS2", "Py-Po"),
    "po-mt": ("fS2", "Po-Mt"),
    "py-mt": ("fS2", "Py-Mt"),
    "iw-h2": ("fH2", "IW-H2"),
    "fmq-h2": ("fH2", "FMQ-H2"),
    "c-h2o": ("fH2", "C-H2O"),
}


def _canonical_buffer_name(name, category=None):
    if name is None:
        return None
    key = str(name).strip()
    if category and key in BUFFER_REACTIONS.get(category, {}):
        return category, key
    lowered = key.lower()
    if lowered in BUFFER_ALIASES:
        return BUFFER_ALIASES[lowered]
    if category and lowered in BUFFER_REACTIONS.get(category, {}):
        return category, lowered
    return None


def _species_g0(db, name, T_K, P_bar):
    P_pa = float(P_bar) * 1e5
    return db.species().get(name).props(float(T_K), P_pa).G0


def _validate_species(db, reaction):
    missing = []
    for name, _ in reaction:
        try:
            db.species().get(name)
        except Exception:
            missing.append(name)
    return missing


def buffer_log10_fugacity_bar(category, buffer_name, T_C, P_bar, db=None):
    """Return log10(f_gas/bar) for a buffer at T (Â°C) and P (bar)."""
    canon = _canonical_buffer_name(buffer_name, category=category)
    if canon is None:
        raise KeyError(f"Unknown buffer '{buffer_name}' in category '{category}'.")
    category_key, buffer_key = canon

    if db is None:
        db = _load_supcrtbl_db()

    if category_key not in BUFFER_REACTIONS:
        raise KeyError(f"Unknown buffer category '{category_key}'.")

    reaction = BUFFER_REACTIONS[category_key][buffer_key]
    missing = _validate_species(db, reaction)
    if missing:
        raise KeyError(
            f"Missing species in SUPCRTBL for buffer '{buffer_key}': {', '.join(missing)}"
        )

    gas_map = {
        "fO2": "O2(g)",
        "fCO2": "CO2(g)",
        "fS2": "S2(g)",
        "fH2": "H2(g)",
    }
    gas_name = gas_map.get(category_key)
    if not gas_name:
        raise KeyError(f"No gas species defined for category '{category_key}'.")

    T_K = float(T_C) + 273.15
    sum_g = 0.0
    nu_gas = None
    g0_gas = None

    for species_name, nu in reaction:
        if species_name == gas_name:
            nu_gas = float(nu)
            g0_gas = _species_g0(db, species_name, T_K, P_bar)
        else:
            sum_g += float(nu) * _species_g0(db, species_name, T_K, P_bar)

    if nu_gas is None or nu_gas == 0.0:
        raise ValueError(f"Reaction for '{buffer_key}' missing {gas_name} term.")

    try:
        R = float(universalGasConstant)
    except Exception:
        R = 8.314462618

    ln_f = (sum_g - nu_gas * g0_gas) / (nu_gas * R * T_K)
    return ln_f / math.log(10.0)


def buffer_fugacity_bar(category, buffer_name, T_C, P_bar, db=None):
    """Return f_gas in bar for a buffer at T (Â°C) and P (bar)."""
    return 10.0 ** buffer_log10_fugacity_bar(category, buffer_name, T_C, P_bar, db=db)


def main():
    T_C = 800.0
    P_kbar = 10.0
    P_bar = P_kbar * 1000.0

    db = _load_supcrtbl_db()
    print(f"T = {T_C:.1f} °C, P = {P_kbar:.2f} kbar")

    for category, buffers in BUFFER_REACTIONS.items():
        print(f"\n{category} buffers:")
        for buf in buffers.keys():
            try:
                logf = buffer_log10_fugacity_bar(category, buf, T_C, P_bar, db=db)
                print(f"  {buf:12s}: log10 f(bar) = {logf: .4f}")
            except Exception as e:
                print(f"  {buf:12s}: unavailable ({e})")


if __name__ == "__main__":
    main()
