import math
import re
from pathlib import Path

import pandas as pd
import yaml

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CAL2J = 4.184  # cal → J
BAR2PA = 1.0e5  # bar → Pa
CM3_TO_M3 = 1.0e-6  # cm³ → m³ (not strictly needed here)
R = 8.31446261815324  # J/(mol·K)
TREF = 298.15  # K

# ---------------------------------------------------------------------------
# Special cases for missing / wrong formulas in DEW sheets
# ---------------------------------------------------------------------------

SPECIAL_AQUEOUS_SPECIES = {
    # Name in "Chemical" column : (Formula, extra_comment)
    "bG,NaCl": ("NaCl", "Background electrolyte (bG,NaCl in DEW)."),
    "Glycinate,aq": ("C2H4NO2-", None),
    "H2CO3,aq": ("H2CO3", None),
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def parse_elements_from_formula(formula: str) -> str:
    """
    Parse a chemical formula like 'CH3COO-' or 'Fe(OH)3' or 'Ar'
    into a Reaktoro 'Elements' string, e.g. '1:C 3:H 2:O'.

    This is a simple parser:
      - strips phase flags after comma (',g', ',aq', etc.)
      - strips charge at the end
      - reads tokens like 'C', 'H2', 'SiO4'
    Good enough for the DEW species formulas.
    """
    # Remove anything after comma (phase, etc.)
    f = re.split(r"[,\s]", formula)[0]

    # Remove parentheses around charges, e.g. 'Ag(+)' → 'Ag'
    f = re.sub(r"\([+-]?\d*\)", "", f)

    # Strip trailing charge, e.g. 'CH3COO-', 'Fe3+', 'SO4-2'
    f = re.sub(r"[\+\-]\d*$", "", f)
    f = re.sub(r"\([\+\-]\d+\)$", "", f)

    # Element tokens: one capital + optional lowercase + optional digits
    tokens = re.findall(r"([A-Z][a-z]?)(\d*)", f)
    elements = {}
    for el, count in tokens:
        n = int(count) if count else 1
        elements[el] = elements.get(el, 0) + n

    # Reaktoro expects e.g. "1:C 3:H 2:O"
    return " ".join(f"{n}:{el}" for el, n in elements.items())


def nasa_from_cp_poly(a_J, b_JK, c_JK2, Hf_J, S_J, Tmin, Tmax, label):
    """
    Build NASA polynomial parameters that reproduce:

        Cp(T) = a_J + b_JK*T + c_JK2*T^2   [J/(mol·K)]

    exactly (for all T in [Tmin, Tmax]) and match Hf° and S° at TREF.

    Returns a dict suitable for StandardThermoModelParamsNasa in Reaktoro.
    """
    # cp(T)/R = a1 + a2*T + a3*T²
    a1 = a_J / R
    a2 = b_JK / R
    a3 = c_JK2 / R
    a4 = 0.0
    a5 = 0.0
    a6 = 0.0
    a7 = 0.0

    T = TREF

    # NASA formulas (with a4=a5=a6=0):
    #   H(T)/(R T) = a1 + a2*T/2 + a3*T²/3 + b1/T
    #   S(T)/R     = a1 ln T + a2*T + a3*T²/2 + b2

    phi_H = a1 + a2 * T / 2.0 + a3 * T * T / 3.0
    phi_S = a1 * math.log(T) + a2 * T + a3 * T * T / 2.0

    b1 = Hf_J / R - phi_H * T  # enforce H(TREF) = Hf_J
    b2 = S_J / R - phi_S  # enforce S(TREF) = S_J

    poly = {
        "Tmin": float(Tmin),
        "Tmax": float(Tmax),
        "label": str(label),
        "state": "Gas",
        "a1": float(a1),
        "a2": float(a2),
        "a3": float(a3),
        "a4": float(a4),
        "a5": float(a5),
        "a6": float(a6),
        "a7": float(a7),
        "b1": float(b1),
        "b2": float(b2),
    }

    nasa_params = {
        "dHf": float(Hf_J),  # formation enthalpy at 298.15 K
        "dH0": 0.0,  # keep 0 for single interval
        "H0": float(Hf_J),  # enthalpy at T0
        "T0": float(TREF),
        "polynomials": [poly],
    }

    return nasa_params


# ---------------------------------------------------------------------------
# Load tables from a DEW Excel
# ---------------------------------------------------------------------------


def load_aqueous_table(xls_path: Path) -> pd.DataFrame:
    """
    Load 'Aqueous Species Table' from a DEW Excel (2019 or 2024).

    Structure:
      - row 2: header names (Chemical, ΔGfo, ΔHfo, So, Vo, ...)
      - row 3: units
      - data start at row 4
    """
    raw = pd.read_excel(xls_path, sheet_name="Aqueous Species Table", header=None)

    header_names = raw.iloc[2]  # row with labels
    df = raw.iloc[4:].copy()  # data below
    df.columns = header_names

    # Drop rows without a Chemical name
    df = df[df["Chemical"].notna()].copy()

    # Drop first dummy column and keep real data
    df = df.iloc[:, 1:].reset_index(drop=True)

    # Column index 1 (currently unnamed) should be the formula/Symbol column
    # (DEW 2019 and 2024 both follow this structure)
    if df.columns[1] != "Symbol":
        df = df.rename(columns={df.columns[1]: "Symbol"})

    return df


def load_gas_table(xls_path: Path) -> pd.DataFrame:
    """
    Load 'Gas Table' from a DEW Excel (2019 or 2024).

    Structure:
      - row 2: header names (Chemical, ΔGfo, ΔHfo, So, a, b x 103, c x 10-5, T)
      - row 3: units
      - data start at row 4
    """
    raw = pd.read_excel(xls_path, sheet_name="Gas Table", header=None)

    header_names = raw.iloc[2]
    df = raw.iloc[4:].copy()
    df.columns = header_names

    # Keep rows with a Chemical name
    df = df[df["Chemical"].notna()].copy()

    # Drop first dummy column
    df = df.iloc[:, 1:].reset_index(drop=True)

    # Column index 1 is the English name (ARGON, METHANE, ...)
    if df.columns[1] != "SpeciesName":
        df = df.rename(columns={df.columns[1]: "SpeciesName"})

    return df


# ---------------------------------------------------------------------------
# Build HKF aqueous species
# ---------------------------------------------------------------------------


def build_aqueous_species(aq: pd.DataFrame):
    """
    Build species dict for aqueous solutes using StandardThermoModelHKF
    with SI units, in the "new" Reaktoro YAML format.

    Expects columns like:
      Chemical, Symbol, ΔGfo , ΔHfo, So, Vo, Cpo, a1 x 10, a2 x 10-2,
      a3, a4 x 10-4, c1, c2 x 10-4, ω x 10-5, Z, Comments
    """
    species = {}

    for _, row in aq.iterrows():
        name = str(row["Chemical"]).strip()
        raw_formula = row["Symbol"]

        # Default: use Symbol from sheet
        formula = ""
        if isinstance(raw_formula, str):
            formula = raw_formula.strip()

        # Apply special-case overrides (bG,NaCl; Glycinate; H2CO3)
        if name in SPECIAL_AQUEOUS_SPECIES:
            formula_override, extra_comment = SPECIAL_AQUEOUS_SPECIES[name]
            formula = formula_override

        # Skip if we still don't have a usable formula
        if not formula or formula.lower() == "nan":
            # You can log or print here if you want to see which ones are skipped
            continue

        elements = parse_elements_from_formula(formula)

        # Charge (if available)
        charge = 0.0
        if "Z" in row and not pd.isna(row["Z"]):
            charge = float(row["Z"])

        # Basic thermodynamic data (cal → J)
        Gf = float(row["ΔGfo "]) * CAL2J
        Hf = float(row["ΔHfo"]) * CAL2J
        Sr = float(row["So"]) * CAL2J  # cal/(mol·K) → J/(mol·K)

        # HKF parameters with scaling from units row
        # a1 x 10 [cal mol-1 bar-1] → a1 [J mol-1 Pa-1]
        a1 = float(row["a1 x 10"]) * 10.0 * CAL2J / BAR2PA

        # a2 x 10-2 [cal mol-1] → a2 [J mol-1]
        a2 = float(row["a2 x 10-2"]) * 1.0e-2 * CAL2J

        # a3 [cal K mol-1 bar-1] → a3 [(J·K) mol-1 Pa-1]
        a3 = float(row["a3"]) * CAL2J / BAR2PA

        # a4 x 10-4 [cal K mol-1] → a4 [(J·K) mol-1]
        a4 = float(row["a4 x 10-4"]) * 1.0e-4 * CAL2J

        # c1 [cal mol-1 K-1] → J mol-1 K-1
        c1 = float(row["c1"]) * CAL2J

        # c2 x 10-4 [cal K mol-1] → (J·K) mol-1
        c2 = float(row["c2 x 10-4"]) * 1.0e-4 * CAL2J

        # ω x 10-5 [cal mol-1] → wref [J mol-1]
        wref = float(row["ω x 10-5"]) * 1.0e-5 * CAL2J

        # Base comment from the DEW sheet (if present)
        base_comment = ""
        if "Comments" in row and not pd.isna(row["Comments"]):
            base_comment = str(row["Comments"]).strip()

        # Add extra comment if this is a special species
        extra_comment = None
        if name in SPECIAL_AQUEOUS_SPECIES:
            _, extra_comment = SPECIAL_AQUEOUS_SPECIES[name]

        comment_field = None
        if base_comment and extra_comment:
            comment_field = f"{base_comment} | {extra_comment}"
        elif extra_comment:
            comment_field = extra_comment
        elif base_comment:
            comment_field = base_comment

        species_block = {
            "Name": name,
            "Formula": formula,
            "Elements": elements,
            "Charge": float(charge),
            "AggregateState": "Aqueous",
            "StandardThermoModel": {
                "HKF": {
                    "Gf": Gf,
                    "Hf": Hf,
                    "Sr": Sr,
                    "a1": a1,
                    "a2": a2,
                    "a3": a3,
                    "a4": a4,
                    "c1": c1,
                    "c2": c2,
                    "wref": wref,
                    "charge": float(charge),
                    # Optional: "Tmax": 1000.0,
                }
            },
        }

        if comment_field:
            # Not used by Reaktoro, but useful documentation
            species_block["Comment"] = comment_field

        species_key = formula  # use the formula as YAML key
        species[species_key] = species_block

    return species


# ---------------------------------------------------------------------------
# Build ideal-gas species (Cp°(T) = a + bT + cT²)
# ---------------------------------------------------------------------------


def build_gas_species(gas: pd.DataFrame):
    """
    Build gas species using StandardThermoModelNasa, with Cp°(T) = a + bT + cT²
    implemented exactly via a NASA polynomial.

    From DEW 'Gas Table':
      - ΔGfo , ΔHfo, So in cal-based units
      - 'a'          [cal mol-1]
      - 'b x 103'    [cal mol-1 K-1]   (scale 1e3)
      - 'c x 10-5'   [cal mol-1 K2]    (scale 1e-5)
      - 'T'          [K] (max T)
    """
    species = {}

    for _, row in gas.iterrows():
        chem_field = str(row["Chemical"]).strip()  # e.g. 'Ar,g', 'CO2,g'
        human_name = str(row["SpeciesName"]).strip()  # 'ARGON', 'CARBON-DIOXIDE'

        if not chem_field or chem_field.lower() == "nan":
            continue

        # e.g. 'Ar,g' -> formula 'Ar'
        formula = chem_field.split(",")[0].strip()
        name = re.sub(r",g$", "(g)", chem_field)

        elements = parse_elements_from_formula(formula)
        charge = 0.0

        # Thermo data (cal → J)
        Gf = float(row["ΔGfo "]) * CAL2J
        Hf = float(row["ΔHfo"]) * CAL2J
        S = float(row["So"]) * CAL2J

        # Cp polynomial from table, with their scaling
        a_cal = float(row["a"])  # cal mol-1
        b_cal = float(row["b x 103"])  # cal mol-1 K-1 (scaled by 1e3)
        c_cal = float(row["c x 10-5"])  # cal mol-1 K2 (scaled by 1e-5)
        Tmax = float(row["T"])  # max T [K] for fit

        # Define:
        #   Cp(T) [cal/(mol·K)] = a_cal + (b_cal * 1e-3)*T + (c_cal * 1e-5)*T²
        a_J = a_cal * CAL2J
        b_JK = b_cal * 1.0e-3 * CAL2J
        c_JK2 = c_cal * 1.0e-5 * CAL2J

        nasa_params = nasa_from_cp_poly(
            a_J=a_J,
            b_JK=b_JK,
            c_JK2=c_JK2,
            Hf_J=Hf,
            S_J=S,
            Tmin=TREF,
            Tmax=Tmax,
            label=name,
        )

        species_block = {
            "Name": name,
            "Formula": formula,
            "Elements": elements,
            "Charge": float(charge),
            "AggregateState": "Gas",
            "StandardThermoModel": {"Nasa": nasa_params},
        }

        # Use 'Formula(g)' as key to avoid collision with aqueous/solid forms
        species_key = f"{formula}(g)"
        species[species_key] = species_block

    return species


# ---------------------------------------------------------------------------
# Build YAML databases for a given Excel file
# ---------------------------------------------------------------------------


def build_yaml_for_excel(xls_path: str, tag: str):
    """
    For a DEW Excel file, generate two YAML databases:
      - {tag}-aqueous.yaml  (all aqueous species, HKF)
      - {tag}-gas.yaml      (all gas species, NASA ideal-gas Cp)
    """
    xls_path = Path(xls_path)

    print(f"Processing: {xls_path} -> tag '{tag}'")

    aq_df = load_aqueous_table(xls_path)
    gas_df = load_gas_table(xls_path)

    aq_species = build_aqueous_species(aq_df)
    gas_species = build_gas_species(gas_df)

    aq_db = {"Species": aq_species}
    gas_db = {"Species": gas_species}

    aq_out = f"{tag}-aqueous.yaml"
    gas_out = f"{tag}-gas.yaml"

    with open(aq_out, "w", encoding="utf-8") as f:
        yaml.safe_dump(aq_db, f, sort_keys=False, allow_unicode=True)
    with open(gas_out, "w", encoding="utf-8") as f:
        yaml.safe_dump(gas_db, f, sort_keys=False, allow_unicode=True)

    print(f"  -> wrote {aq_out} with {len(aq_species)} aqueous species")
    print(f"  -> wrote {gas_out} with {len(gas_species)} gas species")


if __name__ == "__main__":
    # Adjust paths if needed; assumes script is in same folder as the Excel files
    build_yaml_for_excel("Latest_DEW2024.xlsm", tag="dew2024")
    build_yaml_for_excel("DEW_2019.xlsm", tag="dew2019")
