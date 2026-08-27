"""
Cross-database quartz thermodynamic data comparison.

Checks H°, G°, S°, V° for quartz ("q") across all available PerpleX JSON databases
and the Reaktoro SUPCRT/DEW databases to detect convention or unit inconsistencies.

Run with:
    python DEW_Experimental_Benchmark/DatabaseTests/quartz_thermo_comparison.py
"""

import sys
import os
import json
import glob
from pathlib import Path

# ──────────────────────────────────────────────────────────────────────────────
# 1. Scan all PerpleX JSON databases for quartz "q" entry
# ──────────────────────────────────────────────────────────────────────────────

DB_DIR = Path(__file__).parent.parent.parent / "embedded" / "databases" / "perplex"
json_files = sorted(DB_DIR.glob("*-reaktoro.json"))

print("=" * 100)
print(
    f"{'DATABASE':<45} {'H (kJ/mol)':>13} {'GH (kJ/mol)':>13} {'S0 (J/mol/K)':>13} {'V0 (cm3/mol)':>13} {'Notes'}"
)
print("=" * 100)

perplex_rows = []

for jf in json_files:
    db_name = jf.stem.replace("-reaktoro", "")
    try:
        with open(jf, encoding="utf-8") as f:
            data = json.load(f)
    except Exception as e:
        print(f"  [SKIP] {db_name}: cannot load JSON ({e})")
        continue

    species = data.get("Species", {})
    if "q" not in species:
        print(f"  [SKIP] {db_name:<41} — 'q' not found")
        continue

    entry = species["q"]
    meta = entry.get("Metadata", {})
    params = meta.get("PerpleX_Params", {})

    H = params.get("H", float("nan")) / 1000.0  # J → kJ
    GH = params.get("GH", float("nan")) / 1000.0  # J → kJ
    S0 = params.get("S0", float("nan"))  # J/mol/K
    V0 = params.get("V0", float("nan"))  # cm³/mol  (PerpleX native)
    mode = meta.get("GfConversionMode", "—")

    print(f"  {db_name:<43} {H:>13.2f} {GH:>13.2f} {S0:>13.3f} {V0:>13.4f}  [{mode}]")
    perplex_rows.append(dict(db=db_name, H=H, GH=GH, S0=S0, V0=V0, mode=mode))

print()

# ──────────────────────────────────────────────────────────────────────────────
# 2. Reaktoro built-in databases: SUPCRTBL, Phreeqc, DEW
# ──────────────────────────────────────────────────────────────────────────────

try:
    RKTORO_PATH = str(
        Path(__file__).parent.parent.parent / "build" / "Reaktoro" / "Release"
    )
    sys.path.insert(0, RKTORO_PATH)
    import reaktoro as rkt

    HAS_REAKTORO = True
    print("Reaktoro imported OK.")
except Exception as e:
    HAS_REAKTORO = False
    print(f"[WARNING] Could not import reaktoro: {e}")
    print("Skipping SUPCRT/DEW/Phreeqc database checks.")

if HAS_REAKTORO:

    def species_thermo_at_reference(species_obj):
        """Return (G0_kJ, H0_kJ, S0_J, V0_cm3) at 298.15 K, 1 bar via Reaktoro props."""
        T = 298.15  # K
        P = 1e5  # Pa  (1 bar)
        props = species_obj.standardThermoProps(T, P)
        G0 = props.G0 / 1000.0  # J → kJ
        H0 = props.H0 / 1000.0  # J → kJ
        S0 = (H0 - G0) * 1000.0 / T  # derived S° = (H-G)/T  in J/mol/K
        V0 = props.V0 * 1e6  # m³ → cm³
        return G0, H0, S0, V0

    rkt_dbs = {
        "supcrtbl": lambda: rkt.SupcrtDatabase("supcrtbl"),
        "supcrt07": lambda: rkt.SupcrtDatabase("supcrt07"),
        "supcrt16": lambda: rkt.SupcrtDatabase("supcrt16"),
        "dew2024-aqueous": lambda: rkt.DEWDatabase("dew2024-aqueous"),
        "dew2019-aqueous": lambda: rkt.DEWDatabase("dew2019-aqueous"),
        "phreeqc.dat": lambda: rkt.PhreeqcDatabase("phreeqc.dat"),
    }

    # Quartz names to try per database
    QUARTZ_NAMES = ["Quartz", "q", "SiO2", "Quartz(alpha)", "SiO2(a)", "qtz"]

    print()
    print("=" * 100)
    print(
        f"{'DATABASE':<45} {'G0 (kJ/mol)':>13} {'H0 (kJ/mol)':>13} {'S0 (J/mol/K)':>13} {'V0 (cm3/mol)':>13}"
    )
    print("=" * 100)

    for db_label, db_factory in rkt_dbs.items():
        try:
            db = db_factory()
        except Exception as e:
            print(f"  [SKIP] {db_label:<41} — cannot load ({e})")
            continue

        found = False
        for qname in QUARTZ_NAMES:
            try:
                sp = db.species(qname)
                G0, H0, S0_derived, V0 = species_thermo_at_reference(sp)
                print(
                    f"  {db_label:<43} {G0:>13.2f} {H0:>13.2f} {S0_derived:>13.3f} {V0:>13.4f}  [name={qname!r}]"
                )
                found = True
                break
            except Exception:
                continue

        if not found:
            # List available SiO2 species in this DB
            try:
                names = [
                    s.name()
                    for s in db.species()
                    if "Si" in s.formula().str() and "O" in s.formula().str()
                ]
                print(
                    f"  [SKIP] {db_label:<41} — quartz not found. SiO2 species: {names[:8]}"
                )
            except Exception:
                print(f"  [SKIP] {db_label:<41} — quartz not found")

    print()

# ──────────────────────────────────────────────────────────────────────────────
# 3. Summary: check spread across PerpleX databases
# ──────────────────────────────────────────────────────────────────────────────

if perplex_rows:
    import statistics

    Hs = [r["H"] for r in perplex_rows if r["H"] == r["H"]]  # NaN-safe
    GHs = [r["GH"] for r in perplex_rows if r["GH"] == r["GH"]]
    S0s = [r["S0"] for r in perplex_rows if r["S0"] == r["S0"]]
    V0s = [r["V0"] for r in perplex_rows if r["V0"] == r["V0"]]

    print("=" * 100)
    print("PERPLEX JSON SPREAD SUMMARY for quartz 'q'")
    print(
        f"  H  (kJ/mol): min={min(Hs):.2f}  max={max(Hs):.2f}  range={max(Hs) - min(Hs):.3f}"
    )
    print(
        f"  GH (kJ/mol): min={min(GHs):.2f}  max={max(GHs):.2f}  range={max(GHs) - min(GHs):.3f}"
    )
    print(
        f"  S0 (J/K/mol):min={min(S0s):.3f}  max={max(S0s):.3f}  range={max(S0s) - min(S0s):.4f}"
    )
    print(
        f"  V0 (cm³/mol):min={min(V0s):.4f}  max={max(V0s):.4f}  range={max(V0s) - min(V0s):.5f}"
    )
    print()

    # Flag any outliers (>0.5 kJ/mol deviation from median for energy terms)
    H_med = statistics.median(Hs)
    GH_med = statistics.median(GHs)
    print("  Outliers (|deviation| > 0.5 kJ/mol from median):")
    for r in perplex_rows:
        dH = abs(r["H"] - H_med)
        dGH = abs(r["GH"] - GH_med)
        if dH > 0.5 or dGH > 0.5:
            print(f"    {r['db']:<43} ΔH={dH:.3f}  ΔGH={dGH:.3f} kJ/mol  [{r['mode']}]")
    print("=" * 100)
