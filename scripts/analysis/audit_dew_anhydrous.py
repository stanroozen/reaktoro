"""
audit_dew_anhydrous.py
======================
Audits the DEW2024 aqueous database for species where the anhydrous
convention has been applied — i.e., species where the stored formula
has fewer H atoms than the true dissolved molecular species.

Key insight: Reaktoro parses element composition from the FORMULA field
(not the Elements YAML field). We therefore use Reaktoro's DEWDatabase
directly, so formula parsing is identical to what the solver sees.

The DEW convention projects out n H2O units:
    formula_stored = formula_true - n * H2O

Detection method:
  For "hydroxo-type" elements (Si, Al, B, Ge, Ti, Zr, Sn) where the
  true dissolved form is a fully OH-substituted M(OH)_m species:
      n_H2O = (nH_expected - nH_stored) / 2
      where nH_expected = nO_stored + max(0, -charge)

  False-positive filter: exclude "oxyanion-type" elements (P, As, Mo, W,
  S, Cr) where the formula P=O double bonds mean nH < nO even in the
  true molecular formula (e.g. H3PO4 has 3H but 4O due to P=O).

Output: table of all flagged species with formula, n_H2O, and element.
"""

import os
import sys

# --- Reaktoro bootstrap (same as benchmark scripts) ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if os.name == "nt":
    ep = sys.prefix
    import os as _os

    env_paths = [
        ep,
        _os.path.join(ep, "Library", "mingw-w64", "bin"),
        _os.path.join(ep, "Library", "bin"),
        _os.path.join(ep, "Scripts"),
    ]
    sr = _os.environ.get("SystemRoot", r"C:\Windows")
    _os.environ["PATH"] = ";".join(
        [
            p
            for p in env_paths + [_os.path.join(sr, "System32"), sr]
            if _os.path.isdir(p)
        ]
    )

for _build_pkg in [
    os.path.join(SCRIPT_DIR, "build-dew", "python", "package"),
    os.path.join(SCRIPT_DIR, "build-msvc", "python", "package"),
    os.path.join(SCRIPT_DIR, "build", "python", "package"),
]:
    _rkt_inner = os.path.join(_build_pkg, "reaktoro")
    if os.path.isdir(_rkt_inner):
        if _build_pkg not in sys.path:
            sys.path.insert(0, _build_pkg)
        os.environ["PATH"] = _rkt_inner + os.pathsep + os.environ.get("PATH", "")
        break

import autodiff  # noqa

try:
    from reaktoro import *  # noqa
except ModuleNotFoundError:
    import importlib

    for _d in [
        os.path.join(SCRIPT_DIR, "build-dew", "Reaktoro", "Release"),
        os.path.join(SCRIPT_DIR, "build-msvc", "Reaktoro", "Release"),
    ]:
        if os.path.isdir(_d):
            sys.path.insert(0, _d)
            _m = importlib.import_module("reaktoro4py")
            globals().update(
                {k: getattr(_m, k) for k in dir(_m) if not k.startswith("_")}
            )
            break

# ---------------------------------------------------------------------------
# Classification tables
# ---------------------------------------------------------------------------

# Elements where the DEW anhydrous convention applies:
# true dissolved species = fully OH-substituted M(OH)_m
# ALL oxygens should carry a hydrogen in the true hydrated form.
HYDROXO_ELEMENTS = {"Si", "Al", "B", "Ge", "Ti", "Zr", "Sn"}

# Elements where nH < nO is correct even in the true formula (P=O double bonds):
# H3PO4 has 4O but only 3H; MoO4^2- has no H at all. These are FALSE POSITIVES
# for the hydroxo heuristic — do NOT flag them.
OXYANION_ELEMENTS = {"P", "As", "Mo", "W", "S", "Cr", "V", "Nb", "Se", "Te"}

# Known true hydrated formulas. Key = sp.name() as returned by Reaktoro
# (= YAML 'Name:' field with commas replaced by underscores).
# Value = (nH_true, nO_true); n_H2O = (nH_true - nH_stored) / 2.
#
# Derivation for each species:
#   neutral M(OH)_k :  nH_true = nO_true = k
#   "excess-OH" anions (Al(OH)4-, B(OH)4-): nH_true = nO_true = k (all O carry H)
#   "deprotonated oxyacid" anions (H3SiO4-): nH_true = k-|z|, nO_true = k
KNOWN_TRUE_FORMULA = {
    # Si mononuclear
    "SiO2_aq": (4, 4),  # H4SiO4,  n=2  (SiO2 + 2H2O -> H4SiO4)
    "HSiO3-": (3, 4),  # H3SiO4-, n=1  (HSiO3- + H2O -> H3SiO4-)
    # Si oligomers (each SiO2 unit needs 2 H2O)
    "Si2O4_aq": (8, 8),  # 2xH4SiO4, n=4
    "Si3O6_aq": (12, 12),  # 3xH4SiO4, n=6
    # Si + metal cation complexes (silicate ligand same as HSiO3- part, n=1)
    "Ca(HSiO3)+": (3, 4),  # Mg2+ + H3SiO4-, n=1
    "Fe(HSiO3)+": (3, 4),
    "Mg(HSiO3)+": (3, 4),
    "NaHSiO3_aq": (3, 4),  # Na+ ion pair with H3SiO4-, n=1
    # Mixed Al-Si (Al part: 2 H2O; Si part: 2 H2O -> total n=4)
    "AlO2(SiO2)-": (8, 8),  # Al(OH)4- . H4SiO4 complex, n=4
    # Al mononuclear (excess-OH type: Al(OH)4- has all 4 O with H)
    "AlO2-": (4, 4),  # Al(OH)4-, n=2
    "HAlO2_aq": (3, 3),  # Al(OH)3,  n=1
    # B species (excess-OH type: B(OH)4- has all 4 O with H)
    "BO2-": (4, 4),  # B(OH)4-, n=2
    "BO(OH)_aq": (3, 3),  # B(OH)3,   n=1
    # Ti, Zr (neutral M(OH)4 type)
    "TIO2_aq": (4, 4),  # Ti(OH)4, n=2  (YAML Name: TIO2,aq -- Reaktoro capitalises)
    "ZrO2_aq": (4, 4),  # Zr(OH)4, n=2
    # Ge (if present in database)
    "GeO2_aq": (4, 4),  # H4GeO4, n=2
    "HGeO3-": (3, 4),  # H3GeO4-, n=1
}

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def get_elements(species):
    """Return dict {symbol: amount} for a Reaktoro Species object."""
    d = {}
    for item in species.elements():
        if isinstance(item, tuple):
            elem, coeff = item
            d[elem.symbol()] = coeff
        else:
            # ElementAmount or similar object
            d[item.element().symbol()] = item.coefficient()
    return d


def n_h2o_heuristic(nH, nO, charge):
    """
    For a fully-hydroxylated M(OH)_k species, each O carries one H.
    The number of H2O units stripped = nO_stored - nH_stored.
    Derivation: M(OH)_k - n*H2O -> MO_n(OH)_{k-2n}
      nH_stored = k - 2n,  nO_stored = k - n
      => n = nO_stored - nH_stored
    Valid for neutral M(OH)_k and 'excess-OH' anions like Al(OH)4-, B(OH)4-.
    NOT valid for 'deprotonated oxyacid' anions (HSiO3-, NaHSiO3) -- those
    must be in KNOWN_TRUE_FORMULA.
    """
    return float(nO - nH)


# ---------------------------------------------------------------------------
# Main audit
# ---------------------------------------------------------------------------


def main():
    print("Loading DEW2024 database via Reaktoro…")
    db = DEWDatabase("dew2024-aqueous")
    all_species = db.species()

    results = []

    for sp in all_species:
        name = sp.name()
        formula = sp.formula().str()
        charge = sp.charge()
        elems = get_elements(sp)

        nH = elems.get("H", 0.0)
        nO = elems.get("O", 0.0)

        # Identify which oxyacid/hydroxo elements are present
        hydroxo_present = [e for e in HYDROXO_ELEMENTS if e in elems]
        oxyanion_present = [e for e in OXYANION_ELEMENTS if e in elems]

        if not hydroxo_present:
            continue  # not an oxyacid-type metal species

        # Skip pure cations/anions with no oxygen (e.g. Al+3)
        if nO == 0:
            continue

        # --- KNOWN_TRUE lookup (highest confidence) ---
        if name in KNOWN_TRUE_FORMULA:
            nH_true, nO_true = KNOWN_TRUE_FORMULA[name]
            n = (nH_true - nH) / 2.0
            if n > 0:
                results.append(
                    {
                        "name": name,
                        "formula": formula,
                        "nH": nH,
                        "nO": nO,
                        "charge": charge,
                        "n_H2O": n,
                        "elem": ", ".join(hydroxo_present),
                        "method": "known",
                    }
                )
            continue

        # --- Heuristic: hydroxo-type elements (Si, Al, B, Ge, Ti, Zr, Sn) ---
        # Skip if also contains a false-positive oxyanion element (e.g. AlPO4)
        if oxyanion_present:
            continue

        # Skip species where the formula already contains its true H count
        # (i.e. B(OH)3 is already fully protonated — parsed formula gives H=3, O=3)
        # Heuristic: if nH >= nO (ignoring charge) -> no stripping needed
        if nH >= nO:
            continue

        n = n_h2o_heuristic(nH, nO, charge)
        if n >= 0.5:
            results.append(
                {
                    "name": name,
                    "formula": formula,
                    "nH": nH,
                    "nO": nO,
                    "charge": charge,
                    "n_H2O": n,
                    "elem": ", ".join(hydroxo_present),
                    "method": "heuristic",
                }
            )

    # Sort by element, then n_H2O descending
    results.sort(key=lambda r: (r["elem"], -r["n_H2O"]))

    # -----------------------------------------------------------------------
    # Report
    # -----------------------------------------------------------------------
    print()
    print("=" * 95)
    print("DEW2024 Aqueous Database — Anhydrous Convention Audit")
    print("Species where the stored formula has H2O units stripped (n_H2O > 0)")
    print("=" * 95)
    print(f"\nTotal aqueous species in DEW2024: {len(list(all_species))}")
    print(f"Species requiring hydration correction: {len(results)}\n")

    hdr = f"  {'Species name':<26} {'Formula':<22} {'nH':>4} {'nO':>4} {'z':>5}  {'n_H2O':>6}  {'Element(s)':<12}  Source"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    current_elem = None
    for r in results:
        elem = r["elem"]
        if elem != current_elem:
            if current_elem is not None:
                print()
            print(f"\n  -- {elem} --")
            current_elem = elem
        n_str = (
            f"{r['n_H2O']:.1f}"
            if r["n_H2O"] != int(r["n_H2O"])
            else f"{int(r['n_H2O'])}"
        )
        flag = "*" if r["method"] == "known" else " "
        print(
            f"  {flag} {r['name']:<25} {r['formula']:<22} {r['nH']:>4.0f} {r['nO']:>4.0f} "
            f"{r['charge']:>5.1f}  {n_str:>6}  {elem:<12}  {r['method']}"
        )

    print()
    print("  * = verified against known true formula")
    print("    = flagged by hydroxo-type heuristic (nO > nH, no P=O oxyanion elements)")
    print()
    print("=" * 95)
    print("Summary by element:")
    print("=" * 95)
    by_elem = {}
    for r in results:
        by_elem.setdefault(r["elem"], []).append(r)
    for elem, sps in sorted(by_elem.items()):
        n_vals = [s["n_H2O"] for s in sps]
        print(
            f"  {elem:<12}  {len(sps):2d} species  "
            f"n_H2O: {min(n_vals):.1f} – {max(n_vals):.1f}  "
            f"({', '.join(s['name'] for s in sps)})"
        )

    print()
    print("Interpretation:")
    print("  n_H2O = number of H2O units projected out of the true formula.")
    print("  G°_hydrated(T,P) = G°_DEW(T,P) + n · G°_H2O(T,P)")
    print("  This correction must be applied to each flagged species to allow")
    print("  the Gibbs minimiser to propagate water-activity changes automatically.")
    print()
    print("NOT flagged (correct formulas already stored):")
    print("  · Simple cations/anions: Na+, Ca2+, Cl-, SO4^2-, CO3^2-, etc.")
    print("  · Oxyanion-type: PO4^3-, HPO4^2-, H2PO4-, H3PO4, HAsO4^2-,")
    print("    MoO4^2-, WO4^2- (P=O double bonds -> nH < nO is physically correct)")
    print("  · B(OH)3(0): stored formula parses to BH3O3 — already the true formula")
    print("  · B(OH)4-: stored formula parses to BH4O4 — already the true formula")


if __name__ == "__main__":
    main()
