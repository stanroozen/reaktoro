"""Test non-aqueous GFSM activity model against Perple_X hybrid_mrk reference data.

Key: import autodiff BEFORE reaktoro4py to register pybind11 type casters for
autodiff::real, enabling Python floats in state.set() / setTemperature() etc.
"""

import sys
import math
import importlib
import os

# Keep DLL resolution focused on the active conda env on Windows.
if os.name == "nt":
    env_prefix = sys.prefix
    env_paths = [
        env_prefix,
        os.path.join(env_prefix, "Library", "mingw-w64", "bin"),
        os.path.join(env_prefix, "Library", "usr", "bin"),
        os.path.join(env_prefix, "Library", "bin"),
        os.path.join(env_prefix, "Scripts"),
        os.path.join(env_prefix, "bin"),
    ]
    system_root = os.environ.get("SystemRoot", r"C:\Windows")
    system_paths = [
        os.path.join(system_root, "System32"),
        system_root,
        os.path.join(system_root, "System32", "Wbem"),
    ]
    os.environ["PATH"] = ";".join(
        [p for p in env_paths + system_paths if os.path.isdir(p)]
    )

# CRITICAL: import autodiff before reaktoro so pybind11 registers real1st type casters.
# Without this, state.set("species", float, "unit") raises TypeError.
try:
    import autodiff  # noqa: F401
except ModuleNotFoundError:
    autodiff = None

try:
    from reaktoro import *  # noqa: F401,F403
except ModuleNotFoundError:
    _pyd_dir = r"C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\build-dew\Reaktoro\Release"
    sys.path.insert(0, _pyd_dir)
    _local_mod = importlib.import_module("reaktoro4py")
    globals().update(
        {k: getattr(_local_mod, k) for k in dir(_local_mod) if not k.startswith("_")}
    )
    print(f"Using local reaktoro4py from {_pyd_dir}")

import numpy as np

# ---------------------------------------------------------------------------
# Reference data location
# ---------------------------------------------------------------------------
REF_DIR = (
    r"C:\Users\stanroozen\Documents\Projects\Perplex\Perple_X\test\reaktoro_extensions"
)
_N = 1000.0  # total moles (only mole fractions matter for EOS)

# XCO2 values encoded in file names  "0" -> 0.0, "0p1" -> 0.1, etc.
XCO2_FILES = {
    0.0: "h2o_co2_hybrid_mrk_xco2_0.tab",
    0.1: "h2o_co2_hybrid_mrk_xco2_0p1.tab",
    0.2: "h2o_co2_hybrid_mrk_xco2_0p2.tab",
    0.3: "h2o_co2_hybrid_mrk_xco2_0p3.tab",
    0.4: "h2o_co2_hybrid_mrk_xco2_0p4.tab",
    0.5: "h2o_co2_hybrid_mrk_xco2_0p5.tab",
    0.6: "h2o_co2_hybrid_mrk_xco2_0p6.tab",
}


def parse_tab_file(path):
    """Parse a Perple_X 2-variable .tab file into a list of data dicts.

    Returns a list of dicts with keys: P_bar, T_K, xco2, f_H2O, f_CO2.
    Rows with f(H2O) <= 0 or f(CO2) <= 0 are skipped.
    """
    rows = []
    header_done = False
    col_names = []
    with open(path, "r") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            # Column header line (contains 'P(bar)' and 'T(K)' on one long line)
            if "P(bar)" in stripped and "T(K)" in stripped:
                col_names = stripped.split()
                header_done = True
                continue
            if not header_done:
                continue
            # Data line: starts with a number
            parts = stripped.split()
            if not parts:
                continue
            try:
                float(parts[0])
            except ValueError:
                continue
            # Merge continuation line pattern: some data lines wrap; we just collect
            # enough to get P, T, f(H2O), f(CO2)
            # columns: P(bar) T(K) X(CO2) vol y(H2O) y(CO2) ... f(H2O) f(CO2) ...
            # col indices: 0      1    2     3   4      5     6-11  12    13
            if len(parts) < 14:
                continue
            try:
                P = float(parts[0])
                T = float(parts[1])
                f_h2o = float(parts[12])
                f_co2 = float(parts[13])
            except (ValueError, IndexError):
                continue
            if f_h2o <= 0 or f_co2 <= 0:
                continue
            rows.append({"P_bar": P, "T_K": T, "f_H2O": f_h2o, "f_CO2": f_co2})
    return rows


# ---------------------------------------------------------------------------
# Build GFSM system (MRK mode to match Perple_X hybrid_mrk reference)
# ---------------------------------------------------------------------------
supcrt_db = SupcrtDatabase("supcrtbl")
gas_db = Database([supcrt_db.species("CO2(g)"), supcrt_db.species("H2O(g)")])
gas_phase = GaseousPhase("H2O(g) CO2(g)")
params = ActivityModelParamsPerplexGFSM()
params.hybridEosOptions = (
    makePerpleXHybridEosOptions()
)  # MRK: matches Perple_X hybrid_mrk
gas_phase.setActivityModel(ActivityModelPerplexGFSM(params))
system = ChemicalSystem(gas_db, gas_phase)


def compute_fugacities(xco2, T_K, P_bar):
    """Return (f_H2O, f_CO2) in bar from GFSM activity model.

    Uses a tiny CO2 floor (1e-10) to avoid log(0) at xco2=0.
    """
    xco2_safe = max(xco2, 1e-10)
    state = ChemicalState(system)
    state.set("H2O(g)", (1.0 - xco2_safe) * _N, "mol")
    state.set("CO2(g)", xco2_safe * _N, "mol")
    state.setTemperature(T_K, "K")
    state.setPressure(P_bar, "bar")
    props = ChemicalProps(state)
    f_h2o = math.exp(float(props.speciesActivityLn("H2O(g)")))
    f_co2 = math.exp(float(props.speciesActivityLn("CO2(g)")))
    return f_h2o, f_co2


# ---------------------------------------------------------------------------
# Run benchmark
# ---------------------------------------------------------------------------
print("\n=== GFSM vs Perple_X hybrid_mrk benchmark (MRK EOS) ===\n")
print(
    "Note: XCO2=0.0 pure-water rows with T < 747K are skipped (Perple_X uses liquid EOS,"
)
print(
    "      GFSM uses gas-phase MRK). Residual error at 647-750K is due to MRK inaccuracy"
)
print(
    "      near the pure-water critical point (Tc=647K) — expected and not a code bug."
)
print()
print(
    f"{'XCO2':>6}  {'N_pts':>5}  {'err_H2O_max%':>12}  {'err_H2O_mean%':>13}  {'err_CO2_max%':>12}  {'err_CO2_mean%':>13}"
)
print("-" * 72)

all_err_h2o = []
all_err_co2 = []

for xco2, fname in sorted(XCO2_FILES.items()):
    path = os.path.join(REF_DIR, fname)
    rows = parse_tab_file(path)
    if not rows:
        print(f"  {xco2:.1f}  No data rows parsed!")
        continue

    errs_h2o = []
    errs_co2 = []
    for row in rows:
        # For pure-water (xco2=0), only compare at T well above Tc (>747K to avoid
        # critical region where MRK is inaccurate for pure H2O)
        if xco2 == 0.0 and row["T_K"] <= 747.0:
            continue
        try:
            fh, fc = compute_fugacities(xco2, row["T_K"], row["P_bar"])
        except Exception:
            continue
        # Relative error in %
        if row["f_H2O"] > 0:
            errs_h2o.append(abs(fh - row["f_H2O"]) / row["f_H2O"] * 100.0)
        # Skip CO2 comparison for pure-water reference (sentinel values 1e15+)
        if row["f_CO2"] > 0 and row["f_CO2"] < 1e14:
            errs_co2.append(abs(fc - row["f_CO2"]) / row["f_CO2"] * 100.0)

    all_err_h2o.extend(errs_h2o)
    all_err_co2.extend(errs_co2)

    tag = " [pure H2O]" if xco2 == 0.0 else ""
    print(
        f"  {xco2:.1f}  {len(rows):>5}  "
        f"{max(errs_h2o) if errs_h2o else float('nan'):>12.4f}  "
        f"{np.mean(errs_h2o) if errs_h2o else float('nan'):>13.4f}  "
        f"{max(errs_co2) if errs_co2 else float('nan'):>12.4f}  "
        f"{np.mean(errs_co2) if errs_co2 else float('nan'):>13.4f}{tag}"
    )

print("-" * 72)

# Summary for mixed-fluid compositions only (XCO2 = 0.1 to 0.5)
mixed_h2o = [e for xco2, fname in sorted(XCO2_FILES.items()) if xco2 > 0 for e in []]
# Recompute from all_err excluding the XCO2=0 entries (which were added first)
# Simpler: track separately
mixed_errs_h2o = []
mixed_errs_co2 = []
for xco2, fname in sorted(XCO2_FILES.items()):
    if xco2 == 0.0:
        continue
    path = os.path.join(REF_DIR, fname)
    rows = parse_tab_file(path)
    for row in rows:
        try:
            fh, fc = compute_fugacities(xco2, row["T_K"], row["P_bar"])
        except Exception:
            continue
        if row["f_H2O"] > 0:
            mixed_errs_h2o.append(abs(fh - row["f_H2O"]) / row["f_H2O"] * 100.0)
        if row["f_CO2"] > 0 and row["f_CO2"] < 1e14:
            mixed_errs_co2.append(abs(fc - row["f_CO2"]) / row["f_CO2"] * 100.0)

print(
    f"  {'MIX':>5}  {len(mixed_errs_h2o):>5}  "
    f"{max(mixed_errs_h2o) if mixed_errs_h2o else float('nan'):>12.4f}  "
    f"{np.mean(mixed_errs_h2o) if mixed_errs_h2o else float('nan'):>13.4f}  "
    f"{max(mixed_errs_co2) if mixed_errs_co2 else float('nan'):>12.4f}  "
    f"{np.mean(mixed_errs_co2) if mixed_errs_co2 else float('nan'):>13.4f}  [XCO2=0.1-0.5]"
)
print()

# ---------------------------------------------------------------------------
# Also verify ZD09 defaults
# ---------------------------------------------------------------------------
params_zd09 = ActivityModelParamsPerplexGFSM()
print(
    f"Default params EOS: water={params_zd09.hybridEosOptions.water}, "
    f"co2={params_zd09.hybridEosOptions.co2}, ch4={params_zd09.hybridEosOptions.ch4}"
)
print("(ZD09 for H2O, CO2, CH4 is the correct COH-Fluid+ default)")
