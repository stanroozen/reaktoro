"""Compare MRK vs ZD09 water EOS in GFSM against Perple_X liquid water reference."""

import sys, os, math, importlib

if os.name == "nt":
    ep = sys.prefix
    env_paths = [
        ep,
        os.path.join(ep, "Library", "mingw-w64", "bin"),
        os.path.join(ep, "Library", "usr", "bin"),
        os.path.join(ep, "Library", "bin"),
        os.path.join(ep, "Scripts"),
        os.path.join(ep, "bin"),
    ]
    sr = os.environ.get("SystemRoot", r"C:\Windows")
    os.environ["PATH"] = ";".join(
        [p for p in env_paths + [os.path.join(sr, "System32"), sr] if os.path.isdir(p)]
    )
import autodiff  # noqa: F401

try:
    from reaktoro import *  # noqa
except ModuleNotFoundError:
    _pyd = r"C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\build-dew\Reaktoro\Release"
    sys.path.insert(0, _pyd)
    m = importlib.import_module("reaktoro4py")
    globals().update({k: getattr(m, k) for k in dir(m) if not k.startswith("_")})

db = SupcrtDatabase("supcrtbl")
gd = Database([db.species("CO2(g)"), db.species("H2O(g)")])


def make_system(water_eos):
    gp = GaseousPhase("H2O(g) CO2(g)")
    p = ActivityModelParamsPerplexGFSM()
    opts = makePerpleXHybridEosOptions()
    opts.water = water_eos
    p.hybridEosOptions = opts
    gp.setActivityModel(ActivityModelPerplexGFSM(p))
    return ChemicalSystem(gd, gp)


sys_mrk = make_system(PerpleXWaterEos.Mrk)
sys_zd09 = make_system(PerpleXWaterEos.ZhangDuan09)

N = 1000.0
xco2 = 1e-10  # near-pure water

# Perple_X hybrid_mrk reference at XCO2=0, P=1000 bar
# T below Tc: Perple_X uses liquid water EOS; T above Tc: gas-phase MRK
refs = [
    (373.15, 1000, 1.805024),  # liquid, T < Tc
    (473.15, 1000, 24.00007),  # liquid, T < Tc
    (573.15, 1000, 105.5029),  # liquid, T < Tc
    (623.15, 1000, 175.4780),  # just below Tc=647K
    (673.15, 1000, 262.3274),  # just above Tc
    (773.15, 1000, 460.8531),  # supercritical
    (873.15, 1000, 640.3860),  # supercritical
    (1073.15, 1000, 857.6867),  # supercritical
    (1573.15, 1000, 1013.745),  # supercritical
]

print("GFSM water EOS comparison vs Perple_X hybrid_mrk (XCO2=0, P=1000 bar)")
print(
    f"{'T(K)':>7}  {'f_ref':>10}  {'f_MRK':>10}  {'err_MRK%':>9}  {'f_ZD09':>10}  {'err_ZD09%':>10}  note"
)
print("-" * 78)
for T_K, P, f_ref in refs:
    results = {}
    for sys_, key in [(sys_mrk, "mrk"), (sys_zd09, "zd09")]:
        st = ChemicalState(sys_)
        st.set("H2O(g)", (1.0 - xco2) * N, "mol")
        st.set("CO2(g)", xco2 * N, "mol")
        st.setTemperature(T_K, "K")
        st.setPressure(P, "bar")
        pr = ChemicalProps(st)
        results[key] = math.exp(float(pr.speciesActivityLn("H2O(g)")))
    f_mrk = results["mrk"]
    f_zd09 = results["zd09"]
    err_mrk = (f_mrk - f_ref) / f_ref * 100
    err_zd09 = (f_zd09 - f_ref) / f_ref * 100
    note = "liquid" if T_K < 647 else ("near-crit" if T_K < 700 else "supercrit")
    print(
        f"{T_K:>7.2f}  {f_ref:>10.4f}  {f_mrk:>10.4f}  {err_mrk:>+9.2f}%  {f_zd09:>10.4f}  {err_zd09:>+10.2f}%  {note}"
    )

print()
print("ZD09 range: 240-2000 K, 0-10 GPa (Zhang & Duan 2009)")
print("Tc(H2O) = 647.1 K, Pc(H2O) = 220.6 bar")
