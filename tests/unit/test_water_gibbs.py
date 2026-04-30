import sys

sys.path.insert(
    0,
    r"C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\build-msvc\python\package\build\lib.win-amd64-cpython-312",
)

from reaktoro import *

# Test conditions
T_C = 300.0  # °C
P_kb = 5.0  # kilobar
T_K = T_C + 273.15  # K
P_Pa = P_kb * 1e8  # Pa

print(f"Testing water Gibbs at T={T_C}°C, P={P_kb} kb")
print(f"T = {T_K} K, P = {P_Pa} Pa = {P_Pa / 1e5} bar")
print()

# Create water state with high-precision integration
from reaktoro.core import WaterStateOptions, WaterGibbsModelOptions

opts = WaterStateOptions()
opts.thermo.eosModel = 2  # ZhangDuan2005
opts.computeGibbs = True
opts.gibbs.model = 1  # DewIntegral
opts.gibbs.integrationSteps = 5000
opts.gibbs.useExcelIntegration = False

ws = waterState(T_K, P_Pa, opts)

G_J_per_mol = ws.gibbs
G_cal_per_mol = G_J_per_mol / 4.184

print(f"G (high-precision, 5000 steps) = {G_J_per_mol:.4f} J/mol")
print(f"G (high-precision, 5000 steps) = {G_cal_per_mol:.4f} cal/mol")
print()

# Test with Excel mode
opts.gibbs.useExcelIntegration = True
ws_excel = waterState(T_K, P_Pa, opts)
G_excel_J = ws_excel.gibbs
G_excel_cal = G_excel_J / 4.184

print(f"G (Excel mode, ~500 steps)     = {G_excel_J:.4f} J/mol")
print(f"G (Excel mode, ~500 steps)     = {G_excel_cal:.4f} cal/mol")
print()

# Truth from Excel
G_truth_cal = -60673.6416951689
print(f"G (Excel truth from CSV)       = {G_truth_cal:.4f} cal/mol")
print()

print("Differences:")
print(
    f"  High-precision vs Truth: {G_cal_per_mol - G_truth_cal:.4f} cal/mol ({100 * (G_cal_per_mol - G_truth_cal) / G_truth_cal:.4f}%)"
)
print(
    f"  Excel mode vs Truth:     {G_excel_cal - G_truth_cal:.4f} cal/mol ({100 * (G_excel_cal - G_truth_cal) / G_truth_cal:.4f}%)"
)
