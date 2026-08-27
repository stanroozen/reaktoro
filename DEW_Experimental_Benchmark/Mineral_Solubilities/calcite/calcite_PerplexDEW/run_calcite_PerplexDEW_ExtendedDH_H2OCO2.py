"""
Calcite solubility analysis â€” PerplexDEW backend, Extended Debye-HÃ¼ckel model, H2O-CO2 fluid.

Runs calcite_solubility_analysis.py with:
  --backend PerplexDEW --dh-model ExtendedDH --fluid H2OCO2 --xco2 0.1

Requires a local reaktoro4py build with PerplexDEW symbols in one of:
  build/Reaktoro/Release
  build/Reaktoro/Release
  build/Reaktoro/Release

Water EOS: Zhang-Duan 2005.
CO2 phase EOS: Peng-Robinson.

Outputs (in Mineral_Solubilities/):
  calcite_solubility_comparison_high_P_dew24_PerplexDEW_ExtendedDH_H2OCO2.png
  calcite_solubility_residuals_dew24_PerplexDEW_ExtendedDH_H2OCO2.png
  calcite_speciation_dew24_PerplexDEW_ExtendedDH_H2OCO2.png
"""

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

sys.argv = [
    sys.argv[0],
    "--backend",
    "PerplexDEW",
    "--dh-model",
    "ExtendedDH",
    "--fluid",
    "H2OCO2",
    "--xco2",
    "0.1",
]

import calcite_solubility_analysis as _cs  # noqa: E402

_cs.main()

