"""
Calcite solubility analysis — DEW backend, CO2 as dissolved aqueous solute.

CO2 is included as aqueous CO2_aq at a fixed initial molality rather than as a
gas-phase solvent with an EOS. No GaseousPhase is added; the solvent is pure H2O.

Species included via DEW element filter (H, O, Ca, C), after excluding CaO_aq,
CO_aq, H2_aq, O2_aq (redox gases irrelevant to carbonate chemistry):
  WATER,AQ  H+  OH-  CO2_aq  CO3-2  HCO3-  H2CO3_aq
  CaCO3_aq  Ca(HCO3)+  Ca(OH)+  Ca+2

Runs calcite_solubility_analysis.py with
  --backend DEW --fluid H2OCO2aq --co2aq-molality 6.168

Outputs (in Mineral_Solubilities/):
  calcite_solubility_comparison_high_P_dew24_DEW_H2OCO2AQ.png
  calcite_solubility_residuals_dew24_DEW_H2OCO2AQ.png
  calcite_speciation_dew24_DEW_H2OCO2AQ.png
"""

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

sys.argv = [
    sys.argv[0],
    "--backend",
    "DEW",
    "--fluid",
    "H2OCO2aq",
    "--co2aq-molality",
    "6.168",
]

import calcite_solubility_analysis as _cs  # noqa: E402

_cs.main()
