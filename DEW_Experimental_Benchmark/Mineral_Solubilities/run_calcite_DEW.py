"""
Calcite solubility analysis — DEW backend, Zhang-Duan 2005 EOS.

Runs calcite_solubility_analysis.py with --backend DEW.

Outputs (in Mineral_Solubilities/):
  calcite_solubility_comparison_high_P_dew24_DEW.png
  calcite_solubility_residuals_dew24_DEW.png
  calcite_speciation_dew24_DEW.png
"""

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

sys.argv = [sys.argv[0], "--backend", "DEW"]

import calcite_solubility_analysis as _cs  # noqa: E402

_cs.main()
