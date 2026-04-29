"""
Calcite solubility analysis â€” PerplexDEW backend, Extended Debye-HÃ¼ckel (HB) activity model.

Runs calcite_solubility_analysis.py with --backend PerplexDEW --dh-model ExtendedDH.
Requires a local reaktoro4py build with PerplexDEW symbols in one of:
  build/Reaktoro/Release
  build/Reaktoro/Release
  build/Reaktoro/Release

Water EOS: Zhang-Duan 2005.

Outputs (in Mineral_Solubilities/):
  calcite_solubility_comparison_high_P_dew24_PerplexDEW_ExtendedDH.png
  calcite_solubility_residuals_dew24_PerplexDEW_ExtendedDH.png
  calcite_speciation_dew24_PerplexDEW_ExtendedDH.png
"""

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

sys.argv = [sys.argv[0], "--backend", "PerplexDEW", "--dh-model", "ExtendedDH"]

import calcite_solubility_analysis as _cs  # noqa: E402

_cs.main()

