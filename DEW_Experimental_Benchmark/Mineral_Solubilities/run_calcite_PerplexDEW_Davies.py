"""
Calcite solubility analysis â€” PerplexDEW backend, Davies activity model.

Runs calcite_solubility_analysis.py with --backend PerplexDEW --dh-model Davies.
Requires a local reaktoro4py build with PerplexDEW symbols in one of:
  build/Reaktoro/Release
  build/Reaktoro/Release
  build/Reaktoro/Release

Water EOS: Zhang-Duan 2005.

Outputs (in Mineral_Solubilities/):
  calcite_solubility_comparison_high_P_dew24_PerplexDEW_Davies.png
  calcite_solubility_residuals_dew24_PerplexDEW_Davies.png
  calcite_speciation_dew24_PerplexDEW_Davies.png
"""

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

sys.argv = [sys.argv[0], "--backend", "PerplexDEW", "--dh-model", "Davies"]

import calcite_solubility_analysis as _cs  # noqa: E402

_cs.main()

