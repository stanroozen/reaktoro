"""
Quartz solubility analysis â€” PerplexDEW backend, Extended Debye-HÃ¼ckel activity model.

Runs quartz_solubility_analysis_v2_dew24.py with --backend PerplexDEW --dh-model ExtendedDH.
Requires a local reaktoro4py build with PerplexDEW symbols in one of:
  build/Reaktoro/Release
  build/Reaktoro/Release
  build/Reaktoro/Release

Water EOS: Zhang-Duan 2005 (same as DEW backend default).

Outputs:
  quartz_solubility_comparison_low_P_dew24_PerplexDEW_ExtendedDH.png
  quartz_solubility_comparison_high_P_dew24_PerplexDEW_ExtendedDH.png
  quartz_solubility_residuals_dew24_PerplexDEW_ExtendedDH.png
  quartz_residuals_dew24_PerplexDEW_ExtendedDH.csv
  quartz_curves_dew24_PerplexDEW_ExtendedDH.csv
"""

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

sys.argv = [sys.argv[0], "--backend", "PerplexDEW", "--dh-model", "ExtendedDH"]

import quartz_solubility_analysis_v2_dew24 as _qs  # noqa: E402

_qs.main()

