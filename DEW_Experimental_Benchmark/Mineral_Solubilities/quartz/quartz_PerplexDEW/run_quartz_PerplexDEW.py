"""
Quartz solubility analysis â€” PerplexDEW backend, Davies activity model.

Runs quartz_solubility_analysis_v2_dew24.py with --backend PerplexDEW --dh-model Davies.
Requires a local reaktoro4py build with PerplexDEW symbols in one of:
  build/Reaktoro/Release
  build/Reaktoro/Release
  build/Reaktoro/Release

Water EOS: Zhang-Duan 2005 (same as DEW backend default).

Outputs:
  quartz_solubility_comparison_low_P_dew24_PerplexDEW_Davies.png
  quartz_solubility_comparison_high_P_dew24_PerplexDEW_Davies.png
  quartz_solubility_residuals_dew24_PerplexDEW_Davies.png
  quartz_residuals_dew24_PerplexDEW_Davies.csv
  quartz_curves_dew24_PerplexDEW_Davies.csv
"""

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

sys.argv = [sys.argv[0], "--backend", "PerplexDEW", "--dh-model", "Davies"]

import quartz_solubility_analysis_v2_dew24 as _qs  # noqa: E402

_qs.main()

