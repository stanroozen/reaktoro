"""
Quartz solubility analysis — DEW backend, Zhang-Duan 2005 EOS.

Runs quartz_solubility_analysis_v2_dew24.py with --backend DEW.
Outputs:
  quartz_solubility_comparison_low_P_dew24_DEW.png
  quartz_solubility_comparison_high_P_dew24_DEW.png
  quartz_solubility_residuals_dew24_DEW.png
  quartz_residuals_dew24_DEW.csv
  quartz_curves_dew24_DEW.csv
"""

import sys
import os

# Ensure the benchmark directory is on the path so quartz_solubility_analysis_v2_dew24
# can be imported as a module.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

# Patch sys.argv so argparse in the shared script sees the desired backend.
sys.argv = [sys.argv[0], "--backend", "DEW"]

import quartz_solubility_analysis_v2_dew24 as _qs  # noqa: E402

_qs.main()
