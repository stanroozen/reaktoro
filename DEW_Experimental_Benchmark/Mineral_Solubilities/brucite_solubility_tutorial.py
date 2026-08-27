"""
Brucite solubility tutorial entry point.

Delegates to the full uncertainty-enabled implementation in:
DEW_Experimental_Benchmark/Mineral_Solubilities/brucite/
    brucite_solubility_analysis_v2_dew24_uncertainty.py
"""

import os
import runpy

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TARGET = os.path.join(
    SCRIPT_DIR,
    "brucite",
    "brucite_solubility_analysis_v2_dew24_uncertainty.py",
)

if __name__ == "__main__":
    runpy.run_path(TARGET, run_name="__main__")
