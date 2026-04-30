import os
import runpy
import sys
import traceback
import types

for directory in [
    r"C:/Users/stanroozen/Documents/Projects/reaktoro-dev/reaktoro/build-msvc/Reaktoro/Release",
    r"C:/Users/stanroozen/anaconda3",
    r"C:/Users/stanroozen/anaconda3/DLLs",
    r"C:/Users/stanroozen/anaconda3/Library/bin",
    r"C:/Users/stanroozen/anaconda3/envs/reaktoro/Library/bin",
]:
    os.add_dll_directory(directory)

sys.path.insert(
    0,
    r"C:/Users/stanroozen/Documents/Projects/reaktoro-dev/reaktoro/build-msvc/Reaktoro/Release",
)

import reaktoro4py

reaktoro_module = types.ModuleType("reaktoro")
for name in dir(reaktoro4py):
    if name.startswith("_"):
        continue
    setattr(reaktoro_module, name, getattr(reaktoro4py, name))
sys.modules["reaktoro"] = reaktoro_module

sys.argv = [
    r"DEW_Experimental_Benchmark/Mineral_Solubilities/perplex_mixed_fluid_parity.py",
    "--disable-hydrous-species-correction",
    "--disable-water-activity-correction",
    "--out-csv",
    r"DEW_Experimental_Benchmark/Mineral_Solubilities/perplex_mixed_fluid_parity_results_no_hydrous_water_activity.csv",
]

try:
    runpy.run_path(
        r"DEW_Experimental_Benchmark/Mineral_Solubilities/perplex_mixed_fluid_parity.py",
        run_name="__main__",
    )
except BaseException:
    traceback.print_exc()
    raise
