import importlib.util

p = r"DEW_Experimental_Benchmark/Tutorial/willemite_solubility_tutorial_dew17hp622_zn.py"
s = importlib.util.spec_from_file_location("m", p)
m = importlib.util.module_from_spec(s)
s.loader.exec_module(m)

print("module loaded")
system = m.build_tutorial_system(include_gas=False)
print("system built")
specs = m.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
print("specs built")
solver = m.EquilibriumSolver(specs)
print("solver built", bool(solver is not None))
