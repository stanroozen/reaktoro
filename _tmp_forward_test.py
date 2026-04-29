import json
import autodiff
from pathlib import Path
from reaktoro import *

def to_real(v):
    try:
        return autodiff.real(v)
    except Exception:
        return v

d = json.loads(Path('embedded/databases/hollandpowell/tc-ds62-reaktoro.json').read_text())
qparams = d['Species']['q']['StandardThermoModel']['HollandPowell']

p = StandardThermoModelParamsHollandPowell()
for k,v in qparams.items():
    if hasattr(p,k):
        setattr(p, k, to_real(float(v)))
p.Gf = to_real(float(qparams['Gf']) + 1000.0)
model = StandardThermoModelHollandPowell(p)

sup = SupcrtDatabase('supcrtbl')
dew = DEWDatabase('dew2024-aqueous')
q = sup.species('Quartz').withStandardThermoModel(model)

db = Database(dew.species())
db.addSpecies(q)

# Minimal system build and one solve
aq = AqueousPhase('WATER,AQ H+ OH- SiO2_aq H2_aq O2_aq HO2- HSiO3- Si2O4_aq Si3O6_aq')
aq.setActivityModel(ActivityModelPerplexDEW(ActivityDHModel.Davies))
system = ChemicalSystem(db, aq, MineralPhase('Quartz'))

specs = EquilibriumSpecs(system)
specs.temperature(); specs.pressure()
solver = EquilibriumSolver(specs)
conds = EquilibriumConditions(specs)
state = ChemicalState(system)
state.set('WATER,AQ', to_real(1.0), 'kg')
state.set('H+', to_real(1e-8), 'mol')
state.set('OH-', to_real(1e-8), 'mol')
state.set('SiO2_aq', to_real(1e-6), 'mol')
state.set('Quartz', to_real(10.0), 'mol')
conds.temperature(300.0, 'celsius')
conds.pressure(1000.0, 'bar')
res = solver.solve(state, conds)
print('ok', res.succeeded())
if res.succeeded():
    print(float(AqueousProps(state).elementMolality('Si')))
