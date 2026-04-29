import sys, autodiff
sys.path.insert(0, r'C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\build\Reaktoro\Release')
import reaktoro4py as r
globals().update({n: getattr(r,n) for n in dir(r) if not n.startswith('_')})
# Warnings.disable(906)
dew_db = DEWDatabase('dew2024-aqueous')
supcrt_db = SupcrtDatabase('supcrtbl')
combined_db = Database(dew_db.species())
combined_db.addSpecies(supcrt_db.species('Quartz'))
aq = AqueousPhase('WATER,AQ H+ OH- SiO2_aq H2_aq O2_aq HO2- HSiO3- Si2O4_aq Si3O6_aq')
aq.setActivityModel(ActivityModelPerplexDEW())
sys_ = ChemicalSystem(combined_db, aq, MineralPhase('Quartz'))
specs = EquilibriumSpecs(sys_)
specs.temperature()
specs.pressure()
solver = EquilibriumSolver(specs)
conds = EquilibriumConditions(specs)

conditions = [(200, 500), (250, 1000), (300, 1000), (400, 2000), (500, 5000)]
for T_C, P_bar in conditions:
    state = ChemicalState(sys_)
    state.set('WATER,AQ', autodiff.real(1.0), 'kg')
    state.set('H+', autodiff.real(1e-8), 'mol')
    state.set('OH-', autodiff.real(1e-8), 'mol')
    state.set('SiO2_aq', autodiff.real(1e-6), 'mol')
    state.set('Quartz', autodiff.real(10.0), 'mol')
    conds.temperature(float(T_C), 'celsius')
    conds.pressure(float(P_bar), 'bar')
    result = solver.solve(state, conds)
    if result.succeeded():
        ap = AqueousProps(state)
        si_m = float(ap.elementMolality('Si'))
        print(f'T={T_C}C P={P_bar}bar: Si={si_m:.4e}')
    else:
        print(f'T={T_C}C P={P_bar}bar: FAILED')

print('Done baseline test.')


