import os, sys, importlib
root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
cands = [
    os.path.join(root, 'build', 'Reaktoro', 'Release'),
    os.path.join(root, 'build', 'Reaktoro', 'Release'),
    os.path.join(root, 'build', 'Reaktoro', 'Release'),
]
for p in cands:
    if os.path.isdir(p) and p not in sys.path:
        sys.path.insert(0, p)
mod = importlib.import_module('reaktoro4py')
print('loaded:', getattr(mod, '__file__', 'n/a'))
names = ['DEWDatabase','DEW2024','SupcrtDatabase','SUPCRTBL','ActivityModelPerplexDEW','ActivityModelPerplexGFSM','ActivityDHModel','PerpleXWaterEos','PerpleXCO2Eos']
for n in names:
    print(f'{n}={hasattr(mod,n)}')

