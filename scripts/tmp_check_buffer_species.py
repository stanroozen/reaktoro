import ast
import json

buf_path = r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\DEW_Experimental_Benchmark\Mineral_Solubilities\buffer_fO2_from_supcrtbl.py"
db_path = r"c:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\embedded\databases\reaktoro\supcrtbl.json"

mod = ast.parse(open(buf_path, "r", encoding="utf-8").read(), filename=buf_path)
buffer_reactions = None
for node in mod.body:
    if isinstance(node, ast.Assign):
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id == "BUFFER_REACTIONS":
                buffer_reactions = ast.literal_eval(node.value)
                break
    if buffer_reactions is not None:
        break

if buffer_reactions is None:
    raise SystemExit("BUFFER_REACTIONS not found")

species = set()
for category, mapping in buffer_reactions.items():
    for reaction in mapping.values():
        for name, _ in reaction:
            species.add(name)

with open(db_path, "r", encoding="utf-8") as f:
    data = json.load(f)

available = set(data.get("Species", {}).keys())
missing = sorted([s for s in species if s not in available])

print("Total species in buffers:", len(species))
print("Missing species:", len(missing))
for m in missing:
    print("  -", m)
