# DEW / PerplexDEW Integration Backlog

Prioritized implementation backlog — ordered by impact-to-effort ratio, each scoped as a single concrete slice of work.

---

## P1 — Highest impact, lowest risk (unblocks everything else)

| # | Task | Target file(s) | Size | Status |
|---|------|----------------|------|--------|
| 1 | Unify `ActivityModelParamsDEW` and `ActivityModelParamsPerplexDEW` shape — add `ActivityDHModel` enum to DEW params so the DH variant is switchable on both backends with the same field name | `ActivityModelDEW.hpp/.cpp/.py.cxx` | S | ✅ Done |
| 2 | Promote `densityTolerance` and `pSatRelTol` to Python on both models — last internal-only knobs blocking full Python tuning | `WaterModels.py.cxx`, `ActivityModelDEW.py.cxx` | S | ✅ Done |
| 3 | Add `ActivityModelPerplexDEW` to the DEW brucite tutorial as a drop-in alternative (one line swap) with a comment explaining when to prefer each — proves swappability at the tutorial level | `DEW_Experimental_Benchmark/Tutorial/brucite_solubility_tutorial.py` | XS | ✅ Done |

---

## P2 — Kinetics/transport bridge

| # | Task | Target file(s) | Size | Status |
|---|------|----------------|------|--------|
| 4 | Write `test_kinetics_dew.py` — constructs a `KineticsSolver` run with both `ActivityModelDEW` and `ActivityModelPerplexDEW`, asserts both converge and produce matching elemental molalities | `Testing/regression/smoke/` | M | ✅ Done |
| 5 | Write `test_transport_dew.py` — 1D reactive transport step with DEW aqueous phase, asserts `ChemicalProps` fields are populated without `nan` | `Testing/regression/smoke/` | M | ✅ Done |
| 6 | Audit `KineticsSolver.cpp` and `TransportSolver.cpp` for any hard-coded checks on activity model type that would reject non-standard models — file issue or fix inline if trivial | `Reaktoro/Kinetics/KineticsSolver.cpp`, `Reaktoro/Transport/TransportSolver.cpp` | S | ✅ Done (no type checks found — both files are model-agnostic) |

---

## P3 — Constraint/inverse/sensitivity hooks

| # | Task | Target file(s) | Size | Status |
|---|------|----------------|------|--------|
| 7 | Smoke-test `specs.pH()` with `ActivityModelPerplexDEW` — port `test_dew_ph_constraint.py` (now in `Testing/scripts/`) to a proper pytest, assert pH constraint converges at a fixed T/P | `Testing/regression/smoke/` | S | ✅ Done |
| 8 | Smoke-test `EquilibriumSensitivity` with both models — compute `dAmount_dT` for one species, assert it's finite and non-zero | `Testing/regression/smoke/` | S | ✅ Done |
| 9 | Add `warningif` in `EquilibriumSpecs` when `specs.pH()` or `specs.fugacity()` is called while the active aqueous model does not expose `AqueousMixtureState` in `props.extra` — currently silent fallback | `Reaktoro/Equilibrium/EquilibriumSpecs.cpp` | S | ✅ Done |

---

## P4 — Multiphase / GFSM coupling standardization

| # | Task | Target file(s) | Size | Status |
|---|------|----------------|------|--------|
| 10 | Add `warnOnUnmappedGFSMCoupling` to the PerplexDEW tutorial with options comment + physical meaning (parallel to the existing `errorOnConflictingStandardState` block) | `DEW_Experimental_Benchmark/Tutorial/brucite_solubility_tutorial_dew24hp622.py` | XS | ✅ Done |
| 11 | Expose `ActivityModelPerplexGFSM` params to Python (currently C++-only) so the gas-phase GFSM handoff state-ID mechanism can be tested and tuned from scripts | `Reaktoro/Extensions/Perple_X/ActivityModelPerplexGFSM.py.cxx` (create) | M | ✅ Done (expanded existing stub with `HybridEosOptions`, `MrkMixOptions`, `PerpleXPureEosOptions` bindings + nested enums) |
| 12 | Write `test_gfsm_handoff.py` — two-phase system (aqueous + GFSM gas), assert `PerplexGFSM::WaterActivity::StateId` handoff is consumed by the aqueous phase and water activity differs from the pure-water value | `Testing/regression/smoke/` | M | ✅ Done |

---

## P5 — API harmonization and docs

| # | Task | Target file(s) | Size | Status |
|---|------|----------------|------|--------|
| 13 | Update `DEW_PERPLEXDEW_SETTABLE_API.md` with `warnOnUnmappedGFSMCoupling` entry, plus a "Model swapability cheatsheet" section comparing DEW vs PerplexDEW constructor signatures | `DEW_Experimental_Benchmark/DEW_PERPLEXDEW_SETTABLE_API.md` | S | ✅ Done |
| 14 | Add a capability matrix comment block at the top of both `ActivityModelDEW.hpp` and `ActivityModelPerplexDEW.hpp` listing supported/unsupported workflow categories — makes it self-documenting | Both `.hpp` files | S | ✅ Done |
| 15 | Phase 5 mandatory test matrix — a single pytest conftest-level marker (`@pytest.mark.workflow_coverage`) on the most representative test per workflow category, so `pytest -m workflow_coverage` produces the acceptance checklist | `Testing/conftest.py` | S | ✅ Done |

---

## Suggested Execution Order

```
1 → 3 → 2   (Python surfaces harmonized, swap proven at tutorial level)
→ 7 → 8 → 4 → 5   (constraint + kinetics smoke tests catch regressions early)
→ 9 → 6   (diagnostics close silent-failure gaps)
→ 10 → 11 → 12   (GFSM coupling standardized last, highest complexity)
→ 13 → 14 → 15   (docs + coverage matrix as final gate)
```

---

## Notes

- **`ActivityDHModel` enum** is declared in `ActivityModelDEW.hpp` and shared by `ActivityModelParamsPerplexDEW` via include — no duplicate registration in Python bindings.
- **`TransportSolver`** has no Python binding; Task 5 (`test_transport_dew.py`) implements operator splitting directly (equilibrate → advect loop) at the Python level.
- **Size legend**: XS < 30 min · S < 2 h · M < 1 day · L > 1 day
