# DEW and PerplexDEW Settable API

This file lists the current settable API surface for the DEW and PerplexDEW models in this repository.

Scope notes:

- "Settable API" means values a caller can configure through public C++ headers or through the current Python bindings.
- This file is based on the code currently in `Reaktoro/Extensions/DEW`, `Reaktoro/Extensions/Perple_X`, and the related pybind11 exports.
- When a field exists in C++ but is not exposed in Python, that is called out explicitly.
- `PerplexDEW` below uses the spelling in code. Some older notes in the repo use `Perplex` or `Perple_X` interchangeably.

## Python-Exposed DEW API

### 1. `DEWDatabase`

Purpose: load DEW YAML databases and access embedded DEW datasets.

Python entrypoints:

- `DEWDatabase()`
- `DEWDatabase(name: str)`
- `load(database: str)`
- `DEWDatabase.withName(name: str)`
- `DEWDatabase.fromFile(path: str)`
- `DEWDatabase.fromContents(contents: str)`
- `DEWDatabase.contents(name: str)`
- `DEWDatabase.namesEmbeddedDatabases()`

Typical embedded names mentioned in code:

- `dew2024-aqueous`
- `dew2019-aqueous`
- `dew2024-gas`
- `dew2019-gas`

Source:

- `Reaktoro/Extensions/DEW/DEWDatabase.hpp`
- `Reaktoro/Extensions/DEW/DEWDatabase.py.cxx`

### 2. `WaterModelOptions`

Purpose: configure which DEW water-property submodels are used inside `StandardThermoModelDEW`.

Python-assignable fields:

- `eosModel`
- `dielectricModel`
- `gibbsModel`
- `bornModel`
- `usePsatPolynomials`
- `psatRelTol`
- `densityTolerance`

Python helper:

- `makeWaterModelOptionsDEW()`

Python enums used by `WaterModelOptions`:

#### `WaterEosModel`

- `WaterEosModel.WagnerPruss`
- `WaterEosModel.HGK`
- `WaterEosModel.ZhangDuan2005`
- `WaterEosModel.ZhangDuan2009`

#### `WaterDielectricModel`

- `WaterDielectricModel.JohnsonNorton1991`
- `WaterDielectricModel.Franck1990`
- `WaterDielectricModel.Fernandez1997`
- `WaterDielectricModel.PowerFunction`

#### `WaterGibbsModel`

- `WaterGibbsModel.DelaneyHelgeson1978`
- `WaterGibbsModel.DewIntegral`

#### `WaterBornModel`

- `WaterBornModel.None`
- `WaterBornModel.Shock92Dew`

Source:

- `Reaktoro/Extensions/DEW/WaterModelOptions.hpp`
- `Reaktoro/Extensions/DEW/WaterModels.py.cxx`

### 3. `StandardThermoModelParamsDEW`

Purpose: configure the DEW-based HKF standard thermodynamic model for an aqueous species.

Python-assignable fields:

- `Gf`
- `Hf`
- `Sr`
- `a1`
- `a2`
- `a3`
- `a4`
- `c1`
- `c2`
- `wref`
- `charge`
- `Tmax`
- `waterOptions`

Python model factory:

- `StandardThermoModelDEW(params)`

Meaning of the fields:

- `Gf`: apparent standard Gibbs energy of formation at reference conditions.
- `Hf`: apparent standard enthalpy of formation at reference conditions.
- `Sr`: standard entropy at reference conditions.
- `a1`, `a2`, `a3`, `a4`, `c1`, `c2`: HKF coefficients.
- `wref`: Born coefficient at reference conditions.
- `charge`: species charge.
- `Tmax`: maximum valid temperature.
- `waterOptions`: nested `WaterModelOptions` controlling DEW water submodels.

Source:

- `Reaktoro/Models/StandardThermoModels/StandardThermoModelDEW.hpp`
- `Reaktoro/Models/StandardThermoModels/StandardThermoModelDEW.py.cxx`

### 4. `ActivityModelDEW`

Purpose: return the DEW aqueous activity-model generator.

Python entrypoints:

- `ActivityModelDEW()`
- `ActivityModelDEW(params)`

### 5. `ActivityModelParamsDEW`

Purpose: parameter object for the DEW activity model.

Python-assignable fields:

- `dhModel` — Debye-Hückel variant (`ActivityDHModel.ExtendedDH` by default; `ActivityDHModel.Davies` is also supported). The `ActivityDHModel` enum is declared here and shared with `ActivityModelParamsPerplexDEW`.
- `waterOptions`
- `bExtended`

Important:

- `errorOnConflictingStandardState` is not a DEW parameter. It is a PerplexDEW-only safety switch (see section 8).

Notes:

- `dhModel` selects the Debye-Hückel formula. Both variants use DEW-derived A and B parameters computed from the configured water submodels.
- `waterOptions` lets the DEW activity model use the same configurable DEW water-property surface already used by `StandardThermoModelDEW`.
- `bExtended` exposes the extended Debye-Hückel correction term that was previously fixed internally to `0.0`.

Current configurability status:

- The default `ActivityModelDEW()` still uses DEW-style defaults (ExtendedDH variant).
- The parameterized overload `ActivityModelDEW(params)` now allows callers to override:
  - Debye-Hückel variant through `params.dhModel`
  - water EOS and dielectric path through `params.waterOptions`
  - Psat handling and density tolerance through `params.waterOptions`
  - extended-term parameter `b_c,k` through `params.bExtended`
- Effective ionic radii remain fixed internally through the current lookup table with charge-based fallback.

So `ActivityModelDEW` is now both public API and a parameterized settable API in Python.

Source:

- `Reaktoro/Models/ActivityModels/ActivityModelDEW.hpp`
- `Reaktoro/Models/ActivityModels/ActivityModelDEW.cpp`
- `Reaktoro/Models/ActivityModels/ActivityModelDEW.py.cxx`

## Python-Exposed PerplexDEW API

### 6. `StandardThermoModelParamsPerplexDEW`

Purpose: configure the Perple_X-backed standard thermodynamic model that mirrors the DEW HKF-style parameter surface.

Python-assignable fields:

- `Gf`
- `Hf`
- `Sr`
- `a1`
- `a2`
- `a3`
- `a4`
- `c1`
- `c2`
- `wref`
- `charge`
- `Tmax`

Python model factory:

- `StandardThermoModelPerplexDEW(params)`

Notes:

- This surface intentionally matches `StandardThermoModelParamsDEW` closely, but there is no nested `waterOptions` field here.
- Internally the implementation converts `a1` and `a3` from Pa-based units into Perple_X bar-based units.

Source:

- `Reaktoro/Extensions/Perple_X/StandardThermoModelPerplexDEW.hpp`
- `Reaktoro/Models/StandardThermoModels/StandardThermoModelPerplexDEW.py.cxx`

### 7. `ActivityDHModel` (shared enum)

Purpose: select the Debye-Hückel variant.  Used by **both** `ActivityModelParamsDEW` and `ActivityModelParamsPerplexDEW` via the same field name `dhModel`.

Python enum values:

- `ActivityDHModel.ExtendedDH` — extended Debye-Hückel with ion-size parameter `a`; good to ~1 mol/kg.
- `ActivityDHModel.Davies` — Davies equation; no ionic-size parameter; good to ~0.5 mol/kg.

Default:

- `ActivityModelDEW`: `ExtendedDH`
- `ActivityModelPerplexDEW`: `Davies` (matches Perple_X GFSM `aqact` convention exactly)

Source:

- `Reaktoro/Models/ActivityModels/ActivityModelDEW.hpp` (declaration)
- `Reaktoro/Models/ActivityModels/ActivityModelDEW.py.cxx` (Python registration — registered once, shared)

### 8. `ActivityModelParamsPerplexDEW`

Purpose: parameter object for the PerplexDEW activity model.

Python-assignable fields:

- `dhModel` — Debye-Hückel variant (shared `ActivityDHModel` enum; default: `Davies`).
- `errorOnConflictingStandardState` — raise an error instead of a warning when a species is configured as both a GFSM solvent (mole-fraction reference state) and a PerplexDEW HKF solute (molal reference state). Default is `False`, but for coupled GFSM+PerplexDEW workflows use `True`.
- `warnOnUnmappedGFSMCoupling` — emit a one-time warning when a GFSM-phase species cannot be matched to any Reaktoro aqueous species during PerplexDEW evaluation; default `True`.
- `requireCoupledGFSMHandoff` — fail-fast switch for coupled-fluid workflows; when `True`, throw if PerplexDEW does not consume a fresh GFSM handoff for the current `ChemicalProps` state instead of silently falling back to aqueous-only water activity; default `False`.

Recommended strict setting for coupled-fluid runs:

```python
params = ActivityModelParamsPerplexDEW()
params.errorOnConflictingStandardState = True
```

Conflict validation behavior:

- `errorOnConflictingStandardState=False` (default): emit a warning if a species is configured as both a GFSM solvent and a PerplexDEW HKF solute.
- `errorOnConflictingStandardState=True`: throw an error for that conflict to prevent accidental double-counting.

Unmapped coupling behavior:

- `warnOnUnmappedGFSMCoupling=True` (default): print a one-time stderr warning for each GFSM species that cannot be matched. Benign when those species are absent from the aqueous phase.
- `warnOnUnmappedGFSMCoupling=False`: suppress all unmapped-coupling warnings.

Strict coupling behavior:

- `requireCoupledGFSMHandoff=False` (default): allow fallback to aqueous-neutral water activity when no valid GFSM handoff is available.
- `requireCoupledGFSMHandoff=True`: throw a runtime error unless `PerplexGFSM::WaterActivity::StateId` matches the current `Reaktoro::ChemicalProps::StateId` and `PerplexGFSM::WaterActivity::ln_f_ratio_h2o` is present.

Python factory overloads:

- `ActivityModelPerplexDEW(params)`
- `ActivityModelPerplexDEW(model=ActivityDHModel.Davies)`

Source:

- `Reaktoro/Extensions/Perple_X/ActivityModelPerplexDEW.hpp`
- `Reaktoro/Models/ActivityModels/ActivityModelPerplexDEW.py.cxx`

## C++ Public Options Not Currently Exposed to Python

These are still part of the settable public API in C++, but they are not currently bound in Python.

### DEW low-level options

#### `WaterSolventFunctionOptions`

Fields:

- `Psat`
- `densityEquation`

Purpose:

- controls the low-level DEW solvent-function evaluation used by omega calculations.

Source:

- `Reaktoro/Extensions/DEW/WaterSolventFunctionDEW.hpp`

#### `WaterBornOmegaOptions`

Fields:

- `solvent`
- `isHydrogenLike`
- `maxPressureForVariation`

Purpose:

- controls low-level Born omega evaluation.

Source:

- `Reaktoro/Extensions/DEW/WaterBornOmegaDEW.hpp`

#### `WaterStateOptions`

Fields:

- `thermo`
- `dielectric`
- `computeGibbs`
- `gibbs`
- `computeSolventG`
- `solvent`
- `computeOmega`
- `omega`

Purpose:

- controls the full DEW water-state orchestration layer.

Source:

- `Reaktoro/Extensions/DEW/WaterState.hpp`

Note:

- `StandardThermoModelDEW` maps `WaterModelOptions` into this deeper option graph internally.

## Practical Summary

If you are scripting in Python today, the settable surfaces you can actually use directly are:

### DEW options and APIs

- `DEWDatabase(...)`
- `DEWDatabase.load(...)`
- `DEWDatabase.withName(...)`
- `DEWDatabase.fromFile(...)`
- `DEWDatabase.fromContents(...)`
- `makeWaterModelOptionsDEW()`
- `WaterModelOptions.eosModel`
- `WaterModelOptions.dielectricModel`
- `WaterModelOptions.gibbsModel`
- `WaterModelOptions.bornModel`
- `WaterModelOptions.usePsatPolynomials`
- `WaterModelOptions.psatRelTol`
- `WaterModelOptions.densityTolerance`
- `StandardThermoModelParamsDEW.Gf`
- `StandardThermoModelParamsDEW.Hf`
- `StandardThermoModelParamsDEW.Sr`
- `StandardThermoModelParamsDEW.a1`
- `StandardThermoModelParamsDEW.a2`
- `StandardThermoModelParamsDEW.a3`
- `StandardThermoModelParamsDEW.a4`
- `StandardThermoModelParamsDEW.c1`
- `StandardThermoModelParamsDEW.c2`
- `StandardThermoModelParamsDEW.wref`
- `StandardThermoModelParamsDEW.charge`
- `StandardThermoModelParamsDEW.Tmax`
- `StandardThermoModelParamsDEW.waterOptions`
- `StandardThermoModelDEW(params)`
- `ActivityModelDEW()`
- `ActivityModelParamsDEW.dhModel`
- `ActivityModelParamsDEW.waterOptions`
- `ActivityModelParamsDEW.bExtended`
- `ActivityModelDEW(params)`

### PerplexDEW options and APIs

- `StandardThermoModelParamsPerplexDEW.Gf`
- `StandardThermoModelParamsPerplexDEW.Hf`
- `StandardThermoModelParamsPerplexDEW.Sr`
- `StandardThermoModelParamsPerplexDEW.a1`
- `StandardThermoModelParamsPerplexDEW.a2`
- `StandardThermoModelParamsPerplexDEW.a3`
- `StandardThermoModelParamsPerplexDEW.a4`
- `StandardThermoModelParamsPerplexDEW.c1`
- `StandardThermoModelParamsPerplexDEW.c2`
- `StandardThermoModelParamsPerplexDEW.wref`
- `StandardThermoModelParamsPerplexDEW.charge`
- `StandardThermoModelParamsPerplexDEW.Tmax`
- `StandardThermoModelPerplexDEW(params)`
- `ActivityDHModel.Davies`  ← shared with DEW; registered by `exportActivityModelDEW`
- `ActivityDHModel.ExtendedDH`
- `ActivityModelParamsPerplexDEW.dhModel`
- `ActivityModelParamsPerplexDEW.errorOnConflictingStandardState`
- `ActivityModelParamsPerplexDEW.warnOnUnmappedGFSMCoupling`
- `ActivityModelParamsPerplexDEW.requireCoupledGFSMHandoff`
- `ActivityModelPerplexDEW(params)`
- `ActivityModelPerplexDEW(model=...)`

---

## Model Swapability Cheatsheet

Both `ActivityModelDEW` and `ActivityModelPerplexDEW` implement the same `ActivityModelGenerator` interface and are plug-compatible at the `setActivityModel` call site.  The only difference is the parameter type.

| Aspect | `ActivityModelDEW` | `ActivityModelPerplexDEW` |
|--------|-------------------|---------------------------|
| **Parameter struct** | `ActivityModelParamsDEW` | `ActivityModelParamsPerplexDEW` |
| **DH variant field** | `dhModel` (default: `ExtendedDH`) | `dhModel` (default: `Davies`) |
| **Water submodels** | Configurable via `waterOptions` (EOS, dielectric, Gibbs, Born, Psat) | Perple_X ZD05 + Looyenga mixing; not configurable via `waterOptions` |
| **GFSM gas coupling** | Not supported | Supported via `props.extra` StateId handoff |
| **Conflict guard** | — | `errorOnConflictingStandardState` |
| **Unmapped GFSM warn** | — | `warnOnUnmappedGFSMCoupling` |
| **Strict GFSM handoff** | — | `requireCoupledGFSMHandoff` |
| **Best for** | DEW database directly, single-fluid aqueous chemistry | Perple_X-derived database, mixed fluid phases, matching Perple_X `aqact` convention |

### Minimal swap pattern (Python)

```python
# --- Using ActivityModelDEW ---
params = ActivityModelParamsDEW()
params.dhModel = ActivityDHModel.ExtendedDH   # or Davies
aqueous_phase.setActivityModel(ActivityModelDEW(params))

# --- Swap to ActivityModelPerplexDEW (one struct change, rest identical) ---
params = ActivityModelParamsPerplexDEW()
params.dhModel = ActivityDHModel.Davies        # default for PerplexDEW
params.warnOnUnmappedGFSMCoupling = True
params.requireCoupledGFSMHandoff = False       # set True for fail-fast coupled runs
aqueous_phase.setActivityModel(ActivityModelPerplexDEW(params))
```

The `ChemicalSystem`, solver, and results loop are unchanged for either model.

---

## Gaps Worth Noting

The current asymmetries are:

### DEW-specific notes

- DEW standard thermo is configurable in Python through `waterOptions`.
- DEW activity model is configurable through `ActivityModelParamsDEW.waterOptions` and `ActivityModelParamsDEW.bExtended`.
- `bExtended` has no PerplexDEW equivalent (PerplexDEW uses the Davies formula's 0.3I term directly).

### PerplexDEW-specific notes

- PerplexDEW standard thermo remains HKF-parameter based and does not currently expose a nested `waterOptions` object.
- PerplexDEW supports GFSM gas-phase coupling via `props.extra`; DEW does not.
- `errorOnConflictingStandardState`, `warnOnUnmappedGFSMCoupling`, and
  `requireCoupledGFSMHandoff` are PerplexDEW-only diagnostics/guards.

### Shared

- `ActivityDHModel` enum (`Davies`, `ExtendedDH`) is now declared in `ActivityModelDEW.hpp` and shared by both param structs via the unified `dhModel` field name.