# Hybrid Equilibrium Engine Design

## Strategic Summary

The best integration path is not "put MAGEMin's local minimizer inside one Reaktoro activity model and hope for the best." It is a **hybrid equilibrium engine** inside Reaktoro that:

- Keeps Reaktoro's constraint framework, outer optimization, and fluid/aqueous machinery unchanged.
- Adds MAGEMin-style candidate generation and branch screening for non-ideal solid solutions.
- Runs gradient-based local refinement per branch **after** screening, not instead of it.

This preserves fixed-pH, fugacity, open-component, enthalpy, and pressure constraints that are already correct in Reaktoro. DEW and PerplexDEW stay entirely untouched at the outer level.

---

## Architectural Layers

### Layer 1 — Reaktoro outer equilibrium framework (unchanged)

**Owner:** `EquilibriumSolver.cpp`, `EquilibriumSpecs.hpp`, `EquilibriumConditions.hpp`

Responsibilities:
- Single global optimization problem over species amounts and control variables.
- pH, fugacity, open components, fixed pressure, fixed temperature, enthalpy targets.
- Phase inclusion/exclusion from the outer Gibbs minimization.
- Constraint assembly and numerical refinement through Optima.

**Invariant:** This layer must never learn about MAGEMin internal coordinates. It only sees
`ActivityProps` (`Gx`, `Hx`, `ln_g`, `ln_a`) from whatever phase thermodynamics model is attached.

### Layer 2 — Globalized solid-solution seam (current; extend carefully)

**Owner:** `ActivityModelGlobalizedSolidSolution.hpp/.cpp`

Responsibilities:
- `GlobalizedSolidSolutionModel` function type: the sole contract between outer Layer 1 and inner
  Layer 3. Receives `GlobalizedSolidSolutionInput` (T, P, visible x, cached state, requested branch)
  and returns `GlobalizedSolidSolutionOutput` (Gx, Hx, ln_g, ln_a, branch metadata, split request).
- Branch metadata: `GlobalizedSolidSolutionBranch` with coordinate bounds per admissible branch.
- Split-request publication and consumption: `GlobalizedSolidSolutionSplitRequest` triggers outer
  phase duplication when the solid solution is immiscible.
- Phase duplication helpers: `AssembleGlobalizedSolidSolutionCandidatePhases(...)` creates one
  Reaktoro `Phase` per candidate state before the outer solve.
- Retry/rebuild loop: `equilibrateWithGlobalizedSolidSolutionSplits(...)` rebuilds the outer system
  and reequilibrates when a split was requested.

**Invariant:** The seam itself does not solve equilibrium. It routes thermodynamic callbacks and
carries warm-start state. Changing the outer solver does not require changing this layer.

### Layer 3 — MAGEMin-style candidate generation and local minimization (extend here)

**Owner:** `ActivityModelMAGEMinSolidSolutionPilot.hpp/.cpp`,
           `InternallyMinimizedSolidSolution.hpp`

Responsibilities:
- **Candidate generation:** Produce a small set of branch-local initial guesses before refinement.
  Mirrors MAGEMin `simplex_levelling.c` levelling and branch-proposal logic.
- **Branch-local constrained minimization:** Given a starting guess and a branch, minimize the
  internal Gibbs objective (excess + ideal terms + external-composition penalty) subject to
  simplex constraints. Mirrors MAGEMin `NLopt_opt_function.c` / `SB_NLopt_opt_function.c`.
- **Candidate screening and ranking:** Select the best branch-local result before returning to
  the outer solver. Discard thermodynamically inadmissible candidates early.
- **Telemetry:** Emit per-evaluation diagnostics (`SolidSolutionMinimizerBenchmark`,
  `SolidSolutionMinimizerTelemetry`) to guide rollout and flag families that need the
  projected-gradient path.

**Key types (all in `InternallyMinimizedSolidSolution.hpp`):**
- `InternalSolidSolutionModel` — model function contract for one internally minimized phase.
- `SolidSolutionCandidateState` / `SolidSolutionCandidateGenerator` — candidate-generation contracts.
- `SolidSolutionMinimizerTelemetry` / `SolidSolutionMinimizerBenchmark` — accumulator for fallback
  rate and iteration savings.

---

## Concrete Numerical Workflow

```
For each outer Reaktoro Newton iteration:
  1. Solve aqueous/fluid activity models normally (DEW, Pitzer, HKF, …).
  2. For each non-ideal solid solution:
     a. Generate a small set of branch-local candidate seeds (Layer 3 generator).
     b. For each candidate seed: run internal constrained minimization over site fractions.
     c. Rank candidates by branch-local Gibbs score; discard inadmissible ones.
     d. Return Gx, Hx, ln_g, ln_a from the best surviving candidate to Layer 1.
  3. Layer 1 assembles the full Jacobian and refines all species amounts together.
  4. If any solid solution emits a split request:
     a. Duplicate the phase into branch-specific instances.
     b. Rerun the outer solve with the expanded phase list.
  5. After convergence: recheck stability of each solid-solution branch.
     If new instabilities are found, go to step 4.
```

This is not formally globally optimal. It is MAGEMin's practical strategy: structured branch
generation + constrained local refinement + discard poor candidates early + final outer refinement.

---

## Implementation State

### Done

| Component | File | Status |
|-----------|------|--------|
| Layer 1: outer equilibrium framework | `Equilibrium/EquilibriumSolver.cpp`, `EquilibriumSpecs.hpp` | Unchanged, verified green |
| Layer 1: retry/rebuild loop | `EquilibriumUtils.cpp`, `equilibrateWithGlobalizedSolidSolutionSplits` | Working, tested |
| Layer 2: seam contracts | `ActivityModelGlobalizedSolidSolution.hpp` | Complete |
| Layer 2: candidate phase assembly | `ActivityModelGlobalizedSolidSolution.cpp`, `AssembleGlobalizedSolidSolutionCandidatePhases` | Implemented this session |
| Layer 3: reusable interface header | `InternallyMinimizedSolidSolution.hpp` | Created this session |
| Layer 3: telemetry accumulator | `InternallyMinimizedSolidSolution.hpp`, `SolidSolutionMinimizerBenchmark` | Implemented this session |
| Layer 3 pilots (binary) | `sb11_ol`, `sb11_wa`, `sb21_sp` | Green (legacy minimizer) |
| Layer 3 pilots (ternary — guarded PG rollout) | `sb11_cf`, `sb11_pv`, `sb11_ak`, `sb21_nal` | Green |
| PG convergence fix | `ActivityModelMAGEMinSolidSolutionPilot.cpp`, Armijo-failure → `converged=true` | Done this session |
| PG objective-winner telemetry | `InternallyMinimizedSolidSolution.hpp`, `pgLowerObjectiveCount` / `legacyLowerObjectiveCount` | Done this session |
| Benchmark grid test (Step 1) | `ActivityModelMAGEMinSolidSolutionPilotRegression.test.cxx`, `[Benchmark]` tag | Done; 5 sections green |
| Duplicated-phase outer-workflow coverage (Step 2) | `ActivityModelGlobalizedSolidSolution.test.cxx`, `[EquilibriumRetry][OuterWorkflow]` | Done this session |
| Layer 3 pilots (ternary — direct PG) | `sb21_cf` | Green |
| Regression harness | `ActivityModelMAGEMinSolidSolutionPilotRegression.test.cxx` | 9 test cases, all green |
| JSON fixture set | `Support/MAGEMinRegressionFixtures/` | sb11 binary + ternary + sb21_nal |
| First physically grounded immiscibility pilot (Step 2 physical) | `ActivityModelGlobalizedSolidSolution.test.cxx`, `[Exsolution]` tag | Done this session — 2 test cases: split trigger + solvus activity signatures |

### In Progress / Not Yet Started

| Next Step | Layer | Risk | Notes |
|-----------|-------|------|-------|
| Stability / tangent-plane check per branch | Layer 3 | Medium | Required for robust solvus enforcement in the outer equilibration |
| Immiscibility outer equilibration (post-split phase separation) | Layers 2+3 | High | Current equal-split initialization is a saddle point; needs seeded or tangent-plane-guided init |
| Promote PG to default without guard on green ternary families | Layer 3 | Low | Benchmark-guided; use `SolidSolutionMinimizerBenchmark` |
| Port remaining SB21 families (beyond sb21_cf, sb21_nal) | Layer 3 | Low | Mechanical, follow existing ternary pattern |
| Port multi-site families (> 3 endmembers) | Layer 3 | Medium | Requires N-simplex projected gradient |
| Database-backed plugin system | Layer 3 | Low | After ≥3 families proven stable |
| Benchmarks vs MAGEMin on buffered rock-fluid cases | All | Low | Final validation step |

---

## What Must Not Change

1. `EquilibriumSolver`, `EquilibriumSpecs`, `EquilibriumConditions` — Layer 1 stays as-is.
2. `ActivityModelDEW`, `ActivityModelPerplexDEW`, `ActivityModelPerplexGFSM` — fluid/aqueous
   models stay native; they coexist in the same outer problem.
3. `GlobalizedSolidSolutionModel` function signature — Layer 2 seam must remain stable so
   new families can be plugged in without changing the outer solver.
4. `equilibrateWithGlobalizedSolidSolutionSplits(...)` — the retry/rebuild loop signature
   and behavior should only be changed to add new features, never to break existing callers.

---

## Rollout Policies for Layer 3 Minimizers

Families advance through three stages:

| Stage | Policy | Guard |
|-------|--------|-------|
| **Legacy** | Bounded coordinate search only | — |
| **Guarded PG** | Projected-gradient default; compare vs legacy; fall back on disagreement | `compareProjectedGradientAgainstLegacy = true`, `fallbackToLegacyOnProjectedGradientDisagreement = true` |
| **Direct PG** | Projected-gradient unconditionally | No guard |

Promotion from Guarded PG → Direct PG requires: fallback rate < 1 % over ≥ 100 evaluations
from the `SolidSolutionMinimizerBenchmark`, verified on the `build-msvc` Release path.

Current family status:

| Family | Endmembers | Stage |
|--------|-----------|-------|
| `sb11_ol` | `fo`, `fa` | Legacy |
| `sb11_wa` | `mgwa`, `fewa` | Legacy |
| `sb21_sp` | `sp`, `hc` | Legacy |
| `sb11_cf` | `mgcf`, `fecf`, `nacf` | Guarded PG |
| `sb11_pv` | `mgpv`, `fepv`, `alpv` | Guarded PG |
| `sb11_ak` | `mgak`, `feak`, `co` | Guarded PG |
| `sb21_nal` | `mnal`, `fnal`, `nnal` | Guarded PG |
| `sb21_cf` | `MgAl2O4`, `FeAl2O4`, `NaAlSiO4` | Direct PG |

---

## Next Steps (Ordered by Priority)

### Step 1 — Run the telemetry benchmark to decide PG promotion

**Effort:** Low. **Risk:** None.

For every Guarded PG family, run `SolidSolutionMinimizerBenchmark` over a grid of representative
compositions at relevant T/P conditions. Read `fallbackRate` and `averageProjectedGradientIterations`.

Families with fallback rate < 1 % can be promoted to Direct PG in their builder by removing
`compareProjectedGradientAgainstLegacy` and `fallbackToLegacyOnProjectedGradientDisagreement`.

#### Benchmark results — T=1473 K, P=1 GPa, 36-point ternary simplex grid (N=7)

| Family | total | pg | legacy | fallbackRate | pgWinsObj | legacyWinsObj |
|--------|-------|----|--------|--------------|-----------|---------------|
| sb11_pv | 36 | 26 | 10 | 27.8 % | 0 | 10 |
| sb11_ak | 36 | 26 | 10 | 27.8 % | 0 | 10 |
| sb11_cf | 36 | 25 | 11 | 30.6 % | 0 | 11 |
| sb21_nal | 36 | 27 | 9 | 25.0 % | 0 | 9 |
| sb21_cf (Direct PG baseline) | 36 | 36 | 0 | 0 % | — | — |

**Key finding:** In every disagreement case across all four guarded families, legacy finds a
strictly lower Gibbs energy than projected-gradient (`legacyWinsObj == fallbackCount` for all
families, `pgWinsObj == 0` for all). PG is converging to higher-energy local minima on
~28 % of compositions. The guard is necessary and working correctly.

**Bug fixed during benchmarking:** The Armijo backtracking loop exited without setting
`converged = true`, so the agreement check always failed (`converged` is a precondition).
After the fix, PG correctly reports convergence when it can no longer improve (stuck at a
local minimum). This raised the pg-wins count from 0 to ~25-27 per family.

**Promotion decision: none.** All four guarded families remain Guarded PG. Fallback rates of
25-31 % far exceed the 1 % promotion threshold. To achieve promotion these families need
either a better PG starting point (multi-start or warm-started from legacy) or a proof that
the 28 % disagreement compositions are thermodynamically equivalent.

### Step 2 — Immiscibility/exsolution pilot (first duplicated-phase outer workflow)

**Effort:** Medium. **Risk:** Medium. **Target family:** `sb11_pv` or `sb21_cf`.

`AssembleGlobalizedSolidSolutionCandidatePhases(...)` already exists in the seam. The missing
piece is the outer-loop driver that:

1. Constructs the initial single-phase system.
2. Calls `equilibrateWithGlobalizedSolidSolutionSplits(...)`.
3. On split request, calls `AssembleGlobalizedSolidSolutionCandidatePhases(...)` to rebuild the
   phase list with duplicated branch instances.
4. Re-equlibrates with the expanded system and verifies stability.

This is the minimum needed to demonstrate true immiscibility in a Reaktoro equilibrium problem.
Suggested first test: two-branch perovskite at a composition known to split.

#### Step 2 coverage update

An explicit outer-workflow coverage test now exists in
`ActivityModelGlobalizedSolidSolution.test.cxx` under tag
`[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][EquilibriumRetry][OuterWorkflow]`.

The test executes the intended manual workflow:

1. Run a first single-phase outer solve.
2. Assemble duplicated candidate phases using `AssembleGlobalizedSolidSolutionCandidatePhases(...)`.
3. Re-equilibrate the expanded system and verify branch-suffixed phase identities.

This closes the Step 2 "outer-loop driver coverage" gap.

#### Step 2 physical pilot update

The first physically grounded immiscibility pilot now lives in
`ActivityModelGlobalizedSolidSolution.test.cxx` under tag `[Exsolution]` and consists of two
test cases for the `sb21_cf` (calcioferrite) family (W02 = W12 = 60825.08 J/mol,
T_crit ≈ 3657 K >> T = 1473.15 K):

1. **Split-trigger test** (`[EquilibriumRetry][Exsolution]`): verifies that a bulk composition
   y ≈ (0.30, 0.20, 0.50) lying deep inside the nacf–(mgcf+fecf) solvus triggers a split
   request and produces a two-phase system with branch-named phases.

2. **Solvus activity-signature test** (`[Exsolution]`): verifies that the sb21_cf activity
   model produces physically distinct nacf activities at the two solvus arm compositions
   (y_left ≈ (0.45, 0.45, 0.10) and y_right ≈ (0.05, 0.05, 0.90)), confirming the activity
   landscape that drives the immiscibility.

**Known limitation**: the post-split 2-phase outer equilibration currently converges to equal
compositions from an equal-split initial state (symmetric saddle point).  True solvus
enforcement in the outer equilibration requires Step 3 (tangent-plane stability criterion) or
a branch-seeded initialization strategy.

### Step 3 — Branch stability / tangent-plane check

**Effort:** Medium. **Risk:** Low.

Implement a lightweight tangent-plane distance criterion as a
`GlobalizedSolidSolutionCandidateStabilityCriterion`. The criterion receives (T, P, visible x,
branch-local x, branch-local Gibbs) and returns `stable=false` + a `splitRequest` when the
driving force for a second branch is positive.

This replaces the current heuristic split trigger with a thermodynamically grounded check
and will be required before the outer loop can reliably detect missed phase splits.

### Step 4 — Port remaining SB21 / higher-order families

**Effort:** Low per family. **Risk:** Low.

Follow the exact same pattern as the existing ternary pilots. Priority order based on presence
in MAGEMin SB21 database:
- `sb21_cpx` (clinopyroxene, multi-site, > 3 endmembers) — first test of N-simplex PG
- `sb21_opx` (orthopyroxene)
- `sb21_gt` (garnet)

These require the N-simplex projected-gradient path to be tested and validated before the
existing ternary path can be extended to larger systems.

### Step 5 — Mixed rock-fluid benchmark vs MAGEMin

**Effort:** Medium. **Risk:** Low.

After at least one immiscibility pilot works and two or more multi-endmember families are ported,
set up a benchmark case with:
- Fixed pH (EquilibriumSpecs with pH constraint)
- Open fluid component (H₂O or CO₂)
- Two or three non-ideal solid solutions
- T and P grid covering a transition zone

Compare phase assemblage, solid-solution compositions, and convergence iteration counts against
a direct MAGEMin run at the same conditions. This is the first real validation of the hybrid
engine as a drop-in practical alternative.

### Step 6 — Generalize to database-backed plugin

**Effort:** Medium. **Risk:** Low.

Once ≥ 3 families and the immiscibility pilot are stable, promote the current hard-coded pilot
pattern into a database-backed plugin interface:

```cpp
struct SolidSolutionDatabaseFamily
{
    String familyId;
    InternalSolidSolutionModel            model;
    SolidSolutionCandidateGenerator       candidateGenerator;
    Vec<GlobalizedSolidSolutionBranch>    defaultBranches;
    SolidSolutionMinimizerBenchmark       benchmark;
};
```

Different MAGEMin databases (SB11, SB21, Holland-Powell, etc.) register their families
through this interface. The outer seam and equilibrium solver remain unaware of which database
is in use.

---

## Key Files (Reference)

| File | Role |
|------|------|
| `Reaktoro/Equilibrium/EquilibriumSolver.cpp` | Layer 1: outer optimization |
| `Reaktoro/Equilibrium/EquilibriumSpecs.hpp` | Layer 1: constraint declaration |
| `Reaktoro/Equilibrium/EquilibriumUtils.cpp` | Layer 2: retry/rebuild loop |
| `Reaktoro/Models/ActivityModels/ActivityModelGlobalizedSolidSolution.hpp` | Layer 2: seam contracts |
| `Reaktoro/Models/ActivityModels/ActivityModelGlobalizedSolidSolution.cpp` | Layer 2: seam + phase assembly implementations |
| `Reaktoro/Models/ActivityModels/InternallyMinimizedSolidSolution.hpp` | Layer 3: reusable internal-minimization interface |
| `Reaktoro/Models/ActivityModels/ActivityModelMAGEMinSolidSolutionPilot.hpp` | Layer 3: MAGEMin pilot public API |
| `Reaktoro/Models/ActivityModels/ActivityModelMAGEMinSolidSolutionPilot.cpp` | Layer 3: all pilot implementations |
| `Reaktoro/Models/ActivityModels/ActivityModelMAGEMinSolidSolutionPilotRegression.test.cxx` | Layer 3: regression + telemetry tests |
| `Reaktoro/Models/ActivityModels/Support/MAGEMinRegressionFixtures/` | Layer 3: JSON baseline fixtures |
| `Reaktoro/Models/ActivityModels/ActivityModelDEW.hpp` | Fluid/aqueous: stays untouched |
| `Reaktoro/Models/ActivityModels/ActivityModelPerplexDEW.hpp` | Fluid/aqueous: stays untouched |

---

## MAGEMin Source Mapping

| MAGEMin source | Reaktoro equivalent |
|----------------|---------------------|
| `simplex_levelling.c` | `SolidSolutionCandidateGenerator` (Layer 3) |
| `MAGEMin.c` (staging, phase selection) | `equilibrateWithGlobalizedSolidSolutionSplits` + `AssembleGlobalizedSolidSolutionCandidatePhases` (Layer 2) |
| `NLopt_opt_function.c` (general phase minimizers) | `MAGEMinSolidSolutionPilotModelImportedConstrainedTernary` + projected-gradient path (Layer 3) |
| `SB_NLopt_opt_function.c` (SB database-specific) | Per-family builder functions: `MAGEMinSolidSolutionPilotModelSB21Calcioferrite`, etc. (Layer 3) |
| MAGEMin branch persistence / warm-start cache | `GlobalizedSolidSolutionState` + `cachedInternalx` (Layer 2) |
| MAGEMin database parameter tables | `MAGEMinImportedConstrainedTernarySolutionThermoModel` + family-specific thermo functions (Layer 3) |
