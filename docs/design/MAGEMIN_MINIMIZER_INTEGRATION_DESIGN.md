# MAGEMin Minimizer Integration Design

## Goal

Integrate a MAGEMin-style gradient-based inner minimizer for solid solutions into Reaktoro without changing Reaktoro's outer equilibrium workflow, split/retry machinery, or phase-activity API.

The immediate target is the current MAGEMin pilot path, not a repo-wide solver replacement.

## Verified Baseline

The current green baseline already isolates the outer and inner responsibilities well enough to support a narrow integration.

Verified code paths:

- `Reaktoro/Models/ActivityModels/ActivityModelGlobalizedSolidSolution.hpp`
- `Reaktoro/Models/ActivityModels/ActivityModelGlobalizedSolidSolution.cpp`
- `Reaktoro/Models/ActivityModels/ActivityModelMAGEMinSolidSolutionPilot.hpp`
- `Reaktoro/Models/ActivityModels/ActivityModelMAGEMinSolidSolutionPilot.cpp`
- `Reaktoro/Equilibrium/EquilibriumUtils.cpp`

What is already true:

1. Reaktoro consumes reduced solid-solution models through `GlobalizedSolidSolutionModel = Fn<GlobalizedSolidSolutionOutput(GlobalizedSolidSolutionInput)>`.
2. The outer equilibrium logic only sees `GlobalizedSolidSolutionOutput`, especially branch metadata, cached state, and optional split requests.
3. The retry/rebuild loop in `equilibrateWithGlobalizedSolidSolutionSplits(...)` does not depend on the internal minimizer algorithm. It only reacts to split requests and rebuilt phase definitions.
4. The current MAGEMin pilot already behaves like an adapter layer on top of that seam.

This means the correct first integration point is inside the MAGEMin pilot's local minimization path, not in the outer equilibrium solver.

## Current Responsibility Split

### Outer Reaktoro responsibilities

Owned by the globalized solid-solution seam and equilibrium retry layer:

- visible composition input
- branch metadata
- candidate screening and stability policy plumbing
- split request publication
- branch duplication and rebuilt systems
- cached state reuse across activity-model evaluations

Critical code:

- `ActivityModelGlobalizedSolidSolution(...)`
- `ComposeGlobalizedSolidSolutionBranch(...)`
- `ApplyGlobalizedSolidSolutionSplitRequests(...)`
- `equilibrateWithGlobalizedSolidSolutionSplits(...)`

These should remain unchanged for the first minimizer integration experiment.

### Current inner minimizer responsibilities

Owned by the MAGEMin pilot:

- thermodynamic family data and callbacks
- candidate warm-start handling for ternary branches
- local constrained minimization over internal coordinates
- conversion of the minimized internal composition into `Gx`, `Hx`, `ln_g`, `ln_a`

Critical code:

- `selectConstrainedTernaryBranch(...)`
- `minimizeTernaryObjective(...)`
- `MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(...)`

Today, the actual local solve is:

1. build a `GlobalizedSolidSolutionInternalProblem`
2. define an objective from imported MAGEMin excess and ideal terms plus the external-composition penalty
3. call `MinimizeGlobalizedSolidSolutionInternalProblem(...)`

That is the narrowest replaceable unit.

## Important Design Constraint

Do not replace `MinimizeGlobalizedSolidSolutionInternalProblem(...)` globally as the first step.

Reason:

- it is a seam-owned generic bounded search utility
- it is already reused outside the MAGEMin pilot
- changing it broadens the blast radius to non-MAGEMin models
- the user goal is to import MAGEMin's local minimization behavior while preserving Reaktoro functionality, not to replatform every internal solver immediately

The safer design is to add a MAGEMin-specific inner-minimizer hook in the pilot layer first.

## Smallest Viable Integration

### Step 1: Introduce a pilot-local minimizer strategy

Add a callable field to constrained ternary MAGEMin pilot options, conceptually:

```cpp
using MAGEMinConstrainedTernaryMinimizer = Fn<GlobalizedSolidSolutionInternalResult(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef visiblex,
    Optional<ArrayXr> warmstart)>;
```

Then extend `MAGEMinImportedConstrainedTernarySolutionOptions` with:

```cpp
MAGEMinConstrainedTernaryMinimizer minimizer;
```

Behavior:

- if `options.minimizer` is empty, use the current `minimizeTernaryObjective(...)`
- if `options.minimizer` is set, `selectConstrainedTernaryBranch(...)` uses it instead

This keeps the contract unchanged at the globalized seam while letting the MAGEMin pilot choose a different local solver.

### Step 2: Keep the current output contract unchanged

The new minimizer must still return `GlobalizedSolidSolutionInternalResult`, which already contains everything the pilot needs for branch selection diagnostics:

- `x`
- `objective`
- `iterations`
- `converged`

That means the branch screening, stability policy, and output assembly can stay structurally unchanged.

### Step 3: Continue to compute thermodynamic outputs in Reaktoro

For the first experiment, the new minimizer should only replace the search algorithm, not the surrounding thermodynamic bookkeeping.

Keep this in Reaktoro:

- `regularTernaryExcessChemicalPotentials(...)`
- `idealGibbs`
- `idealLnActivities`
- final `Gx`, `Hx`, `ln_g`, `ln_a` assembly

Reason:

- it limits the first experiment to "new local minimizer, same thermodynamic model"
- it separates solver validation from model-ingestion validation
- it preserves the current regression surface and fixture semantics

## First Implementation Sketch

Minimal code path change:

1. Add a minimizer function type and option field in `ActivityModelMAGEMinSolidSolutionPilot.hpp`.
2. Rename the current `minimizeTernaryObjective(...)` to something like `minimizeTernaryObjectiveWithDefaultSearch(...)`.
3. Add a small dispatcher:

```cpp
auto solveConstrainedTernaryInternalProblem(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult
{
    if(options.minimizer)
        return options.minimizer(options, T, x, warmstart);

    return minimizeTernaryObjectiveWithDefaultSearch(options, T, x, warmstart);
}
```

4. In `selectConstrainedTernaryBranch(...)`, replace the direct call to `minimizeTernaryObjective(...)` with that dispatcher.
5. Leave binary pilots unchanged for the first pass.

This is the smallest experiment that exercises the intended seam without disturbing outer equilibrium behavior.

## Why This Preserves Reaktoro Functionality

The following Reaktoro behaviors remain owned by existing code and therefore are not endangered by the first integration step:

- phase activity-model callback shape
- `ChemicalProps` caching behavior via `GlobalizedSolidSolutionState`
- branch ambiguity screening
- split request emission
- duplicate-phase rebuild workflow
- equilibrium retry and rebuild orchestration

Only the internal coordinate search changes.

That directly matches the project goal: import MAGEMin's local solution-phase minimization style while preserving Reaktoro's outer problem, diagnostics, and split-aware handling of immiscible solutions.

## What Not To Do First

Avoid these wider changes in the first experiment:

1. Do not replace Reaktoro's global equilibrium minimizer.
2. Do not change `GlobalizedSolidSolutionOutput` or `GlobalizedSolidSolutionInput` yet.
3. Do not collapse the globalized branch/split machinery into a MAGEMin-specific path.
4. Do not couple the first experiment to full MAGEMin database ingestion.
5. Do not generalize binary and ternary pilots simultaneously unless the ternary hook is already validated.

Each of those would confound solver integration with broader architectural change.

## Recommended Phased Plan

### Phase A: Solver hook only

Deliverable:

- pilot-local minimizer callback for constrained ternary imported MAGEMin models

Validation target:

- existing MAGEMin pilot regression fixtures remain green with the default minimizer path unchanged

### Phase B: One gradient-backed implementation

Deliverable:

- a new minimizer implementation for one constrained ternary family, preferably `sb11_cf` or `sb21_cf`

Why these are good first targets:

- ternary behavior already exists
- branch handling is meaningful
- fixture coverage already exists
- they exercise the local solver more honestly than a binary test case

Validation target:

- compare convergence, selected branch, and internal composition against the current baseline
- confirm retry/split tests still behave identically or explain any deliberate behavior shift

### Phase C: Generalize the pilot contract

If Phase B succeeds, then decide whether the minimizer callback should evolve into a fuller MAGEMin-backed local-model interface, for example if later work needs:

- analytic gradients
- family-specific constraints beyond simplex bounds
- richer solver diagnostics
- direct MAGEMin thermodynamic evaluations instead of imported callback closures

That is the point where a more general "local reduced model" abstraction may be justified.

## Proposed Success Criteria

The first integration should be considered successful if all of the following hold:

1. `ActivityModelGlobalizedSolidSolution(...)` call sites do not change.
2. `equilibrateWithGlobalizedSolidSolutionSplits(...)` does not change.
3. Existing MAGEMin pilot tests still pass with the default minimizer path.
4. A new MAGEMin-backed local minimizer can be selected through pilot options without changing the outer API.
5. The new minimizer returns the same shape of diagnostics currently consumed by branch selection.

## Suggested Initial Code Changes

If implementation starts from this note, the first edit slice should be:

- `Reaktoro/Models/ActivityModels/ActivityModelMAGEMinSolidSolutionPilot.hpp`
- `Reaktoro/Models/ActivityModels/ActivityModelMAGEMinSolidSolutionPilot.cpp`

Nothing else needs to change for the first solver-hook experiment.

## Test Surface To Reuse

Keep the current green baseline as the acceptance surface:

- `ActivityModelMAGEMinSolidSolutionPilotRegression.test.cxx`
- `ActivityModelGlobalizedSolidSolution.test.cxx`
- retry-related MAGEMin pilot filters on the trusted `build-msvc` path

The regression fixtures are especially valuable because they pin:

- selected thermodynamic family metadata
- internal composition
- `Gx`
- `Hx`
- `ln_g`
- `ln_a`
- split-request behavior

That is exactly the surface most likely to move if the inner minimizer changes behavior.

## Systematic Execution Status (2026-04-29)

This section records the current status of the proposed to-do sequence and the latest validation evidence.

### Completed

1. Phase A solver hook in pilot layer
    - Implemented pilot-local ternary minimizer callback type and option field.
    - Implemented dispatcher behavior: custom minimizer when provided, default strategy otherwise.
    - Wired constrained ternary branch selection through dispatcher path.

2. Outer seam preservation
    - Outer equilibrium split/retry orchestration was not changed for minimizer-hook integration.
    - Globalized seam contract remained unchanged.

3. Thermodynamic output contract preservation
    - Pilot path continues to assemble `Gx`, `Hx`, `ln_g`, and `ln_a` using existing bookkeeping.

4. Validation for no-regression and parity
    - Focused custom-hook regression passed:
      - filter: `*MAGEMin imported pilot custom ternary minimizer hook*`
      - result: 8 assertions, 1 test case
    - Projected-gradient versus legacy parity suite passed:
      - filter: `*projected-gradient and legacy minimizers agree*`
      - result: 47 assertions, 4 test cases
    - Retry/split fixture surface passed:
      - filter: `*MAGEMin equilibrium retry regression fixtures*`
      - result: 81 assertions, 1 test case
    - Globalized seam regression surface passed:
      - filter: `*ActivityModelGlobalizedSolidSolution*`
      - result: 23 assertions, 1 test case
    - Broad MAGEMin slice passed:
      - filter: `*MAGEMin*`
      - result: 392 assertions, 32 test cases

5. Optional hardening
    - Added dedicated negative-path regression for guarded projected-gradient disagreement and fallback selection.
    - Validation target passed:
      - filter: `*guarded projected-gradient falls back to legacy on forced disagreement*`
      - result: 9 assertions, 1 test case

6. Phase C minimal richer local-model contract
    - Added opt-in local-model contract type for constrained ternary pilots (`MAGEMinConstrainedTernaryLocalModel`), including objective and gradient callbacks plus bounds and solver settings.
    - Added opt-in local-model minimizer callback (`localModelMinimizer`) in constrained ternary options.
    - Preserved backward compatibility:
      - existing `minimizer` callback path remains unchanged
      - built-in legacy/projected-gradient behavior remains unchanged
      - new path is only used when explicitly configured
    - Added regression coverage for local-model callback selection and diagnostics strategy tagging (`custom-local-model`).

7. Extended hardening coverage
    - Added forced-disagreement fallback regressions for additional guarded families:
      - `sb11_pv`
      - `sb11_cf`
      - `sb21_nal`
    - Existing `sb11_ak` forced-disagreement fallback regression retained.

  8. Phase C analytic gradient path (public utility)
    - Exposed `MAGEMinProjectedGradientLocalModelMinimizer` as a public free function in the header.
    - Implemented the projected-gradient algorithm operating directly on `MAGEMinConstrainedTernaryLocalModel` using its `objective` and `gradient` callbacks.
    - External `localModelMinimizer` implementations can now delegate to this utility while injecting analytic or custom gradient callbacks.
    - Added regression coverage: `localModelMinimizer` lambda wraps the model gradient with a counter and delegates to the utility; verifies gradient invocation, strategy tag, and composition agreement with the default baseline.
    - Validation passed:
      - filter: `*MAGEMinProjectedGradientLocalModelMinimizer utility*`
      - result: 4 assertions, 1 test case
    - Broad MAGEMin slice still clean:
      - filter: `*MAGEMin*`
      - result: 442 assertions, 38 test cases

9. Extended hardening stress coverage over full ternary grids
    - Added benchmark-style negative-path grid-sweep fallback regression for guarded ternary families:
      - test: `Benchmarking forced-disagreement fallback over guarded ternary grids`
      - families covered: `sb11_pv`, `sb11_ak`, `sb11_cf`, `sb21_nal`
    - For each family, the test forces disagreement via `minimizerMaxIterations = 0` and asserts:
      - `comparedCount == totalEvaluations`
      - `fallbackCount == totalEvaluations`
      - `fallbackRate == 1.0`
      - `legacyCount == totalEvaluations`
      - `projectedGradientCount == 0`
    - Build validation in this session:
      - `reaktoro-cpptests` target rebuilt successfully on `build-msvc-pure`.
    - Runtime validation with corrected DLL/PATH setup:
      - PATH prepended with `C:\Users\stanroozen\anaconda3\envs\reaktoro\Library\bin` to satisfy `ThermoFun.dll` dependency.
      - filter: `*Benchmarking forced-disagreement fallback over guarded ternary grids*`
      - result: 28 assertions, 1 test case
      - filter: `*MAGEMin*`
      - result: 442 assertions, 38 test cases

10. Local-model contract compliance in projected-gradient utility
    - Updated `MAGEMinProjectedGradientLocalModelMinimizer` to respect `MAGEMinConstrainedTernaryLocalModel` contract fields:
      - `enforceUnityConstraint`
      - `lowerBounds`
      - `upperBounds`
    - Added bounded-simplex projection support for unity-constrained solves and box-clamping behavior for unconstrained-sum solves.
    - Added targeted regression coverage to lock the expected behavior for:
      - unconstrained-sum + box bounds
      - unity-constrained + bounded simplex
    - Validation status in this session:
      - compile reached source-complete stage for modified translation units
      - full link/test execution currently blocked by local environment dependency (`yaml-cpp.lib` not found in active CMake configuration)

11. Nonlinear constraint hooks and Hessian-ready fields on local-model contract
    - Extended `MAGEMinConstrainedTernaryLocalModel` struct with:
      - **Nonlinear constraint callbacks**: `constraints(y)` → m-vector of constraint values, `constraintJacobian(y)` → m×n matrix of constraint gradients
      - **Constraint bounds**: `constraintLowerBounds`, `constraintUpperBounds` for bound-constrained nonlinear solves
      - **Second-order information**: `objectiveHessian(y, multipliers)` → n×n Hessian of Lagrangian, `useSecondOrderInfo` flag for Newton/quasi-Newton methods
    - Added constraint helper utilities:
      - `hasNonlinearConstraints(model)` — detects presence of constraint callbacks
      - `validateConstraintCallbacks(model)` — ensures Jacobian is provided with constraints, validates Hessian consistency
      - `evaluateConstraintFeasibility(model, y)` — checks if composition satisfies all nonlinear constraints
    - Updated `MAGEMinProjectedGradientLocalModelMinimizer` to validate constraint callbacks at startup
    - Added comprehensive regression tests:
      - Constraint callback presence and validation logic
      - Constraint evaluation and feasibility checking
      - Hessian callback metadata validation
      - Mixed constraint scenarios (linear + nonlinear)
    - **Status**: Implementation complete; awaiting full test execution (current blocker: CMake build environment configuration)

12. Constraint-aware merit search and trust-region clipping in the projected-gradient utility
    - Checked upstream MAGEMin before implementation:
      - `MAGEMin/src/SB_database/SB_NLopt_opt_function.c` uses `NLOPT_LD_SLSQP` with declared lower/upper bounds plus a shared equality constraint `sum(x)=1` for SB models.
      - `MAGEMin/src/TC_database/NLopt_opt_function.c` extends that pattern with explicit NLopt inequality multi-constraints for site-fraction style constraints.
      - Upstream does **not** expose a custom trust-region kernel or bespoke backtracking line search in these local solvers; step control is delegated to NLopt.
      - Upstream diagnostics are stored back onto phase/global structs (`status`, `LM_time`, `df_raw`, `xeos`, `sf_ok`) rather than attached through family-specific callback payloads.
    - Reaktoro implementation choice for the projected-gradient local-model minimizer:
      - added a nonlinear-constraint merit function `objective + penalty*violation + barrier`
      - added optional feasible-only trial rejection during backtracking
      - added optional trust-region clipping on raw search displacements before projection
      - added an active-set SQP-style local search direction when nonlinear constraints are present:
        - builds a linearized active-constraint system from `constraintJacobian(y)` plus the unity constraint when enabled
        - solves a constrained steepest-descent step by default
        - upgrades to a KKT Newton step when `useSecondOrderInfo=true` and `objectiveHessian(y, multipliers)` is available
        - maps solved active multipliers back to the nonlinear constraint set for subsequent Hessian evaluations
    - Extended local-model contract fields:
      - `constraintPenaltyWeight`
      - `constraintBarrierWeight`
      - `requireFeasibleTrialPoints`
      - `trustRegionRadius`
    - Added focused regression coverage for:
      - constraint-aware backtracking that respects nonlinear constraint boundaries
      - Jacobian/Hessian callback usage on an active nonlinear constraint
      - trust-region clipping of projected-gradient trial steps
    - Family-specific diagnostics payload attachment required no new seam in this slice because `localModelDiagnostics` was already available and covered by regression.

### Remaining

1. Phase C contract evolution (optional)
  - Decide whether to extend the minimal local-model contract with additional hooks (for example: explicit constraints, Hessian information, or family-specific diagnostics payloads).
  - Only proceed if upcoming requirements need these capabilities.
  - **Status**: Constraint hooks, Hessian-ready fields, and a constraint-aware active-set SQP local step now exist; remaining work would be a more complete SQP globalization strategy or an NLopt-backed local solver, not more contract plumbing.
