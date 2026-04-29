# File Structure & Dependencies (GFSM-Only)

## Directory Structure

```
Reaktoro/Extensions/Perple_X/
├── Pure-Species EOS for GFSM (Explicit Speciation Space)
│   ├── PerpleXPureEos.hpp          [14 pure EOS functions]
│   ├── PerpleXPureEos.cpp          [HSMRK, CORK, PSEOS, Haar, ZD05, ZD09]
│   ├── PerpleXHybridEos.hpp        [Selection framework]
│   ├── PerpleXHybridEos.cpp        [hybEos() implementation]
│   │
│   └── Dependencies:
│       ├── PerpleXMrkPure.hpp/cpp  [MRK base for all]
│       └── PerpleXMrkMixture.hpp/cpp [Mixing framework]
│
├── GFSM Solution Model (Uses Pure EOS)
│   ├── PerpleXGFSMModel.hpp
│   ├── PerpleXGFSMModel.cpp
│   │
│   └── Dependencies:
│       ├── PerpleXPureEos.hpp/cpp          [Pure EOS selection]
│       ├── PerpleXHybridEos.hpp/cpp        [Apply pure correction]
│       ├── PerpleXFluidModel.hpp/cpp       [Fluid properties]
│       └── PerpleXMrkMixture.hpp/cpp       [MRK base]
│
├── Support Components
│   ├── PerpleXMrkPure.hpp/cpp         [Pure species MRK]
│   ├── PerpleXMrkMixture.hpp/cpp      [Mixture MRK law]
│   ├── PerpleXMrkParameters.hpp/cpp   [MRK parameters]
│   ├── PerpleXFluidModel.hpp/cpp      [Fluid properties]
│   ├── PerpleXSpecies.hpp             [Species definitions]
│   ├── PerpleXHKF.hpp/cpp             [Aqueous HKF model]
│   └── PerpleXElectrolyte.hpp/cpp     [Electrolyte model]
│
└── Documentation
    ├── GFSM_SINGLE_SUMMARY.md                [Primary GFSM summary]
    ├── GFSM_EOS_COMPLETE_INVENTORY.md        [Detailed EOS inventory]
    └── [other documentation files]
```

## Dependency Graph

### Pure-Species EOS for GFSM (Explicit Speciation Space)

```
User Code (GFSM Calculations)
    ↓
    ├─→ PerpleXGFSMModel
    │       ├─→ PerpleXPureEos [HSMRK, CORK, PSEOS, Haar, ZD05, ZD09]
    │       ├─→ PerpleXHybridEos [hybEos() selection framework]
    │       │       └─→ PerpleXMrkPure [Base MRK pure]
    │       ├─→ PerpleXMrkMixture [MRK mixing law]
    │       └─→ PerpleXFluidModel [Fluid properties]
    │
    └─→ All depend on:
        ├─ PerpleXMrkPure
        ├─ PerpleXMrkParameters
        └─ PerpleXMrkMixture
```

## Header Include Structure

### PerpleXPureEos.hpp includes:
```cpp
#include "PerpleXHybridEos.hpp"
```
→ Provides `HybridEosOptions` struct for pure EOS selection

### PerpleXHybridEos.hpp includes:
```cpp
#include <array>
#include <functional>
#include <vector>
```
→ No circular dependency.

### PerpleXGFSMModel.hpp includes:
```cpp
#include "PerpleXPureEos.hpp"
#include "PerpleXHybridEos.hpp"
#include "PerpleXFluidModel.hpp"
#include "PerpleXMrkMixture.hpp"
```
→ Uses pure EOS framework.

---

## Compilation Units

### Always Compile
- PerpleXMrkPure.cpp - Base MRK
- PerpleXMrkMixture.cpp - Mixing law
- PerpleXMrkParameters.cpp - Parameters

### For GFSM Support
- PerpleXPureEos.cpp - Pure EOS implementations
- PerpleXHybridEos.cpp - Selection framework
- PerpleXGFSMModel.cpp - GFSM model implementation
- PerpleXFluidModel.cpp - Fluid properties

---

## Testing Considerations

### GFSM Pure-Species Tests
- Test each pure EOS option (H2O, CO2, CH4)
- Test mixed EOS combinations
- Verify against Perple_X reference data

---

## Migration Checklist for Developers

- Include PerpleXPureEos.hpp and PerpleXHybridEos.hpp
- Read PerpleXHybridEos.hpp header for architecture
- Read PerpleXPureEos.hpp header for pure EOS options
