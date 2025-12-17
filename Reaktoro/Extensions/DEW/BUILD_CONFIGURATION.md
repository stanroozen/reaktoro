# DEW Water Model - Build Configuration

## Build System

- **CMake Version**: 4.2.0
- **Generator**: Visual Studio 16 2019 (MSBuild)
- **Build Directory**: `C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\build-dew`
- **Source Directory**: `C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\Reaktoro\Extensions\DEW`

## Compiler & Toolchain

- **Compiler**: Microsoft Visual C++ (MSVC) 19.29.30153
- **Toolset**: MSVC 14.29.30133 (Visual Studio 2019 BuildTools)
- **Architecture**: x64 (64-bit)
- **C++ Standard**: C++17
- **Standard Required**: ON (strict)
- **Extensions**: OFF (disabled)

## Compiler Flags

- **Common**: `/DWIN32 /D_WINDOWS /GR /EHsc`
- **Debug**: `/Zi /Ob0 /Od /RTC1`
- **Release**: `/O2 /Ob2 /DNDEBUG` (optimization level 2)
- **MinSizeRel**: `/O1 /Ob1 /DNDEBUG`
- **RelWithDebInfo**: `/Zi /O2 /Ob1 /DNDEBUG`

## Dependencies (from Conda Environment `reaktoro`)

Located in: `C:\Users\stanroozen\anaconda3\envs\reaktoro\Library`

### 1. **autodiff** 1.1.2 (header-only)
   - Path: `${CONDA_PREFIX}/Library/include`
   - Configuration: `AUTODIFF_ENABLE_IMPLICIT_CONVERSION_REAL=1`

### 2. **Eigen** 3.4.0 (header-only)
   - Found via CMake's `find_package(Eigen3 REQUIRED)`

### 3. **yaml-cpp** 0.8.0 (compiled library)
   - Linked for tests
   - Path: `${CONDA_PREFIX}/Library/lib`

### 4. **tsl-ordered-map** 1.1.0 (header-only)
   - Fetched via CMake `FetchContent` from GitHub
   - Repository: https://github.com/Tessil/ordered-map.git

### 5. **Catch2** (testing framework)
   - Fetched automatically by test CMakeLists.txt

## Build Commands

### Configure (initial setup)
```powershell
cd C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro
mkdir build-dew
cd build-dew
cmake ..\Reaktoro\Extensions\DEW
```

### Build (compile)
```powershell
cd C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\build-dew
cmake --build . --config Release --target test_dew_water
```

### Run Tests
```powershell
cd C:\Users\stanroozen\Documents\Projects\reaktoro-dev\reaktoro\build-dew\tests\Release
.\test_dew_water.exe --reporter compact
```

## Key Configuration Details

### Library Target: `dew_water` (static library)
- Contains all DEW water model implementations
- Links: tsl::ordered_map, autodiff::autodiff, Eigen3::Eigen

### Test Target: `test_dew_water` (executable)
- Links: dew_water library, Catch2WithMain, yaml-cpp
- Location: `build-dew\tests\Release\test_dew_water.exe`

### Include Paths
- DEW source: `Reaktoro/Extensions/DEW`
- Repository root: `Reaktoro/...` (for Common utilities)
- Conda includes: autodiff, Eigen headers

### Critical Settings
- Testing enabled: `BUILD_TESTING=ON`
- Database path: `embedded/databases/DEW/dew2024-aqueous.yaml`
- Conda environment MUST be activated before building

## To Reproduce This Build

### 1. Activate conda environment
```powershell
conda activate reaktoro
```

### 2. Configure with CMake
```powershell
cmake -S Reaktoro/Extensions/DEW -B build-dew -G "Visual Studio 16 2019"
```

### 3. Build Release
```powershell
cmake --build build-dew --config Release
```

### 4. Run tests
```powershell
.\build-dew\tests\Release\test_dew_water.exe
```

## Test Results

✅ **All 22/22 tests passing (100%)**
- Total assertions validated: 16,348
- Configuration: Release build with optimization level O2

This configuration successfully compiles and passes all DEW water model tests!
