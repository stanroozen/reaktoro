# Reaktoro Build Requirements

This document records the build and installation requirements for setting up this repository on another computer, including the exact dependency versions currently working in the local `reaktoro` Conda environment and the minimum versions enforced by CMake.

It covers two use cases:

1. Building the core Reaktoro C++ library and Python bindings.
2. Building the DEW and PerplexDEW-enabled local Python extension used by the benchmark and mineral solubility scripts in `DEW_Experimental_Benchmark/`.

## Scope And Current Baseline

- Repository version: `Reaktoro 2.13.0`
- Minimum CMake version: `3.17.0`
- Active development Python version: `3.12.12`
- Preferred environment file: `environment.devenv.yml`
- Minimal environment file: `environment.yml`

For the DEW benchmark workflows in this repository, the scripts currently expect a local build artifact in one of these folders:

- `build/Reaktoro/Release`

Several benchmark scripts also look for the locally built Python package in:

- `build/python/package`

## Operating System Requirements

### Windows

Required:

- 64-bit Windows
- Conda or Miniconda/Anaconda
- Visual Studio 2019 Build Tools or Visual Studio 2019 with C++ workload
- CMake
- Ninja (recommended)
- Git

Recommended Visual Studio components:

- MSVC v142 build tools
- Windows 10 or Windows 11 SDK
- C++ CMake tools for Windows

Notes:

- The Conda environment includes `vs2019_win-64 19.29.30139`, but you still need a working MSVC toolchain installed on the machine.
- In the current shell snapshot, `cl` was not available on `PATH`, so building on a fresh machine should be done from a Visual Studio Developer Prompt or a shell where MSVC has been initialized.
- For a reproducible Windows setup that creates the Conda environment and performs the build from one entry point, use `scripts/setup-windows.ps1`.

Scripted Windows setup:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/setup-windows.ps1
```

To build the DEW-oriented local artifact layout instead:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/setup-windows.ps1 -Mode dew
```

### Linux

Required:

- GCC toolchain or Clang toolchain
- Make or Ninja
- CMake
- Git

The repository environment file currently requests these Linux-specific packages:

- `gxx_linux-64`
- `clangxx`
- `lld`
- `make`
- `doxygen 1.9.1`
- `graphviz`
- `valgrind`

### macOS

Required:

- Xcode Command Line Tools

The repository environment file currently requests these macOS-specific packages:

- `clangxx_osx-64`
- `clangxx_osx-arm64`

## CMake-Enforced Dependency Requirements

These are the versions explicitly required by `cmake/ReaktoroFindDeps.cmake`.

### Required Build Dependencies
| Dependency | CMake requirement |
|---|---:|
| autodiff | `>= 1.1.1` |
| Eigen3 | `>= 3.4` |
| nlohmann_json | `>= 3.6.1` |
| Optima | `= 0.6.0` |
| phreeqc4rkt | `>= 3.6.2.1` |
| tabulate | `>= 1.4.0` |
| ThermoFun | `= 0.4.5` |
| tsl-ordered-map | `>= 1.0.0` |
| yaml-cpp | `>= 0.6.3` |
| Python | `>= 3.7` |
| pybind11 | `>= 2.10.0` |

### Optional Dependencies

| Dependency | CMake requirement | Purpose |
|---|---:|---|
| Catch2 | `>= 2.6.2` | C++ test build |
| reaktplot | `>= 0.4.1` | plotting helpers |
| openlibm | optional | alternative math library when `REAKTORO_ENABLE_OPENLIBM=ON` |
| pybind11-stubgen | optional | generate Python stubs |

## Exact Versions In The Working `reaktoro` Conda Environment

These are the exact versions installed in `C:\Users\stanroozen\anaconda3\envs\reaktoro` at the time of this update.

### Core Toolchain

| Package | Installed version |
|---|---:|
| python | `3.12.12` |
| cmake | `4.2.0` |
| ninja | `1.13.2` |
| ccache | `4.11.3` |
| vs2019_win-64 | `19.29.30139` |

Notes:

- `git --version` on `PATH` reported `2.51.0.windows.1`, while the Conda environment package version is `2.52.0`.
- `cl` was not available in the checked shell, so install and initialize MSVC separately on Windows.

### C++ / CMake Dependencies

| Package | Installed version |
|---|---:|
| autodiff | `1.1.2` |
| eigen | `3.4.0` |
| nlohmann_json | `3.12.0` |
| optima | `0.6.0` |
| phreeqc4rkt | `3.6.2.2` |
| thermofun | `0.4.5` |
| yaml-cpp | `0.8.0` |
| tsl_ordered_map | `1.1.0` |
| cpp-tabulate | `1.5` |
| openlibm | `0.6.0` |
| spdlog | `1.12.0` |

### Python Build / Test Dependencies

| Package | Installed version |
|---|---:|
| pybind11 | `2.13.6` |
| pybind11-abi | `4` |
| pybind11-stubgen | `2.0.2` |
| catch2 | `2.13.9` |
| pytest | `9.0.1` |
| pytest-regressions | `2.8.3` |
| pytest-xdist | `3.8.0` |
| reaktplot | `0.4.1` |

### Python Runtime Dependencies Used By Analysis And Benchmark Scripts

| Package | Installed version |
|---|---:|
| numpy | `2.3.5` |
| pandas | `2.3.3` |
| matplotlib-base | `3.10.8` |
| openpyxl | `3.1.5` |
| fire | `0.7.1` |
| ipykernel | `7.1.0` |
| nbformat | `5.10.4` |
| oyaml | `1.0` |

## Environment Files In The Repository

### `environment.devenv.yml`

This should be treated as the authoritative development environment file for this repository. It includes the full toolchain, test packages, notebook packages, and platform-specific compiler packages.

Notable pinned entries:

- `python = 3.12` via `PYTHON_VERSION` default
- `autodiff >= 1.1.1`
- `catch2 = 2`
- `doxygen = 1.9.1` on Linux
- `optima = 0.6.0`
- `pybind11 >= 2.10.0`
- `pybind11-stubgen = 2.0.2`
- `reaktplot >= 0.4.1`
- `thermofun = 0.4.5`
- `spdlog = 1.12`
- `yaml-cpp = 0.8.0`

### `environment.yml`

This file is the direct `conda env create -f environment.yml` fallback for machines where `conda-devenv` is not available. It should stay synchronized with the common packages in `environment.devenv.yml` that are needed for building and running the repository workflows.

`environment.devenv.yml` remains the authoritative source because it also carries platform selectors and a few extra development-only tools.

## Fresh Machine Setup

### Option 1: Recommended Conda Development Environment

Install `conda-devenv` in base if needed:

```bash
conda install -n base -c conda-forge conda-devenv=3.5.0
```

Create the environment from the repository root:

```bash
conda devenv -f environment.devenv.yml
conda activate reaktoro
```

### Option 2: Minimal Environment File

If `conda-devenv` is not available:

```bash
conda env create -f environment.yml
conda activate reaktoro
```

This is likely enough for a basic build, but not necessarily for all docs, tests, stubs, or benchmark workflows.

## Build Configurations

### Standard C++ And Python Build

From the repository root:

```bash
cmake -S . -B build \
  -DREAKTORO_BUILD_PYTHON=ON \
  -DREAKTORO_BUILD_TESTS=ON \
  -DREAKTORO_BUILD_DOCS=OFF

cmake --build build --config Release --parallel
cmake --install build --config Release
```

### Windows MSVC Build

Recommended from a Developer Prompt:

```bash
cmake -S . -B build -G "Visual Studio 16 2019" -A x64 \
  -DREAKTORO_BUILD_PYTHON=ON \
  -DREAKTORO_BUILD_TESTS=ON \
  -DREAKTORO_BUILD_DOCS=OFF

cmake --build build --config Release --parallel
cmake --install build --config Release
```

Equivalent checked-in script:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/setup-windows.ps1
```

### DEW / PerplexDEW Local Build

The benchmark scripts in `DEW_Experimental_Benchmark/` and `DEW_Experimental_Benchmark/Mineral_Solubilities/` currently rely on a local build folder named `build`.

Use:

```bash
cmake -S . -B build -G "Visual Studio 16 2019" -A x64 \
  -DREAKTORO_BUILD_PYTHON=ON \
  -DREAKTORO_BUILD_TESTS=ON \
  -DREAKTORO_BUILD_DOCS=OFF

cmake --build build --config Release --parallel
```

Why `build` matters:

- Many DEW scripts prepend `build/Reaktoro/Release` to `sys.path` when importing `reaktoro4py`.
- Some scripts also prepend `build/python/package`.
- The benchmark workflow now assumes a single canonical local build directory: `build`.

Equivalent checked-in script:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/setup-windows.ps1 -Mode dew
```

To run build, install, C++ tests, and the structured regression suite in one pass:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/setup-windows.ps1 -RunRegressionSuite -RegressionTags smoke
```

## Install Commands

The repository root includes an `install` CMake script that configures and installs into a build directory.

Default usage:

```bash
cmake -P install
```

Custom build directory:

```bash
cmake -DBUILD_PATH=build -P install
```

Custom install prefix:

```bash
cmake -DPREFIX=C:/local/reaktoro -P install
```

The `install` script:

- defaults to `Release`
- prefers `ninja` if found
- otherwise uses `make` if found

On Windows, direct `cmake -S/-B` usage is usually clearer than relying on the script.

## Validation After Install

### C++ Validation

If tests were built:

```bash
ctest --test-dir build --output-on-failure
```

For MSVC multi-config builds:

```bash
ctest --test-dir build -C Release --output-on-failure
```

### Python Validation

From the activated environment:

```bash
python -c "import reaktoro; print(reaktoro.__file__)"
```

For local extension validation:

```bash
python -c "import sys; sys.path.insert(0, r'build/Reaktoro/Release'); import reaktoro4py; print(reaktoro4py.__file__)"
```

### DEW Validation

From the repository root or benchmark folder:

```bash
python DEW_Experimental_Benchmark/run_quartz_PerplexDEW.py
```

Expected behavior:

- the script should report that it found a local `reaktoro4py` build
- it should load `DEWDatabase("dew2024-aqueous")`
- it should successfully configure `ActivityModelPerplexDEW`

### Structured Regression Validation

Run all smoke regression cases from a single file:

```bash
python DEW_Experimental_Benchmark/regression_suite.py --tags smoke --continue-on-fail
```

Run the broader DEW/PerplexDEW regression set:

```bash
python DEW_Experimental_Benchmark/regression_suite.py --tags regression --continue-on-fail
```

The suite writes per-case logs and a JSON summary to:

- `DEW_Experimental_Benchmark/regression_results/<timestamp>/summary.json`

## Optional External Requirement For Perple_X Comparison Workflows

The Reaktoro build itself does not require Perple_X.

However, if you want to reproduce the parity and comparison workflows that generate temporary `perplex_*` benchmark directories, you also need a working Perple_X installation and any local scripts or wrappers used to launch it.

That requirement is optional for building Reaktoro itself and optional for the normal Reaktoro Python package.

## Known Documentation Gaps This File Is Intended To Close

This file is needed because the repository currently has these gaps:

- the root `README.md` is intentionally brief and does not describe a fresh-machine build
- the `install` script does not explain the DEW-specific `build` workflow
- `environment.yml` is not the full development environment
- the benchmark scripts depend on local build folder conventions that are not otherwise documented

## Maintenance Rule For Future Updates

Whenever the active development environment changes, update all of the following together:

1. `environment.devenv.yml`
2. `environment.yml` if it is still intended to work
3. `cmake/ReaktoroFindDeps.cmake`
4. this file

If those files drift apart, the repository becomes harder to reproduce on a clean machine.
