#==============================================================================
#
# THREE-TIER TEST ORGANIZATION WITH CTEST INTEGRATION
#
# This file orchestrates Reaktoro's comprehensive test suite organized into
# three tiers for different purposes:
#
#   TIER 1 (CTest/C++): Core unit and deterministic integration tests
#     - Compiled Catch2-based C++ tests (fast, ~minutes)
#     - Perple_X regression tests vs. committed baseline data (deterministic)
#     - CI Usage: `ctest -L "core"` for PR gating (fast)
#
#   TIER 2 (Python): Python API and binding tests
#     - Python wrapper/binding tests (medium, ~seconds to minutes)
#     - Run via pytest through CTest with proper environment setup
#     - CI Usage: `ctest -L "python"` to verify binding layer only
#
#   TIER 3 (External): Long-running parity and regression workflows
#     - DEW benchmark regression suite (slow, ~hours)
#     - Requires external Perple_X tools or generated reference data
#     - CI Usage: `ctest -L "external"` or exclude with `-L "^external$"`
#       on limited-resource environments
#
# CTest Label Categories:
#   "core"       : Fast unit tests essential for PR gating (C++ + Python)
#   "unit"       : Pure unit tests (deterministic, no I/O)
#   "integration": Integration tests (may have external dependencies)
#   "python"     : Python API and binding tests
#   "regression" : Regression tests (verify against baselines)
#   "parity"     : Parity checks (verify consistency between implementations)
#   "external"   : Tests requiring external tools/data not in repo
#   "optional"   : Tests that are informational but not gating
#   "nightly"    : Long-running tests suitable for scheduled runs
#
#==============================================================================

# Configure local pytest.ini file to permit pytest execution without path specification
configure_file(pytest.ini.in ${CMAKE_BINARY_DIR}/pytest.ini)

# Compose a dedicated runtime PATH for tests so CTest can locate external DLL dependencies.
set(REAKTORO_TEST_RUNTIME_DIRS "${REAKTORO_LIBRARY_DIR}")

foreach(_dep_dir
    ThermoFun_DIR
    Optima_DIR
    phreeqc4rkt_DIR
    yaml-cpp_DIR
    spdlog_DIR
    fmt_DIR)
    if(DEFINED ${_dep_dir} AND NOT "${${_dep_dir}}" STREQUAL "")
        get_filename_component(_dep_prefix "${${_dep_dir}}/../../.." ABSOLUTE)
        if(WIN32 AND EXISTS "${_dep_prefix}/bin")
            list(APPEND REAKTORO_TEST_RUNTIME_DIRS "${_dep_prefix}/bin")
        elseif(EXISTS "${_dep_prefix}/lib")
            list(APPEND REAKTORO_TEST_RUNTIME_DIRS "${_dep_prefix}/lib")
        endif()
    endif()
endforeach()

list(REMOVE_DUPLICATES REAKTORO_TEST_RUNTIME_DIRS)

if(WIN32)
    string(JOIN ";" REAKTORO_TEST_RUNTIME_PATH ${REAKTORO_TEST_RUNTIME_DIRS})
    # Keep runtime resolution deterministic under VS multi-config CTest by using the Release output dir.
    set(REAKTORO_TEST_PATH "$<SHELL_PATH:${PROJECT_BINARY_DIR}/Reaktoro/Release>\;${REAKTORO_TEST_RUNTIME_PATH}\;$ENV{PATH}")
else()
    string(JOIN ":" REAKTORO_TEST_RUNTIME_PATH ${REAKTORO_TEST_RUNTIME_DIRS})
    set(REAKTORO_TEST_PATH "${REAKTORO_TEST_RUNTIME_PATH}:$ENV{PATH}")
endif()

# Some script-style tests import the extension as top-level `reaktoro4py`.
# Ensure build output dir, package root, and package subdir are discoverable.
if(WIN32)
    set(REAKTORO_CTEST_PYTHONPATH
        "$<SHELL_PATH:${PROJECT_BINARY_DIR}/Reaktoro/Release>\;${REAKTORO_PYTHONPATH}\;${REAKTORO_PYTHONPATH}/reaktoro")
else()
    set(REAKTORO_CTEST_PYTHONPATH
        "${PROJECT_BINARY_DIR}/Reaktoro/Release:${REAKTORO_PYTHONPATH}:${REAKTORO_PYTHONPATH}/reaktoro")
endif()

# Resolve the Python interpreter used by tests so CTest does not accidentally
# pick an incompatible interpreter from PATH (e.g., Python 3.13 with cp312 module).
if(DEFINED PYTHON_EXECUTABLE AND NOT "${PYTHON_EXECUTABLE}" STREQUAL "")
    set(REAKTORO_TEST_PYTHON "${PYTHON_EXECUTABLE}")
elseif(DEFINED Python_EXECUTABLE AND NOT "${Python_EXECUTABLE}" STREQUAL "")
    set(REAKTORO_TEST_PYTHON "${Python_EXECUTABLE}")
elseif(DEFINED Python3_EXECUTABLE AND NOT "${Python3_EXECUTABLE}" STREQUAL "")
    set(REAKTORO_TEST_PYTHON "${Python3_EXECUTABLE}")
else()
    set(REAKTORO_TEST_PYTHON python)
endif()

# Create cache variable with the number of parallel jobs to be used when executing pytest
set(REAKTORO_PYTEST_PARALLEL_JOBS 1 CACHE STRING "The number of parallel jobs when executing pytest.")

# Create auxiliary variable with parallel specification for pytest command
if(REAKTORO_PYTEST_PARALLEL_JOBS MATCHES 1)
    set(REAKTORO_PYTEST_PARALLEL_OPTION "")  # avoid passing "-n 1" to pytest to avoid overhead in creating a single parallel job
else()
    set(REAKTORO_PYTEST_PARALLEL_OPTION "-n ${REAKTORO_PYTEST_PARALLEL_JOBS}")
endif()

#==============================================================================
#
# TIER 1: CTest/C++ TESTS (Core Unit and Deterministic Integration)
#
#==============================================================================

# Create target `tests-cpp` to execute C++ tests
add_custom_target(tests-cpp
    DEPENDS reaktoro-cpptests
    COMMENT "Running C++ tests..."
    COMMAND ${CMAKE_COMMAND} -E env
        "PATH=${REAKTORO_TEST_PATH}"
            $<TARGET_FILE:reaktoro-cpptests>
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

# Register the compiled Catch-based C++ test executable with CTest.
# These are core unit tests: fast, deterministic, essential for PR gating.
add_test(
    NAME reaktoro-cpptests
    COMMAND ${CMAKE_COMMAND} -E env
        "PATH=${REAKTORO_TEST_PATH}"
        "$<SHELL_PATH:${PROJECT_BINARY_DIR}/Reaktoro/Release/reaktoro-cpptests${CMAKE_EXECUTABLE_SUFFIX}>"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

set_tests_properties(reaktoro-cpptests PROPERTIES
    LABELS "unit;core"
    TIMEOUT 600)  # 10 minutes for full C++ test suite

# Optional: Add Perple_X C++ regression tests (Tier 1 - deterministic parity)
# These tests verify consistency against committed baseline data (CSV files in test/gfsm/)
# They use pure/binary/ternary fluid models and are completely self-contained.
set(PERPLEX_EXTENSION_DIR "${PROJECT_SOURCE_DIR}/Reaktoro/Extensions/Perple_X")
if(EXISTS "${PERPLEX_EXTENSION_DIR}/CMakeLists.txt")
    set(PERPLEX_TEST_DATA_DIR "${PERPLEX_EXTENSION_DIR}/test")
    if(EXISTS "${PERPLEX_TEST_DATA_DIR}/gfsm")
        message(STATUS "Perple_X baseline data found at ${PERPLEX_TEST_DATA_DIR}; enabling parity tests")

        # Note: Perple_X extension is a standalone project (not integrated into main build)
        # If you want to include it in the build, uncomment:
        # add_subdirectory(${PERPLEX_EXTENSION_DIR})
        # Then add:
        # add_test(NAME perplex-parity-gfsm
        #     COMMAND ${CMAKE_COMMAND} -E env "PATH=${REAKTORO_TEST_PATH}"
        #             $<TARGET_FILE:test_regression>
        #     WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
        # set_tests_properties(perplex-parity-gfsm PROPERTIES
        #     LABELS "parity;core;regression"
        #     TIMEOUT 300)
    endif()
endif()

#==============================================================================
#
# TIER 2: PYTHON TESTS (API and Binding Tests)
#
#==============================================================================

# Create target `tests-py` to execute Python tests (for backward compatibility with manual builds)
add_custom_target(tests-py
    DEPENDS reaktoro-setuptools
    COMMENT "Running Python tests..."
    COMMENT "Starting parallel execution of pytest..."
    COMMAND ${CMAKE_COMMAND} -E env
        "PATH=${REAKTORO_PATH}"
        "PYTHONPATH=${REAKTORO_PYTHONPATH}"
        "PYTHONUTF8=1"
        "PYTHONIOENCODING=utf-8"
            ${REAKTORO_TEST_PYTHON} -m pytest -x ${REAKTORO_PYTEST_PARALLEL_OPTION}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

# Register Python binding tests with CTest for better CI integration
# These are core Python API tests: verify bindings and high-level functionality
add_test(
    NAME pytest-bindings-dew
    COMMAND ${CMAKE_COMMAND} -E env
        "PATH=${REAKTORO_PATH}"
        "PYTHONPATH=${REAKTORO_CTEST_PYTHONPATH}"
        "PYTHONUTF8=1"
        "PYTHONIOENCODING=utf-8"
            ${REAKTORO_TEST_PYTHON} -m pytest "${CMAKE_CURRENT_SOURCE_DIR}/bindings/test_dew_python.py" -v
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

set_tests_properties(pytest-bindings-dew PROPERTIES
    LABELS "python;core"
    TIMEOUT 120)

add_test(
    NAME pytest-bindings-standardthermo
    COMMAND ${CMAKE_COMMAND} -E env
        "PATH=${REAKTORO_PATH}"
        "PYTHONPATH=${REAKTORO_CTEST_PYTHONPATH}"
        "PYTHONUTF8=1"
        "PYTHONIOENCODING=utf-8"
            ${REAKTORO_TEST_PYTHON} "${CMAKE_CURRENT_SOURCE_DIR}/bindings/test_standardthermo_dew.py"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

set_tests_properties(pytest-bindings-standardthermo PROPERTIES
    LABELS "python;optional;regression"
    TIMEOUT 120)

add_test(
    NAME pytest-bindings-water-models
    COMMAND ${CMAKE_COMMAND} -E env
        "PATH=${REAKTORO_PATH}"
        "PYTHONPATH=${REAKTORO_CTEST_PYTHONPATH}"
        "PYTHONUTF8=1"
        "PYTHONIOENCODING=utf-8"
            ${REAKTORO_TEST_PYTHON} "${CMAKE_CURRENT_SOURCE_DIR}/bindings/test_water_model_combinations.py"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

set_tests_properties(pytest-bindings-water-models PROPERTIES
    LABELS "python;optional;regression"
    TIMEOUT 120)

add_test(
    NAME pytest-bindings-none-options
    COMMAND ${CMAKE_COMMAND} -E env
        "PATH=${REAKTORO_PATH}"
        "PYTHONPATH=${REAKTORO_CTEST_PYTHONPATH}"
        "PYTHONUTF8=1"
        "PYTHONIOENCODING=utf-8"
            ${REAKTORO_TEST_PYTHON} "${CMAKE_CURRENT_SOURCE_DIR}/bindings/test_none_options.py"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

set_tests_properties(pytest-bindings-none-options PROPERTIES
    LABELS "python;optional;regression"
    TIMEOUT 60)

# Optional: Full pytest suite covering all tests in the repository
# This can be used for comprehensive validation but may take longer
add_test(
    NAME pytest-full-suite
    COMMAND ${CMAKE_COMMAND} -E env
        "PATH=${REAKTORO_PATH}"
        "PYTHONPATH=${REAKTORO_CTEST_PYTHONPATH}"
        "PYTHONUTF8=1"
        "PYTHONIOENCODING=utf-8"
            ${REAKTORO_TEST_PYTHON} -m pytest ${REAKTORO_PYTEST_PARALLEL_OPTION}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

set_tests_properties(pytest-full-suite PROPERTIES
    LABELS "python;integration;optional;nightly"
    TIMEOUT 600)

#==============================================================================
#
# TIER 3: EXTERNAL PARITY AND REGRESSION WORKFLOWS (Optional/Nightly)
#
#==============================================================================

# DEW Benchmark Regression Suite
# These are long-running benchmark comparisons against DEW database predictions.
# They are optional and typically run in scheduled nightly CI builds.
set(DEW_BENCHMARK_DIR "${PROJECT_SOURCE_DIR}/DEW_Experimental_Benchmark")
if(NOT EXISTS "${DEW_BENCHMARK_DIR}/regression_suite.py")
    # Fallback to the in-repo Testing location used by the current layout.
    set(DEW_BENCHMARK_DIR "${PROJECT_SOURCE_DIR}/Testing")
endif()

if(EXISTS "${DEW_BENCHMARK_DIR}/regression_suite.py")

    # Smoke test: Quick subset of DEW regression cases (good for quick validation)
    add_test(
        NAME dew-regression-smoke
        COMMAND ${CMAKE_COMMAND} -E env
            "PATH=${REAKTORO_PATH}"
            "PYTHONPATH=${REAKTORO_CTEST_PYTHONPATH}"
            "PYTHONUTF8=1"
                ${Python_EXECUTABLE} "${DEW_BENCHMARK_DIR}/regression_suite.py"
                    --tags smoke --continue-on-fail
                    --python-exe ${Python_EXECUTABLE}
                    --reaktoro-lib-dir "${REAKTORO_LIBRARY_DIR}"
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

    set_tests_properties(dew-regression-smoke PROPERTIES
        LABELS "regression;dew;external;optional"
        TIMEOUT 1800)  # 30 minutes

    # Full DEW regression: All tagged benchmark cases
    add_test(
        NAME dew-regression-full
        COMMAND ${CMAKE_COMMAND} -E env
            "PATH=${REAKTORO_PATH}"
            "PYTHONPATH=${REAKTORO_CTEST_PYTHONPATH}"
            "PYTHONUTF8=1"
                ${Python_EXECUTABLE} "${DEW_BENCHMARK_DIR}/regression_suite.py"
                    --continue-on-fail
                    --python-exe ${Python_EXECUTABLE}
                    --reaktoro-lib-dir "${REAKTORO_LIBRARY_DIR}"
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

    set_tests_properties(dew-regression-full PROPERTIES
        LABELS "regression;dew;external;nightly"
        TIMEOUT 7200)  # 2 hours

    # Parity-only DEW/PerplexDEW workflows (optional, external, usually nightly)
    add_test(
        NAME dew-parity-smoke
        COMMAND ${CMAKE_COMMAND} -E env
            "PATH=${REAKTORO_PATH}"
            "PYTHONPATH=${REAKTORO_CTEST_PYTHONPATH}"
            "PYTHONUTF8=1"
                ${Python_EXECUTABLE} "${DEW_BENCHMARK_DIR}/regression_suite.py"
                    --tags parity --continue-on-fail --max-failures 1 --max-cases 1
                    --python-exe ${Python_EXECUTABLE}
                    --reaktoro-lib-dir "${REAKTORO_LIBRARY_DIR}"
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

    set_tests_properties(dew-parity-smoke PROPERTIES
        LABELS "parity;dew;external;optional"
        TIMEOUT 1800)

    add_test(
        NAME dew-parity-full
        COMMAND ${CMAKE_COMMAND} -E env
            "PATH=${REAKTORO_PATH}"
            "PYTHONPATH=${REAKTORO_CTEST_PYTHONPATH}"
            "PYTHONUTF8=1"
                ${Python_EXECUTABLE} "${DEW_BENCHMARK_DIR}/regression_suite.py"
                    --tags parity --continue-on-fail
                    --python-exe ${Python_EXECUTABLE}
                    --reaktoro-lib-dir "${REAKTORO_LIBRARY_DIR}"
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

    set_tests_properties(dew-parity-full PROPERTIES
        LABELS "parity;dew;external;nightly"
        TIMEOUT 7200)

else()
    message(STATUS "DEW regression suite not found at ${DEW_BENCHMARK_DIR}; skipping external regression tests")
endif()

#==============================================================================
#
# BACKWARD COMPATIBILITY: Custom Targets for Manual Builds
#
#==============================================================================

# Combine both targets `tests-cpp` and `tests-py` to create target `tests`
# for backward compatibility with existing build scripts and local development
add_custom_target(tests DEPENDS tests-py tests-cpp)

#==============================================================================
#
# C++ EXTENSION TESTS (tests/cpp/)
# Perple_X C++ regression tests — standalone, no external library dependency.
# DEW C++ tests are registered via Reaktoro/Extensions/DEW/CMakeLists.txt
# which owns the dew_water library.
#
#==============================================================================
add_subdirectory(cpp)
