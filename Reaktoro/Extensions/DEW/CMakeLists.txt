cmake_minimum_required(VERSION 3.16)

project(DEW_Water_Standalone LANGUAGES CXX)

# ----------------------------------------------------------------------
# Basic compiler settings
# ----------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to Release as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)
endif()

# Enable parallel build if MSVC is used
add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/MP>)

# Set compiler-specific flags
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_compile_options(
        /D_SCL_SECURE_NO_WARNINGS
        /D_CRT_SECURE_NO_WARNINGS=1
        /EHsc
        /DNOMINMAX
        /Zc:__cplusplus  # Ensure __cplusplus reports correct C++ version
    )
endif()

# Generate compile_commands.json for non-MSVC compilers (used by IDEs/clangd)
if(NOT MSVC)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# ----------------------------------------------------------------------
# Enable testing (CTest)
# ----------------------------------------------------------------------
include(CTest)

# ----------------------------------------------------------------------
# Fetch tsl-ordered-map dependency (header-only library)
# ----------------------------------------------------------------------
include(FetchContent)

FetchContent_Declare(
    tsl-ordered-map
    GIT_REPOSITORY https://github.com/Tessil/ordered-map.git
    GIT_TAG        v1.1.0
)

# Download and populate without including their CMakeLists.txt (avoids CMake version conflict)
FetchContent_GetProperties(tsl-ordered-map)
if(NOT tsl-ordered-map_POPULATED)
    FetchContent_Populate(tsl-ordered-map)

    # Create interface library manually (it's header-only)
    add_library(tsl-ordered-map INTERFACE)
    target_include_directories(tsl-ordered-map INTERFACE
        ${tsl-ordered-map_SOURCE_DIR}/include
    )
    add_library(tsl::ordered_map ALIAS tsl-ordered-map)
endif()

# ----------------------------------------------------------------------
# Find or setup required dependencies
# ----------------------------------------------------------------------
# Find autodiff (header-only library from conda environment)
if(DEFINED ENV{CONDA_PREFIX})
    set(CONDA_PREFIX $ENV{CONDA_PREFIX})
    message(STATUS "Using conda environment: ${CONDA_PREFIX}")
endif()

find_path(AUTODIFF_INCLUDE_DIR
    NAMES autodiff/forward/real.hpp
    PATHS
        ${CONDA_PREFIX}/include
        ${CONDA_PREFIX}/Library/include
        $ENV{CONDA_PREFIX}/include
        $ENV{CONDA_PREFIX}/Library/include
    DOC "Path to autodiff headers"
)

if(NOT AUTODIFF_INCLUDE_DIR)
    message(FATAL_ERROR
        "autodiff headers not found. Please ensure:\n"
        "  1. The conda environment with autodiff is activated, OR\n"
        "  2. Set AUTODIFF_INCLUDE_DIR to the autodiff include directory"
    )
endif()

message(STATUS "Found autodiff headers: ${AUTODIFF_INCLUDE_DIR}")

add_library(autodiff INTERFACE)
target_include_directories(autodiff SYSTEM INTERFACE ${AUTODIFF_INCLUDE_DIR})
add_library(autodiff::autodiff ALIAS autodiff)

# Find Eigen3 from conda environment or system
if(DEFINED CONDA_PREFIX)
    list(PREPEND CMAKE_PREFIX_PATH ${CONDA_PREFIX} ${CONDA_PREFIX}/Library)
endif()
find_package(Eigen3 REQUIRED)# ----------------------------------------------------------------------
# Path to DEW aqueous DB (relative to this DEW directory):
#   DEW dir: Reaktoro/Extensions/DEW
#   DB dir : embedded/databases/DEW/dew2024-aqueous.yaml
# ----------------------------------------------------------------------
set(DEW_AQUEOUS_DB_PATH
    "${CMAKE_CURRENT_SOURCE_DIR}/../../../embedded/databases/DEW/dew2024-aqueous.yaml"
    CACHE FILEPATH "Path to dew2024-aqueous.yaml")

# ----------------------------------------------------------------------
# Get repository root path
# ----------------------------------------------------------------------
get_filename_component(REPO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)

# ----------------------------------------------------------------------
# Build a DEW water library from all .cpp files in this folder
# ----------------------------------------------------------------------
add_library(dew_water
    WaterBornOmegaDEW.cpp
    WaterDielectricFernandez1997.cpp
    WaterDielectricFranck1990.cpp
    WaterDielectricJohnsonNorton.cpp
    WaterDielectricPowerFunction.cpp
    WaterDielectricModel.cpp
    WaterElectroModel.cpp
    WaterEosZhangDuan2005.cpp
    WaterEosZhangDuan2009.cpp
    WaterGibbsModel.cpp
    WaterModelOptions.cpp
    WaterPsatPolynomialsDEW.cpp
    WaterSolventFunctionDEW.cpp
    WaterState.cpp
    WaterThermoModel.cpp
)

target_include_directories(dew_water PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}            # DEW/*.hpp
    ${REPO_ROOT}                           # repo root -> Reaktoro/Common/..., Reaktoro/Extensions/DEW/...
)

target_link_libraries(dew_water PUBLIC
    tsl::ordered_map
    autodiff::autodiff
    Eigen3::Eigen
)

# Enable implicit conversion of autodiff::real to double (same as main Reaktoro build)
target_compile_definitions(dew_water PUBLIC AUTODIFF_ENABLE_IMPLICIT_CONVERSION_REAL=1)

# Set target properties for cross-platform compatibility
set_target_properties(dew_water PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CXX_EXTENSIONS OFF
)

# Apply MSVC-specific compile definition for Eigen (prevents long compile times)
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    target_compile_definitions(dew_water PRIVATE EIGEN_STRONG_INLINE=inline)
endif()


# ----------------------------------------------------------------------
# Ensure yaml-cpp can be found for tests (from conda environment)
# ----------------------------------------------------------------------
if(DEFINED ENV{CONDA_PREFIX})
    set(CONDA_PREFIX $ENV{CONDA_PREFIX})
    list(PREPEND CMAKE_PREFIX_PATH ${CONDA_PREFIX} ${CONDA_PREFIX}/Library)
    # Also add library search path for linking
    link_directories(${CONDA_PREFIX}/Library/lib)
endif()

# ----------------------------------------------------------------------
# Add tests subdirectory for standalone builds
# (When built as part of main Reaktoro, WaterDEW.test.cxx is auto-included)
# ----------------------------------------------------------------------
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Standalone build of DEW - add separate test executable
    add_subdirectory(tests)
else()
    # Built as part of Reaktoro - WaterDEW.test.cxx automatically included
    message(STATUS "DEW tests: WaterDEW.test.cxx will be included in main Reaktoro test suite")
endif()
