cmake_minimum_required(VERSION 3.16)

# This directory is added from the parent DEW CMakeLists.txt
# so we don't need a separate project() here.

# ----------------------------------------------------------------------
# Catch2 via FetchContent
# ----------------------------------------------------------------------
include(FetchContent)

FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.7.1  # or newer if you like
)

FetchContent_MakeAvailable(catch2)

# ----------------------------------------------------------------------
# Test executables
# ----------------------------------------------------------------------
# Water properties regression tests
add_executable(test_dew_water
  test_dew_water.cpp
  WaterTestAdapters.cpp
)

# Reaction thermodynamics regression test
add_executable(test_dew_reaction
  test_dew_reaction.cpp
)

# Species parameters verification (quick diagnostic)
add_executable(verify_species_params
  verify_species_params.cpp
)

# ----------------------------------------------------------------------
# Include directories for the test targets
# ----------------------------------------------------------------------
target_include_directories(test_dew_water PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}          # tests/
  ${CMAKE_CURRENT_SOURCE_DIR}/..       # DEW/
  ${CMAKE_CURRENT_SOURCE_DIR}/../../.. # repo root -> Reaktoro/...
)

target_include_directories(test_dew_reaction PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}          # tests/
  ${CMAKE_CURRENT_SOURCE_DIR}/..       # DEW/
  ${CMAKE_CURRENT_SOURCE_DIR}/../../.. # repo root -> Reaktoro/...
)

target_include_directories(verify_species_params PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}          # tests/
  ${CMAKE_CURRENT_SOURCE_DIR}/..       # DEW/
  ${CMAKE_CURRENT_SOURCE_DIR}/../../.. # repo root -> Reaktoro/...
)

# ----------------------------------------------------------------------
# Compile definitions: pass path to the aqueous DB from parent
# ----------------------------------------------------------------------
target_compile_definitions(test_dew_water PRIVATE
  DEW_AQUEOUS_DB_PATH=\"${DEW_AQUEOUS_DB_PATH}\"
)

target_compile_definitions(test_dew_reaction PRIVATE
  DEW_AQUEOUS_DB_PATH=\"${DEW_AQUEOUS_DB_PATH}\"
)

target_compile_definitions(verify_species_params PRIVATE
  DEW_AQUEOUS_DB_PATH=\"${DEW_AQUEOUS_DB_PATH}\"
)

# ----------------------------------------------------------------------
# Link against DEW water library, Catch2, and yaml-cpp
# yaml-cpp is in conda environment, add its library directory
# ----------------------------------------------------------------------
if(DEFINED ENV{CONDA_PREFIX})
  target_link_directories(test_dew_water PRIVATE
    $ENV{CONDA_PREFIX}/Library/lib
  )
  target_include_directories(test_dew_water PRIVATE
    $ENV{CONDA_PREFIX}/Library/include
  )

  target_link_directories(test_dew_reaction PRIVATE
    $ENV{CONDA_PREFIX}/Library/lib
  )
  target_include_directories(test_dew_reaction PRIVATE
    $ENV{CONDA_PREFIX}/Library/include
  )

  target_link_directories(verify_species_params PRIVATE
    $ENV{CONDA_PREFIX}/Library/lib
  )
  target_include_directories(verify_species_params PRIVATE
    $ENV{CONDA_PREFIX}/Library/include
  )
endif()

# Check if main Reaktoro library exists, use it if available
# Try build-msvc first (Windows), then build (generic)
set(MAIN_REAKTORO_LIB_MSVC "${CMAKE_CURRENT_SOURCE_DIR}/../../../../build-msvc/Reaktoro/Release/Reaktoro.lib")
set(MAIN_REAKTORO_LIB "${CMAKE_CURRENT_SOURCE_DIR}/../../../../build/Reaktoro/Release/Reaktoro.lib")

if(EXISTS "${MAIN_REAKTORO_LIB_MSVC}")
  set(REAKTORO_LIB_PATH "${MAIN_REAKTORO_LIB_MSVC}")
  message(STATUS "Using main Reaktoro library: ${REAKTORO_LIB_PATH}")
elseif(EXISTS "${MAIN_REAKTORO_LIB}")
  set(REAKTORO_LIB_PATH "${MAIN_REAKTORO_LIB}")
  message(STATUS "Using main Reaktoro library: ${REAKTORO_LIB_PATH}")
else()
  set(REAKTORO_LIB_PATH "")
  message(WARNING "Main Reaktoro library not found. Tests may have linking errors.")
endif()

if(REAKTORO_LIB_PATH)
  target_link_libraries(test_dew_water PRIVATE
    dew_water
    Catch2::Catch2WithMain
    yaml-cpp
    "${REAKTORO_LIB_PATH}"
  )
  target_link_libraries(test_dew_reaction PRIVATE
    dew_water
    Catch2::Catch2WithMain
    yaml-cpp
    "${REAKTORO_LIB_PATH}"
  )
  target_link_libraries(verify_species_params PRIVATE
    dew_water
    yaml-cpp
    "${REAKTORO_LIB_PATH}"
  )
else()
  target_link_libraries(test_dew_water PRIVATE
    dew_water
    Catch2::Catch2WithMain
    yaml-cpp
  )
  target_link_libraries(test_dew_reaction PRIVATE
    dew_water
    Catch2::Catch2WithMain
    yaml-cpp
  )
  target_link_libraries(verify_species_params PRIVATE
    dew_water
    yaml-cpp
  )
endif()# ----------------------------------------------------------------------
# Copy required DLLs to the test executable directories after build
# ----------------------------------------------------------------------
# Copy Reaktoro.dll from main build (try build-msvc first, then build)
set(MAIN_REAKTORO_DLL_MSVC "${CMAKE_CURRENT_SOURCE_DIR}/../../../../build-msvc/Reaktoro/Release/Reaktoro.dll")
set(MAIN_REAKTORO_DLL "${CMAKE_CURRENT_SOURCE_DIR}/../../../../build/Reaktoro/Release/Reaktoro.dll")

if(EXISTS "${MAIN_REAKTORO_DLL_MSVC}")
  set(REAKTORO_DLL_PATH "${MAIN_REAKTORO_DLL_MSVC}")
elseif(EXISTS "${MAIN_REAKTORO_DLL}")
  set(REAKTORO_DLL_PATH "${MAIN_REAKTORO_DLL}")
else()
  set(REAKTORO_DLL_PATH "")
endif()

if(REAKTORO_DLL_PATH)
  add_custom_command(TARGET test_dew_water POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${REAKTORO_DLL_PATH}"
      "$<TARGET_FILE_DIR:test_dew_water>"
    COMMENT "Copying Reaktoro.dll to test_dew_water directory"
  )
  add_custom_command(TARGET test_dew_reaction POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${REAKTORO_DLL_PATH}"
      "$<TARGET_FILE_DIR:test_dew_reaction>"
    COMMENT "Copying Reaktoro.dll to test_dew_reaction directory"
  )
  add_custom_command(TARGET verify_species_params POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${REAKTORO_DLL_PATH}"
      "$<TARGET_FILE_DIR:verify_species_params>"
    COMMENT "Copying Reaktoro.dll to verify_species_params directory"
  )
endif()

# Copy yaml-cpp.dll from conda environment
if(DEFINED ENV{CONDA_PREFIX})
  if(WIN32)
    set(YAML_CPP_DLL "$ENV{CONDA_PREFIX}/Library/bin/yaml-cpp.dll")
  else()
    # On Linux/Mac, shared libraries are typically found via rpath
    set(YAML_CPP_DLL "$ENV{CONDA_PREFIX}/lib/libyaml-cpp.so")
  endif()

  if(EXISTS "${YAML_CPP_DLL}")
    add_custom_command(TARGET test_dew_water POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${YAML_CPP_DLL}"
        "$<TARGET_FILE_DIR:test_dew_water>"
      COMMENT "Copying yaml-cpp library to test_dew_water directory"
    )
    add_custom_command(TARGET test_dew_reaction POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${YAML_CPP_DLL}"
        "$<TARGET_FILE_DIR:test_dew_reaction>"
      COMMENT "Copying yaml-cpp library to test_dew_reaction directory"
    )
  endif()
endif()

# ----------------------------------------------------------------------
# Register the tests with CTest
# We run them in the tests/ folder so CSV files are found directly.
# ----------------------------------------------------------------------
add_test(
  NAME dew_water_regression
  COMMAND test_dew_water
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
  NAME dew_reaction_regression
  COMMAND test_dew_reaction
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# ----------------------------------------------------------------------
# Copy test data (CSV / XLSX) to test executable directories so running
# the exe directly from the build tree still finds the truth files.
# This copies all .csv and .xlsx files present in the source tests dir
# into the target output directories for the test executables.
# ----------------------------------------------------------------------
file(GLOB DEW_TEST_DATA
  "${CMAKE_CURRENT_SOURCE_DIR}/*.csv"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.xlsx"
)

foreach(_data_file IN LISTS DEW_TEST_DATA)
  get_filename_component(_data_name ${_data_file} NAME)

  # Copy for test_dew_water
  add_custom_command(TARGET test_dew_water POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${_data_file}"
      "$<TARGET_FILE_DIR:test_dew_water>/${_data_name}"
    COMMENT "Copying DEW test data ${_data_name} to test_dew_water dir"
  )

  # Copy for test_dew_reaction
  add_custom_command(TARGET test_dew_reaction POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${_data_file}"
      "$<TARGET_FILE_DIR:test_dew_reaction>/${_data_name}"
    COMMENT "Copying DEW test data ${_data_name} to test_dew_reaction dir"
  )
endforeach()
