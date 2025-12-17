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
# Test executable
# ----------------------------------------------------------------------
add_executable(test_dew_water
  test_dew_water.cpp
  WaterTestAdapters.cpp
)

# ----------------------------------------------------------------------
# Include directories for the test target
# ----------------------------------------------------------------------
target_include_directories(test_dew_water PRIVATE
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
endif()

# Check if main Reaktoro library exists, use it if available
set(MAIN_REAKTORO_LIB "${CMAKE_CURRENT_SOURCE_DIR}/../../../../build/Reaktoro/Release/Reaktoro.lib")
if(EXISTS "${MAIN_REAKTORO_LIB}")
  message(STATUS "Using main Reaktoro library: ${MAIN_REAKTORO_LIB}")
  target_link_libraries(test_dew_water PRIVATE
    dew_water
    Catch2::Catch2WithMain
    yaml-cpp
    "${MAIN_REAKTORO_LIB}"
  )
else()
  message(WARNING "Main Reaktoro library not found at ${MAIN_REAKTORO_LIB}. Tests may have linking errors.")
  target_link_libraries(test_dew_water PRIVATE
    dew_water
    Catch2::Catch2WithMain
    yaml-cpp
  )
endif()# ----------------------------------------------------------------------
# Copy required DLLs to the test executable directory after build
# ----------------------------------------------------------------------
# Copy Reaktoro.dll from main build
set(MAIN_REAKTORO_DLL "${CMAKE_CURRENT_SOURCE_DIR}/../../../../build/Reaktoro/Release/Reaktoro.dll")
if(EXISTS "${MAIN_REAKTORO_DLL}")
  add_custom_command(TARGET test_dew_water POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${MAIN_REAKTORO_DLL}"
      "$<TARGET_FILE_DIR:test_dew_water>"
    COMMENT "Copying Reaktoro.dll to test directory"
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
      COMMENT "Copying yaml-cpp library to test directory"
    )
  endif()
endif()

# ----------------------------------------------------------------------
# Register the test with CTest
# We run it in the tests/ folder so "truth_*.csv" are found directly.
# ----------------------------------------------------------------------
add_test(
  NAME dew_water_regression
  COMMAND test_dew_water
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
