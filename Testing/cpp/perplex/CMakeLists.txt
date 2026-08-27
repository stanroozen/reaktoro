cmake_minimum_required(VERSION 3.16)

set(PERPLEX_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../../../Reaktoro/Extensions/Perple_X")

add_executable(test_regression
    test_regression.cpp
    ${PERPLEX_SRC}/PerpleXFluidModel.cpp
    ${PERPLEX_SRC}/PerpleXGFSMModel.cpp
    ${PERPLEX_SRC}/PerpleXHybridEos.cpp
    ${PERPLEX_SRC}/PerpleXMrkMixture.cpp
    ${PERPLEX_SRC}/PerpleXMrkParameters.cpp
    ${PERPLEX_SRC}/PerpleXMrkPure.cpp
    ${PERPLEX_SRC}/PerpleXPureEos.cpp
    ${PERPLEX_SRC}/PerpleXElectrolyte.cpp
    ${PERPLEX_SRC}/PerpleXHKF.cpp
)

add_executable(test_hkf
    test_hkf.cpp
    ${PERPLEX_SRC}/PerpleXHKF.cpp
    ${PERPLEX_SRC}/PerpleXElectrolyte.cpp
)

add_executable(test_electrolyte
    test_electrolyte.cpp
    ${PERPLEX_SRC}/PerpleXFluidModel.cpp
    ${PERPLEX_SRC}/PerpleXGFSMModel.cpp
    ${PERPLEX_SRC}/PerpleXHybridEos.cpp
    ${PERPLEX_SRC}/PerpleXMrkMixture.cpp
    ${PERPLEX_SRC}/PerpleXMrkParameters.cpp
    ${PERPLEX_SRC}/PerpleXMrkPure.cpp
    ${PERPLEX_SRC}/PerpleXPureEos.cpp
    ${PERPLEX_SRC}/PerpleXElectrolyte.cpp
    ${PERPLEX_SRC}/PerpleXHKF.cpp
)

foreach(_target test_regression test_hkf test_electrolyte)
    target_include_directories(${_target} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}  # tests/cpp/perplex/
        ${PERPLEX_SRC}               # source headers
        ${PROJECT_SOURCE_DIR}        # Reaktoro/Common/... headers
    )
    target_compile_definitions(${_target} PRIVATE PERPLEX_STANDALONE=1)
    set_target_properties(${_target} PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
    )
endforeach()

# Register with CTest
add_test(
    NAME perplex-test-regression
    COMMAND test_regression
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
    NAME perplex-test-hkf
    COMMAND test_hkf
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
    NAME perplex-test-electrolyte
    COMMAND test_electrolyte
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(perplex-test-regression perplex-test-hkf perplex-test-electrolyte
    PROPERTIES LABELS "cpp;perplex;regression")
