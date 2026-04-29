cmake_minimum_required(VERSION 3.16)
project(PerpleXRegression CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(test_regression
    test_regression.cpp
    PerpleXFluidModel.cpp
    PerpleXGFSMModel.cpp
    PerpleXHybridEos.cpp
    PerpleXMrkMixture.cpp
    PerpleXMrkParameters.cpp
    PerpleXMrkPure.cpp
    PerpleXPureEos.cpp
    PerpleXElectrolyte.cpp
    PerpleXHKF.cpp
)

target_include_directories(test_regression PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(test_regression PRIVATE PERPLEX_STANDALONE=1)
