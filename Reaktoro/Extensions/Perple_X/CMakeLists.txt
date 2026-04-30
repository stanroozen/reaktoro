cmake_minimum_required(VERSION 3.16)
project(PerpleXRegression CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Note: test_regression.cpp, test_hkf.cpp, test_electrolyte.cpp have been moved
# to tests/cpp/perplex/ for consolidated test organization.
# To build/run them, use the tests/cpp/perplex/CMakeLists.txt directly or
# via the main build's add_subdirectory(tests/cpp).
