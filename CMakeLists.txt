# string(CONCAT MULTI "${CONDA_PREFIX}" "/Library/bin/g++.exe")
# string(CONCAT MULTI2 "${CONDA_PREFIX}" "/Library/bin/gcc.exe")
# message(STATUS "MULTI: $ENV{CONDA_PREFIX}")

set(CMAKE_CXX_COMPILER "$ENV{CONDA_PREFIX}/Library/bin/g++.exe")
set(CMAKE_C_COMPILER "$ENV{CONDA_PREFIX}/Library/bin/gcc.exe")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

### Example
# cmake_minimum_required(VERSION 3.15)
# project(example LANGUAGES CXX)

# find_package(Python 3.9 COMPONENTS Interpreter Development REQUIRED)
# find_package(pybind11 3.0.1 EXACT CONFIG REQUIRED)
# pybind11_add_module(example SHARED example.cpp)

cmake_minimum_required(VERSION 3.15)
project(cppFuncs LANGUAGES CXX)

find_package(Python 3.9 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 3.0.1 EXACT CONFIG REQUIRED)

pybind11_add_module(cppFuncs SHARED cppFuncs/bindings.cpp cppFuncs/solvers.cpp cppFuncs/funcs.hpp)

### Include Eigen headers
set(EIGEN3_INCLUDE_DIRS "$ENV{CONDA_PREFIX}/Library/include/eigen3")
target_include_directories(cppFuncs PRIVATE ${EIGEN3_INCLUDE_DIRS})


set(ALL_INCLUDE_DIRS "$ENV{CONDA_PREFIX}/Library/include")
# target_include_directories(cppFuncs PRIVATE ${OMP_INCLUDE_DIRS})

# set(OMP_INCLUDE_DIRS "$ENV{CONDA_PREFIX}/Library/include/suitesparse")
# target_include_directories(cppFuncs PRIVATE ${OMP_INCLUDE_DIRS})

# ### Include OMP
find_package(OpenMP REQUIRED)
target_link_libraries(cppFuncs PRIVATE ${OpenMP_LIBRARIES} ${OpenMP_DEPENDENCIES})

# ### Include SPQR
find_library(SPQRlib SPQR PATHS "$ENV{CONDA_PREFIX}/Library/lib")
find_library(cholmodLib cholmod PATHS "$ENV{CONDA_PREFIX}/Library/lib")
set(SPQRlib_INCLUDE_DIRS "$ENV{CONDA_PREFIX}/Library/include/suitesparse")
set(cholmod_INCLUDE_DIRS "$ENV{CONDA_PREFIX}/Library/include/suitesparse/cholmod")

target_include_directories(cppFuncs PRIVATE ${ALL_INCLUDE_DIRS} ${SPQRlib_INCLUDE_DIRS} ${cholmod_INCLUDE_DIRS})
target_link_libraries(cppFuncs PRIVATE m ${OpenMP_LIBRARIES} ${OpenMP_DEPENDENCIES} ${SPQRlib} ${cholmodLib})

### Finding extrernal packages (if not found by CMake)
# find_library(MYLIB mylib PATHS /path/to/mylib/lib)
# set(MYLIB_INCLUDE_DIRS /path/to/mylib/include)
# target_include_directories(cppFuncs PRIVATE ${MYLIB_INCLUDE_DIRS})
# target_link_libraries(cppFuncs PRIVATE ${MYLIB})