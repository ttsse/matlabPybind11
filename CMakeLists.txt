# CMakeLists.txt -- Cmake source file to build shared C++ library 'cppFuncs'
# Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.

cmake_minimum_required(VERSION 4.0)
## Set path to conda environment and build directory
set(CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")

### Necessary for win64.  See https://cmake.org/cmake/help/latest/variable/CMAKE_MSVC_RUNTIME_LIBRARY.html
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")

### -fPIC flag (position independent code to build shared libraries)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

### Make sure clang is used on macOS
IF (APPLE)
    set(CMAKE_C_COMPILER "$ENV{CONDA_PREFIX}/bin/clang")
    set(CMAKE_CXX_COMPILER "$ENV{CONDA_PREFIX}/bin/clang++")
ENDIF()

### Initialize the project (C required for SuiteSparse)
project(cppFuncs LANGUAGES CXX C)

find_package(Python 3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 3.0.1 EXACT CONFIG REQUIRED)

pybind11_add_module(cppFuncs src/cppFuncs/bindings.cpp src/cppFuncs/solvers.cpp src/cppFuncs/funcs.hpp)

### Necessary for macOS. See https://pybind11.readthedocs.io/en/stable/compiling.html
IF (APPLE)
    target_link_options(cppFuncs PRIVATE
        "-undefined"
        "dynamic_lookup"
    ) 
ENDIF()

### Include Eigen headers
find_package(Eigen3 QUIET)
IF (Eigen3_FOUND)
    message(STATUS "Found Eigen3 version: ${Eigen3_VERSION}.")
    target_link_libraries(cppFuncs PRIVATE Eigen3::Eigen)
ELSE()
    message(STATUS "Eigen3 not found by find_package, trying manual instead.")
    find_path(EIGEN3_INCLUDE_DIRS
        NAMES "Eigen"
        PATH_SUFFIXES include Library/include eigen3
    )
    target_include_directories(cppFuncs PRIVATE ${EIGEN3_INCLUDE_DIRS})
ENDIF()

## Include OpenMP
find_package(OpenMP QUIET COMPONENTS C CXX)
IF (OpenMP_FOUND)
    message(STATUS "Found OpenMP")
    target_link_libraries(cppFuncs PRIVATE OpenMP::OpenMP_CXX OpenMP::OpenMP_C)
ELSE()
    find_library(OpenMP 
                 NAMES omp   
                 PATH_SUFFIXES Library/lib lib 
                 )
    find_PATH(OpenMP_INCLUDE_DIRS 
                 NAMES "omp.h"     
                 PATH_SUFFIXES Library/include include
                 )
    target_include_directories(cppFuncs PRIVATE ${OpenMP_INCLUDE_DIRS})
    target_link_libraries(cppFuncs PRIVATE ${OpenMP})
    message(STATUS "OpenMP include dirs: ${OpenMP_INCLUDE_DIRS}")
    message(STATUS "OpenMP: ${OpenMP}")
ENDIF()

### Include SuiteSparse 
find_package(SuiteSparse QUIET NAMES SuiteSparse_config SuiteSparse suitesparse)
IF (SUITESPARSE_CONFIG_VERSION)
    message(STATUS "Found SuiteSparse_config version: ${SUITESPARSE_CONFIG_VERSION}")
    target_link_libraries(cppFuncs PRIVATE SuiteSparse::SuiteSparseConfig)
    set(SuiteSparse_config_FOUND TRUE)
ELSE()
    message(STATUS "SuiteSparse_config not found by find_package, trying manual instead.")
    find_path(SUITESPARSE_CONFIG_INCLUDE_DIR
        NAMES "suitesparse_config.h"
        PATH_SUFFIXES include/suitesparse Library/include/suitesparse
    )
    find_library(SUITESPARSE_CONFIG_LIBRARY
        NAMES suitesparseconfig SuiteSparse_config
        PATH_SUFFIXES lib Library/lib
    )
    target_include_directories(cppFuncs PRIVATE ${SUITESPARSE_CONFIG_INCLUDE_DIR})
    target_link_libraries(cppFuncs PRIVATE ${SUITESPARSE_CONFIG_LIBRARY})
    set(SuiteSparse_config_FOUND TRUE)
    message(STATUS "SuiteSparse_config found: ${SUITESPARSE_CONFIG_LIBRARY}")
ENDIF()

### Include BLAS and LAPACK
find_path(SUITESPARSE_BLAS_DIR
    NAMES "SuiteSparseBLAS.cmake"
         "SuiteSparseLAPACK.cmake"
    PATH_SUFFIXES Library/lib/cmake/SuiteSparse lib/cmake/SuiteSparse
)
set(BLA_VENDOR "OpenBLAS") # Look specifically for OpenBLAS installed in the conda env
list(APPEND CMAKE_MODULE_PATH "${SUITESPARSE_BLAS_DIR}")
include (  "${SUITESPARSE_BLAS_DIR}/SuiteSparseBLAS.cmake" )     # requires cmake 3.22
include ( "${SUITESPARSE_BLAS_DIR}/SuiteSparseLAPACK.cmake" )   # requires cmake 3.22

IF (LAPACK_FOUND)
    target_link_libraries ( cppFuncs PRIVATE ${LAPACK_LIBRARIES} )
    target_include_directories ( cppFuncs PRIVATE ${LAPACK_INCLUDE_DIR} )
ELSE()
    find_path(LAPACK_INCLUDE_DIR
        NAMES "lapack.h"
        PATH_SUFFIXES include/openblas Library/include/openblas
    )
    find_library(LAPACK_LIBRARY
        NAMES openblas lapack
        PATH_SUFFIXES lib Library/lib
    )
    target_include_directories(cppFuncs PRIVATE ${LAPACK_INCLUDE_DIR})
    target_link_libraries(cppFuncs PRIVATE ${LAPACK_LIBRARY})
    set(LAPACK_FOUND TRUE)
ENDIF()

IF (BLAS_FOUND)
    target_link_libraries ( cppFuncs PRIVATE ${BLAS_LIBRARIES} )
    target_include_directories ( cppFuncs PRIVATE ${BLAS_INCLUDE_DIR} )
    message(STATUS "Found BLAS version: ${BLAS_INCLUDE_DIR}")
ELSE()
    find_path(BLAS_INCLUDE_DIR
        NAMES "cblas.h"
        PATH_SUFFIXES include/openblas Library/include/openblas
    )
    find_library(BLAS_LIBRARIES
        NAMES openblas blas
        PATH_SUFFIXES lib Library/lib
    )
    target_include_directories(cppFuncs PRIVATE ${BLAS_INCLUDE_DIR})
    target_link_libraries(cppFuncs PRIVATE ${BLAS_LIBRARIES})
    set(BLAS_FOUND TRUE)
ENDIF()

### Include SPQR
find_package(SPQR QUIET NAMES SPQR spqr)
IF (SPQR_FOUND AND COLAMD_FOUND AND AMD_FOUND AND CHOLMOD_FOUND AND CAMD_FOUND AND CCOLAMD_FOUND)
    message(STATUS "Found SPQR version: ${SPQR_VERSION}")
    target_link_libraries(cppFuncs PRIVATE SuiteSparse::SPQR SuiteSparse::CHOLMOD SuiteSparse::AMD SuiteSparse::COLAMD SuiteSparse::CAMD SuiteSparse::CCOLAMD)
ELSE()
    message(STATUS "SPQR, CHOLAMD, or AMD not found by find_package, trying manual instead.")
    find_path(SPQR_INCLUDE_DIR
        NAMES "SuiteSparseQR.hpp"
        PATH_SUFFIXES include/suitesparse Library/include/suitesparse
    )
    find_library(SPQR_LIBRARY
        NAMES spqr SPQR
        PATH_SUFFIXES lib Library/lib
    )
    message(STATUS "SPQR found: ${SPQR_LIBRARY}")
    find_library(COLAMD_LIBRARY
        NAMES colamd COLAMD
        PATH_SUFFIXES lib Library/lib
    )
    message(STATUS "COLAMD found: ${COLAMD_LIBRARY}")
    find_library(CHOLMOD_LIBRARY
        NAMES cholmod CHOLMOD
        PATH_SUFFIXES lib Library/lib
    )
    message(STATUS "CHOLMOD found: ${CHOLMOD_LIBRARY}")
    find_library(AMD_LIBRARY
        NAMES amd AMD
        PATH_SUFFIXES lib Library/lib
    )
    message(STATUS "AMD found: ${AMD_LIBRARY}")
    find_library(CAMD_LIBRARY
        NAMES camd CAMD
        PATH_SUFFIXES lib Library/lib
    )
    message(STATUS "CAMD found: ${CAMD_LIBRARY}")
    find_library(CCOLAMD_LIBRARY
        NAMES ccolamd CCHOLAMD
        PATH_SUFFIXES lib Library/lib
    )
    message(STATUS "CCOLAMD found: ${CCOLAMD_LIBRARY}")
    target_include_directories(cppFuncs PRIVATE ${SPQR_INCLUDE_DIR})
    target_link_libraries(cppFuncs PRIVATE ${SPQR_LIBRARY} ${COLAMD_LIBRARY} ${AMD_LIBRARY} ${CHOLMOD_LIBRARY} ${CAMD_LIBRARY} ${CCOLAMD_LIBRARY})
ENDIF()

### Remove and manually copy the built libary if crashes on build
add_custom_command(TARGET cppFuncs POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:cppFuncs> "${PROJECT_SOURCE_DIR}"
)
