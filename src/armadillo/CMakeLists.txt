# CMakeLists.txt -- Cmake source file to build shared C++ library 'armaFuncs' including Armadillo functions
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

### Initialize the project (C required for SuiteSparse)
project(armaFuncs LANGUAGES CXX C)

find_package(Python 3.9 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 3.0.1 EXACT CONFIG REQUIRED)

pybind11_add_module(armaFuncs SHARED armaBindings.cpp)

### Include Eigen headers
find_package(Armadillo QUIET)
IF (Armadillo_FOUND AND NOT TARGET Armadillo::Armadillo)
    message(STATUS "Found Armadillo version: ${Armadillo_VERSION}.")
    add_library(Armadillo::Armadillo INTERFACE IMPORTED)
    set_target_properties(
        Armadillo::Armadillo
        PROPERTIES
        INTERFACE_LINK_LIBRARIES "${ARMADILLO_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${ARMADILLO_INCLUDE_DIRS}"
    )
    target_link_libraries(armaFuncs PRIVATE Armadillo::Armadillo)
ELSE()
    message(STATUS "Armadillo not found by find_package, trying manual instead.")
    find_path(ARMA_INCLUDE_DIRS
        NAMES "armadillo"
        PATH_SUFFIXES include Library/include
    )
    find_library(ARMA_LIB 
                 NAMES armadillo   
                 PATH_SUFFIXES Library/lib lib 
                 )
    target_include_directories(armaFuncs PRIVATE ${ARMA_INCLUDE_DIRS})
    target_link_libraries(armaFuncs PRIVATE ${ARMA_LIB})
ENDIF()

### Remove and manually copy the built libary if crashes on build
add_custom_command(TARGET armaFuncs POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:armaFuncs> "${PROJECT_SOURCE_DIR}"
)
