/*
    bindings.cpp -- pybind11 bindings for all functions declared in funcs.hpp
    Copyright (c) 2023 Andreas Michael <andreas.michael@it.uu.se>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <omp.h>
#include <math.h>
#include <iostream>
#include <chrono>
#include <Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/OrderingMethods>
#include <eigen3/Eigen/SPQRSupport>
#include "funcs.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(cppFuncs, m) {
    m.doc() = "Bindings for Eigen functions";  // optional module docstring

    m.def("denseSolve", &denseSolve,
          "A"_a.noconvert(),
          "b"_a.noconvert(),
          "Solves dense systems Ax = b", py::return_value_policy::reference_internal);
    m.def("sparseSolve", &sparseSolve, "Solves sparse systems Ax = b",
          "A"_a,
          "b"_a.noconvert(),
          "Solves sparse systems Ax = b", py::return_value_policy::reference_internal);
}