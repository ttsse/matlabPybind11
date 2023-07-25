/*
    solvers.cpp -- creating pybind11 bindings for Eigen library
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

using namespace std;
using namespace std::chrono;
namespace py = pybind11;

Eigen::VectorXd denseSolve(Eigen::MatrixXd A, Eigen::Ref<Eigen::VectorXd> b) {
    const int m = A.rows(), n = A.cols();
    if (m == n) {
        Eigen::PartialPivLU<Eigen::MatrixXd> lu;
        lu.compute(A);
        return lu.solve(b);
    } else {
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
        return qr.solve(b);
    }
}

Eigen::VectorXd sparseSolve(Eigen::SparseMatrix<double> A, Eigen::Ref<Eigen::VectorXd> b) {
    omp_set_num_threads(8);
    auto start = high_resolution_clock::now();
    Eigen::SPQR<Eigen::SparseMatrix<double>> spqr;
    spqr.compute(A);
    Eigen::VectorXd x = spqr.solve(b);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
 
    py::print("Time taken for sparse solve: ", static_cast<double>(duration.count()*1e-6));
    return x; 
}