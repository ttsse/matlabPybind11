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
#include <eigen3/Eigen/SparseLU>
#include "funcs.hpp"
#include <thread>

using namespace std;
using namespace std::chrono;
namespace py = pybind11;

Eigen::VectorXd denseSolve(Eigen::MatrixXd A, Eigen::Ref<Eigen::VectorXd> b, bool debug=false) {
    // Use all available threads
    const auto threadNum = std::thread::hardware_concurrency();
    omp_set_num_threads(threadNum);
    Eigen::VectorXd x;
    const int m = A.rows(), n = A.cols();
    double timeTaken;
    if (m == n) {
        Eigen::PartialPivLU<Eigen::MatrixXd> lu;
        auto start = high_resolution_clock::now();
        lu.compute(A);
        x = lu.solve(b);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        timeTaken = static_cast<double>(duration.count()*1e-6);
    } else {
        auto start = high_resolution_clock::now();
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
        x = qr.solve(b);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        timeTaken = static_cast<double>(duration.count()*1e-6);
    };
    if (debug) {
        x.conservativeResize(x.size() + 1);
        x(x.size() - 1) = timeTaken;
    };
    return x;
}

Eigen::VectorXd sparseSolve(Eigen::SparseMatrix<double> A, Eigen::Ref<Eigen::VectorXd> b, bool debug=false) {
    // Use all available threads
    const auto threadNum = std::thread::hardware_concurrency();
    omp_set_num_threads(threadNum);
    Eigen::VectorXd x;
    const int m = A.rows(), n = A.cols();
    double timeTaken;
    if (m == n) {
        Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
        auto start = high_resolution_clock::now();
        lu.compute(A);
        x = lu.solve(b);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        timeTaken = static_cast<double>(duration.count()*1e-6);
    } else {
        Eigen::SPQR<Eigen::SparseMatrix<double>> spqr;
        auto start = high_resolution_clock::now();
        spqr.compute(A);
        x = spqr.solve(b);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        timeTaken = static_cast<double>(duration.count()*1e-6);
    };
    if (debug) {
        x.conservativeResize(x.size() + 1);
        x(x.size() - 1) = timeTaken;
    };
    return x;
}