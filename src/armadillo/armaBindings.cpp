/*
    armaBindings.cpp -- Bindings for armadillo library 
    Copyright (c) 2023 Andreas Michael <andreas.michael@it.uu.se>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <omp.h>
#include <armadillo>
#include <vector>
#include <iostream>

PYBIND11_MAKE_OPAQUE(arma::mat);
PYBIND11_MAKE_OPAQUE(arma::vec);

namespace py = pybind11;

arma::mat np2arma(py::array_t<double, py::array::f_style | py::array::forcecast> pyArray) {
    // auto r = pyArray.mutable_unchecked<2>();
    // void * m_ptrA = A.request().ptr;
    uint n_rows, n_cols;
    auto buffer_info = pyArray.request();
    double *m_ptrA = static_cast<double *>(buffer_info.ptr);
    if (buffer_info.ndim != 1 && buffer_info.ndim != 2)
        throw std::runtime_error("Number of dimensions must be one or two");
    n_rows = pyArray.shape(0);
    n_cols = 1;
    if (buffer_info.ndim == 2) {
        n_rows = pyArray.shape(0);
        n_cols = pyArray.shape(1);
    }
    arma::Mat<double> A(m_ptrA, n_rows, n_cols, false, true);
    // check if pointers match
    // std::cout << "PtrA: " << static_cast<double *>(pyArray.request().ptr) << ", PtrB: " << A.memptr() << "\n";
    return A;
}

arma::vec solveLSystem(py::array_t<double, py::array::f_style | py::array::forcecast> pyArray1,
py::array_t<double, py::array::f_style | py::array::forcecast> pyArray2) {
    if (pyArray1.ndim() != 2 || pyArray2.ndim() != 1) {
        throw std::runtime_error("Incorrect arguments! Expected an MxN matrix and a Mx1 vector!");
    }
    if (pyArray2.shape(0) != pyArray1.shape(0)) {
        throw std::runtime_error("Incorrect arguments! Expected an MxN matrix and a Mx1 vector!");
    }
    arma::mat A = np2arma(pyArray1);
    arma::mat B = np2arma(pyArray2);
    return solve(A,B);
}

arma::vec solveSparseSystem(py::array_t<double, py::array::f_style | py::array::forcecast> pyArray1,
py::array_t<double, py::array::f_style | py::array::forcecast> pyArray2) {
    if (pyArray1.ndim() != 2 || pyArray2.ndim() != 1) {
        throw std::runtime_error("Incorrect arguments! Expected an NXN matrix and a Nx1 vector!");
    }
    if (pyArray2.shape(0) != pyArray1.shape(0)) {
        throw std::runtime_error("Incorrect arguments! Expected an NxN matrix and a Nx1 vector!");
    }
    arma::mat A = np2arma(pyArray1);
    arma::mat B = np2arma(pyArray2);
    return spsolve(arma::sp_mat(A),B);
}

PYBIND11_MODULE(armaFuncs, m) {
    m.doc() = "pybind11 example plugin";  // optional module docstring

    py::class_<arma::mat>(m, "mat", py::buffer_protocol())
        .def(py::init<const int, const int>())
        /* Acceptable input only numpy array with fortran ordering (column major), 
        np.array([[...][...]], order = 'F') */
        .def(py::init([](py::buffer b) {
            /* Request buffer from python */
            py::buffer_info info = b.request();

            /* Validation checks */
            if (info.format != py::format_descriptor<double>::format())
                throw std::runtime_error("Incompatible format: expected a double array!");

            if (info.ndim != 2)
                throw std::runtime_error("Incompatible buffer dimension!");

            double *m_ptrA = static_cast<double *>(info.ptr);
            arma::Mat<double> A(m_ptrA, info.shape[0], info.shape[1], false, true);

            return A;
        }))
        /* If A was converted from a row-major numpy array the converted numpy array will be incorrect. If it is a arma::mat object created here
        the conversion will work correctly. To be certain construct numpy arrays in column-major ordering always. */
        .def_buffer([](arma::mat &A) -> py::buffer_info {
            return py::buffer_info(
                A.memptr(),                             /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                2,                                      /* Number of dimensions */
                {   A.n_rows, A.n_cols},  /* Buffer dimensions */
                { sizeof(double),  /* Strides (in bytes) for each index */
                  sizeof(double)* static_cast<size_t>((A.n_rows))});
        })
        .def_readonly("n_rows", &arma::mat::n_rows)
        .def_readonly("n_cols", &arma::mat::n_cols)
        .def_readonly("n_elem", &arma::mat::n_elem)
        .def("__call__", [](arma::mat A, const int i, const int j) -> double {return A(i, j);}, py::is_operator())
        .def("__add__", [](const arma::mat A, const arma::mat B) -> arma::mat {return A + B;}, py::is_operator())
        .def("randu", [](arma::mat &A) -> void {A.randu();})
        .def("testFunc", [](const arma::mat &A, const arma::vec B) -> arma::vec {return B;}, py::is_operator())
        .def("__mul__", [](const arma::mat A, const arma::mat B) -> arma::mat {return A * B;}, py::is_operator())
        .def("__mul__", [](const arma::mat A, const arma::vec B) -> arma::vec {return A * B;}, py::is_operator())
        .def("print", [](arma::mat &A) -> void {A.print();}, "Print array")
        .def("print", [](arma::mat &A, std::string S) -> void {A.print(S);}, "Print array");

    py::class_<arma::vec>(m, "vec", py::buffer_protocol())
        .def(py::init<const int>())
        .def(py::init([](py::buffer b) {
            /* Request buffer from python */
            py::buffer_info info = b.request();

            /* Validation checks */
            if (info.format != py::format_descriptor<double>::format())
                throw std::runtime_error("Incompatible format: expected a double array!");

            if (info.ndim != 1)
                throw std::runtime_error("Incompatible buffer dimension!");

            double *m_ptrA = static_cast<double *>(info.ptr);
            arma::Mat<double> A(m_ptrA, info.shape[0], 1, false, true);

            return A;
        }))
        .def_buffer([](arma::vec &A) -> py::buffer_info {
            return py::buffer_info(
                A.memptr(),                             /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                1,                                      /* Number of dimensions */
                {   A.n_rows },  /* Buffer dimensions */
                { sizeof(double)});   /* Strides (in bytes) for each index */
        })
        .def("__call__", [](arma::vec A, const int i) -> double {return A(i);}, py::is_operator())
        .def("randu", [](arma::vec &A) -> void {A.randu();})
        .def("__add__", [](const arma::vec A, const arma::vec B) -> arma::vec {return A + B;}, py::is_operator())
        .def("print", [](arma::vec &A) -> void {A.print();}, "Print array")
        .def("print", [](arma::vec &A, std::string S) -> void {A.print(S);}, "Print array");

    m.def("np2arma", &np2arma, "A function that converts",py::return_value_policy::reference_internal);
    m.def("solveLSystem", &solveLSystem, "A function that solves a system",py::return_value_policy::reference_internal);
    m.def("solveSparseSystem", &solveSparseSystem, "A function that solves a sparse system",py::return_value_policy::reference_internal);
}
