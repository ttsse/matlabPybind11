/*
    funcs.hpp -- C++ function declarations (fully inclusive)
    Copyright (c) 2023 Andreas Michael <andreas.michael@it.uu.se>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

Eigen::VectorXd denseSolve(Eigen::MatrixXd, Eigen::Ref<Eigen::VectorXd>);
Eigen::VectorXd sparseSolve(Eigen::SparseMatrix<double>, Eigen::Ref<Eigen::VectorXd>);