#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

Eigen::VectorXd denseSolve(Eigen::MatrixXd, Eigen::Ref<Eigen::VectorXd>);
Eigen::VectorXd sparseSolve(Eigen::SparseMatrix<double>, Eigen::Ref<Eigen::VectorXd>);