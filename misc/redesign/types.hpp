#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace rspl {
    using Vector = Eigen::Matrix<double, -1, 1>; // dynamic Vector.
    using JacobianDense = Eigen::Matrix<double, -1, -1>; // dynamic Matrix.
    using Jacobian = Eigen::SparseMatrix<double, Eigen::RowMajor>; // dynamic Sparse Matrix.
    using Time = double;
    using SplineIndex = size_t;
} // namespace rspl