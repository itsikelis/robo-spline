#pragma once

#include <Eigen/Dense>

// Interface for the polynomial class.

namespace cspl {
    template <unsigned int D>
    class ICubicHermitePolynomial {
    public:
        using VectorD = Eigen::Matrix<double, D, 1>; // D dimensional Vector.
        using VectorX = Eigen::Matrix<double, -1, 1>; // X dimensional Vector.
        using MatrixX = Eigen::Matrix<double, -1, -1>; // M by N dimensional Matrix.

        virtual VectorD position(double) const = 0; // Get position at normalised time t in [0, 1].
        virtual VectorD velocity(double) const = 0; // Get velocity at normalised time t in [0, 1].
        virtual VectorD acceleration(double) const = 0; // Get acceleration at normalised time t in [0, 1].

        virtual VectorX coeffs() const = 0; // Get polynomial coefficients.

        virtual VectorX nodes() const = 0; // Get polynomial parameters (initial, final).                             /
        virtual VectorX nodes_initial() const = 0; // Get initial polynomial parameters (initial, final).             /
        virtual VectorX nodes_target() const = 0; // Get final polynomial parameters (initial, final).                /

        virtual void set_nodes(const VectorX&) = 0; // Set the polynomial parameters (initial, final).
        virtual void set_coeffs(const VectorX&) = 0; // Set the polynomial coefficients.                      /

        virtual MatrixX jac_pos(double) const = 0; // Get Jacobian of position.
        virtual MatrixX jac_vel(double) const = 0; // Get Jacobian of velocity.
        virtual MatrixX jac_acc(double) const = 0; // Get Jacobian of acceleration.
    };
} // namespace cspl