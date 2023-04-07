#pragma once

#include <Eigen/Dense>

namespace cspl {
    // Interface for the polynomial class.
    template <unsigned int D, typename Polynomial>
    class ICubicHermitePolynomial {
    public:
        using VectorD = Eigen::Matrix<double, D, 1>; // D dimensional Vector.
        using VectorX = Eigen::Matrix<double, -1, 1>; // X dimensional Vector.
        using MatrixX = Eigen::Matrix<double, -1, -1>; // M by N dimensional Matrix.

        VectorD position(double t) const // Get position at normalised time t in [0, 1].
        {
            return static_cast<const Polynomial*>(this)->position(t);
        }

        VectorD velocity(double t) const // Get velocity at normalised time t in [0, 1].
        {
            return static_cast<const Polynomial*>(this)->velocity(t);
        }

        VectorD acceleration(double t) const // Get acceleration at normalised time t in [0, 1].
        {
            return static_cast<const Polynomial*>(this)->acceleration(t);
        }

        VectorX coeffs() const // Get polynomial coefficients.
        {
            return static_cast<const Polynomial*>(this)->coeffs();
        }

        VectorX points_all() const // Get polynomial points (initial, final).
        {
            return static_cast<const Polynomial*>(this)->points_all();
        }

        VectorX points_initial() const // Get initial polynomial points (initial, final).
        {
            return static_cast<const Polynomial*>(this)->points_initial();
        }

        VectorX points_target() const // Get final polynomial points (initial, final).
        {
            return static_cast<const Polynomial*>(this)->points_target();
        }

        void set_points(const VectorX& x) // Set the polynomial points (initial, final).
        {
            static_cast<Polynomial*>(this)->set_points(x);
        }

        void set_coeffs(const VectorX& x) // Set the polynomial coefficients.
        {
            static_cast<Polynomial*>(this)->set_coeffs(x);
        }

        MatrixX jac_pos(double t) const // Get Jacobian of position.
        {
            return static_cast<const Polynomial*>(this)->jac_pos(t);
        }

        MatrixX jac_vel(double t) const // Get Jacobian of velocity.
        {
            return static_cast<const Polynomial*>(this)->jac_vel(t);
        }

        MatrixX jac_acc(double t) const // Get Jacobian of acceleration.
        {
            return static_cast<const Polynomial*>(this)->jac_acc(t);
        }
    };
} // namespace cspl