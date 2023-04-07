#pragma once

#include <Eigen/Dense>

namespace cspl {
    // Interface for the polynomial class.
    template <unsigned int D>
    class CubicHermitePolynomial {
    public:
        using VectorD = Eigen::Matrix<double, D, 1>; // D dimensional Vector.
        using VectorX = Eigen::Matrix<double, -1, 1>; // X dimensional Vector.
        using MatrixX = Eigen::Matrix<double, -1, -1>; // M by N dimensional Matrix.

        // p0, v0 initial position and velocity, p1, v2: final position and velocity
        CubicHermitePolynomial(const VectorD& p0, const VectorD& v0, const VectorD& p1, const VectorD& v1)
        {
            VectorX p(D * 4);
            p << p0, v0, p1, v1;
            set_points(p);
        }

        // Get position at normalised time in [0-1].
        VectorD position(double t) const
        {
            return _c0 + (_c1 * t) + (_c2 * t * t) + (_c3 * t * t * t);
        }

        // Get velocity at normalised time in [0-1].
        VectorD velocity(double t) const
        {
            return _c1 + (2 * _c2 * t) + (3 * _c3 * t * t);
        }

        // Get acceleration at normalised time in [0-1].
        VectorD acceleration(double t) const
        {
            return (2 * _c2) + (6 * _c3 * t);
        }

        // Get polynomial coefficients.
        VectorX coeffs() const
        {
            VectorX coeffs(D * 4);
            coeffs << _c0, _c1, _c2, _c3;
            return coeffs;
        }

        // Get polynomial parameters (initial, final).
        virtual VectorX points_all() const
        {
            VectorX points(D * 4);
            points << _p0, _v0, _p1, _v1;
            return points;
        }

        // Get initial polynomial parameters (initial, final).
        virtual VectorX points_initial() const
        {
            VectorX points(D * 2);
            points << _p0, _v0;
            return points;
        }

        // Get final polynomial parameters (initial, final).
        virtual VectorX points_target() const
        {
            VectorX points(D * 2);
            points << _p1, _v1;
            return points;
        }

        // Set polynomial initial and final positions and velocities.
        virtual void set_points(const VectorX& x)
        {
            // assume x.size() == D*4
            _p0 = x.head(D); // initial position
            _v0 = x.segment(D, D); // initial velocity
            _p1 = x.segment(2 * D, D); // final position
            _v1 = x.tail(D); // final velocity

            _c0 = _p0;
            _c1 = _v0;
            _c2 = 3. * _p1 - 3. * _p0 - 2. * _v0 - _v1;
            _c3 = -2. * _p1 + 2. * _p0 + _v0 + _v1;
        }

        // Set polynomial parameters manually.
        virtual void set_coeffs(const VectorX& x)
        {
            // assume x.size() == D*4
            _c0 = x.head(D);
            _c1 = x.segment(D, D);
            _c2 = x.segment(2 * D, D);
            _c3 = x.tail(D);

            _p0 = position(0.);
            _v0 = velocity(0.);
            _p1 = position(1.);
            _v1 = velocity(1.);
        }

        // get position jacobian
        virtual MatrixX jac_pos(double t) const
        {
            const double t2 = t * t;
            const double t3 = t * t2;
            MatrixX jac = MatrixX::Zero(D, D * 4);

            for (unsigned int i = 0; i < D; i++) {
                // initial position
                jac(i, i) = 1. - 3. * t2 + 2. * t3;
                // initial velocity
                jac(i, D + i) = 1. * t - 2. * t2 + 1. * t3;
                // final position
                jac(i, 2 * D + i) = 3. * t2 - 2. * t3;
                // final velocity
                jac(i, 3 * D + i) = -1. * t2 + 1. * t3;
            }

            return jac;
        }

        // get velocity jacobian
        virtual MatrixX jac_vel(double t) const
        {
            const double t2 = t * t;
            MatrixX jac = MatrixX::Zero(D, D * 4);

            for (unsigned int i = 0; i < D; i++) {
                // initial position
                jac(i, i) = -6. * t + 6. * t2;
                // initial velocity
                jac(i, D + i) = 1. - 4. * t + 3. * t2;
                // final position
                jac(i, 2 * D + i) = 6. * t - 6. * t2;
                // final velocity
                jac(i, 3 * D + i) = -2. * t + 3. * t2;
            }

            return jac;
        }

        // get acceleration jacobian
        virtual MatrixX jac_acc(double t) const
        {
            MatrixX jac = MatrixX::Zero(D, D * 4);

            for (unsigned int i = 0; i < D; i++) {
                // initial position
                jac(i, i) = -6. + 12. * t;
                // initial velocity
                jac(i, D + i) = -4. + 6. * t;
                // final position
                jac(i, 2 * D + i) = 6. - 12. * t;
                // final velocity
                jac(i, 3 * D + i) = -2. + 6. * t;
            }

            return jac;
        }

    protected:
        VectorD _c0, _c1, _c2, _c3; // polynomial coefficients
        VectorD _p0, _v0, _p1, _v1; // points

        CubicHermitePolynomial() {}
    };

    using CubicHermitePolynomialReg2D = CubicHermitePolynomial<2>;
    using CubicHermitePolynomialReg3D = CubicHermitePolynomial<3>;
} // namespace cspl