#pragma once

#include <cspl/i_cubic_hermite_polynomial.hpp>

namespace cspl {
    template <unsigned int D>
    class CubicHermitePolynomialAcc : ICubicHermitePolynomial<D, CubicHermitePolynomialAcc<D>> {
    public:
        using VectorD = typename ICubicHermitePolynomial<D, CubicHermitePolynomialAcc<D>>::VectorD;
        using VectorX = typename ICubicHermitePolynomial<D, CubicHermitePolynomialAcc<D>>::VectorX;
        using MatrixX = typename ICubicHermitePolynomial<D, CubicHermitePolynomialAcc<D>>::MatrixX;

        // p0, v0, a0 initial position, velocity and acceleration, p1: final position
        CubicHermitePolynomialAcc(const VectorD& p0, const VectorD& v0, const VectorD& a0, const VectorD& p1)
        {
            VectorX p(D * 4);
            p << p0, v0, a0, p1;
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

        // Get initial polynomial parameters (initial, final).
        VectorX points_initial() const
        {
            VectorX points(D * 3);
            points << _p0, _v0, _a0;
            return points;
        }

        // Get final polynomial parameters (initial, final).
        VectorX points_target() const
        {
            VectorX points(D);
            points << _p1;
            return points;
        }

        // Get polynomial parameters (initial, final).
        VectorX points_all() const
        {
            VectorX points(D * 4);
            points << _p0, _v0, _a0, _p1;
            return points;
        }

        // Set polynomial parameters manually.
        void set_coeffs(const VectorX& x)
        {
            // assume x.size() == D*4
            _c0 = x.head(D);
            _c1 = x.segment(D, D);
            _c2 = x.segment(2 * D, D);
            _c3 = x.tail(D);

            _p0 = position(0.);
            _v0 = velocity(0.);
            _a0 = acceleration(1.);
            _p1 = position(1.);
        }

        // Set polynomial initial position, velocity, acceleration and final position.
        void set_points(const VectorX& x)
        {
            // assume x.size() == D*4
            _p0 = x.head(D);
            _v0 = x.segment(D, D);
            _a0 = x.segment(2 * D, D);
            _p1 = x.tail(D);

            _c0 = _p0;
            _c1 = _v0;
            _c2 = _a0 / 2.;
            _c3 = _p1 - _c2 - _c1 - _c0;
        }

        // TODO : Jacobians.
        MatrixX jac_pos(double t) const
        {
            const double t2 = t * t;
            const double t3 = t * t2;
            MatrixX jac = MatrixX::Zero(D, D * 4);

            for (unsigned int i = 0; i < D; i++) {
                // initial position
                jac(i, i) = 1. - 1. * t3;
                // initial velocity
                jac(i, D + i) = 1. * t - 1. * t3;
                // initial acceleration
                jac(i, 2 * D + i) = 0.5 * t2 - 0.5 * t3;
                // final position
                jac(i, 3 * D + i) = 1. * t3;
            }

            return jac;
        }

        MatrixX jac_vel(double t) const
        {
            const double t2 = t * t;
            MatrixX jac = MatrixX::Zero(D, D * 4);

            for (unsigned int i = 0; i < D; i++) {
                // initial position
                jac(i, i) = -3. * t2;
                // initial velocity
                jac(i, D + i) = 1. - 3. * t2;
                // initial acceleration
                jac(i, 2 * D + i) = t - 1.5 * t2;
                // final position
                jac(i, 3 * D + i) = 3. * t2;
            }

            return jac;
        }

        MatrixX jac_acc(double t) const
        {
            MatrixX jac = MatrixX::Zero(D, D * 4);

            for (unsigned int i = 0; i < D; i++) {
                // initial position
                jac(i, i) = -6. * t;
                // initial velocity
                jac(i, D + i) = -6. * t;
                // initial acceleration
                jac(i, 2 * D + i) = 1. - 3. * t;
                // final position
                jac(i, 3 * D + i) = 6. * t;
            }

            return jac;
        }

    protected:
        VectorD _c0, _c1, _c2, _c3; // polynomial coefficients
        VectorD _p0, _v0, _a0, _p1; // initial and final positions and velocities
    };

    using CubicHermitePolynomialAcc2D = CubicHermitePolynomialAcc<2>;
    using CubicHermitePolynomialAcc3D = CubicHermitePolynomialAcc<3>;
} // namespace cspl
