#pragma once

#include <cspl/cubic_hermite_polynomial.hpp>

namespace cspl {
    template <unsigned int D>
    class CubicHermitePolynomialAcc : public CubicHermitePolynomial<D> {
    public:
        using VectorD = typename CubicHermitePolynomial<D>::VectorD;
        using VectorX = typename CubicHermitePolynomial<D>::VectorX;
        using MatrixX = typename CubicHermitePolynomial<D>::MatrixX;

        // p0, v0, a0 initial position, velocity and acceleration, p1: final position
        CubicHermitePolynomialAcc(const VectorD& p0, const VectorD& v0, const VectorD& a0, const VectorD& p1) : CubicHermitePolynomial<D>()
        {
            VectorX p(D * 4);
            p << p0, v0, a0, p1;
            this->set_points(p);
        }

        // Get initial polynomial parameters (initial, final).
        VectorX points_initial() const override
        {
            VectorX points(D * 3);
            points << this->_p0, this->_v0, this->_v1; // _v1 is treated as a0
            return points;
        }

        // Get final polynomial parameters (initial, final).
        VectorX points_target() const override
        {
            VectorX points(D);
            points << this->_p1;
            return points;
        }

        // Get polynomial parameters (initial, final).
        VectorX points_all() const override
        {
            VectorX points(D * 4);
            points << this->_p0, this->_v0, this->_v1, this->_p1; // _v1 is treated as a0
            return points;
        }

        // Set polynomial parameters manually.
        void set_coeffs(const VectorX& x) override
        {
            // assume x.size() == D*4
            this->_c0 = x.head(D);
            this->_c1 = x.segment(D, D);
            this->_c2 = x.segment(2 * D, D);
            this->_c3 = x.tail(D);

            this->_p0 = this->position(0.);
            this->_v0 = this->velocity(0.);
            this->_v1 = this->acceleration(1.); // _v1 is treated as a0
            this->_p1 = this->position(1.);
        }

        // Set polynomial initial position, velocity, acceleration and final position.
        void set_points(const VectorX& x) override
        {
            // assume x.size() == D*4
            this->_p0 = x.head(D);
            this->_v0 = x.segment(D, D);
            this->_v1 = x.segment(2 * D, D); // _v1 is treated as a0
            this->_p1 = x.tail(D);

            this->_c0 = this->_p0;
            this->_c1 = this->_v0;
            this->_c2 = this->_v1 / 2.; // _v1 is treated as a0
            this->_c3 = this->_p1 - this->_c2 - this->_c1 - this->_c0;
        }

        // TODO : Jacobians.
        MatrixX jac_pos(double t) const override
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

        MatrixX jac_vel(double t) const override
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

        MatrixX jac_acc(double t) const override
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
    };

    using CubicHermitePolynomialAcc2D = CubicHermitePolynomialAcc<2>;
    using CubicHermitePolynomialAcc3D = CubicHermitePolynomialAcc<3>;
} // namespace cspl
