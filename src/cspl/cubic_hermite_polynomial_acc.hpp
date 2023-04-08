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
            this->_v1 = this->acceleration(0.); // _v1 is treated as a0
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

        // get position derivative
        VectorX deriv_pos(double t) const override
        {
            const double t2 = t * t;
            const double t3 = t * t2;
            VectorX deriv = VectorX::Zero(4);

            // initial position
            deriv[0] = 1. - 1. * t3;
            // initial velocity
            deriv[1] = 1. * t - 1. * t3;
            // initial acceleration
            deriv[2] = 0.5 * t2 - 0.5 * t3;
            // final position
            deriv[3] = 1. * t3;

            return deriv;
        }

        // get velocity derivative
        VectorX deriv_vel(double t) const override
        {
            const double t2 = t * t;
            VectorX deriv = VectorX::Zero(4);

            // initial position
            deriv[0] = -3. * t2;
            // initial velocity
            deriv[1] = 1. - 3. * t2;
            // initial acceleration
            deriv[2] = t - 1.5 * t2;
            // final position
            deriv[3] = 3. * t2;

            return deriv;
        }

        // get acceleration derivative
        VectorX deriv_acc(double t) const override
        {
            VectorX deriv = VectorX::Zero(4);

            // initial position
            deriv[0] = -6. * t;
            // initial velocity
            deriv[1] = -6. * t;
            // initial acceleration
            deriv[2] = 1. - 3. * t;
            // final position
            deriv[3] = 6. * t;

            return deriv;
        }
    };

    using CubicHermitePolynomialAcc2D = CubicHermitePolynomialAcc<2>;
    using CubicHermitePolynomialAcc3D = CubicHermitePolynomialAcc<3>;
} // namespace cspl
