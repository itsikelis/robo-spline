#pragma once

#include <cspl/cubic_hermite_polynomial.hpp>

namespace cspl {
    /**
     * @brief Acceleration Cubic Hermite Polynomial class.
     *
     * @tparam D The dimensionality of the trajectory.
     */
    template <unsigned int D>
    class CubicHermitePolynomialAcc : public CubicHermitePolynomial<D> {
    public:
        using VecD = typename CubicHermitePolynomial<D>::VecD;
        using Vector = typename CubicHermitePolynomial<D>::Vector;

        /**
         * @brief Constructs a cubic Hermite polynomial object given initial position, velocity and acceleration and the final position.
         *
         * @param p0 The initial position vector.
         * @param v0 The initial velocity vector.
         * @param a0 The initial acceleration vector.
         * @param p1 The final position vector.
         */
        CubicHermitePolynomialAcc(const VecD& p0, const VecD& v0, const VecD& a0, const VecD& p1) : CubicHermitePolynomial<D>()
        {
            Vector p(D * 4);
            p << p0, v0, a0, p1;
            this->set_points(p);
        }

        /**
         * @brief Returns a Vector object representing all parameters (initial and final) of the cubic Hermite polynomial.
         *
         * @return A Vector object containing all parameters (initial and final) of the cubic Hermite polynomial.
         */
        Vector points_all() const override
        {
            Vector points(D * 4);
            points << this->_p0, this->_v0, this->_v1, this->_p1; // _v1 is treated as a0
            return points;
        }

        /**
         * @brief Returns a Vector object representing the initial parameters (position, velocity, acceleration) of the cubic Hermite polynomial.
         *
         * @return A Vector object containing the initial parameters (position, velocity, acceleration) of the cubic Hermite polynomial.
         */
        Vector points_initial() const override
        {
            Vector points(D * 3);
            points << this->_p0, this->_v0, this->_v1; // _v1 is treated as a0
            return points;
        }

        /**
         * @brief Returns a Vector object representing the final parameters (position, velocity, acceleration) of the cubic Hermite polynomial.
         *
         * @return A Vector object containing the final parameters (position, velocity, acceleration) of the cubic Hermite polynomial.
         */
        Vector points_target() const override
        {
            Vector points(D);
            points << this->_p1;
            return points;
        }

        /**
         * @brief Set polynomial initial and final positions and initial velocity and acceleration.
         *
         * @param x Parameters vector (initial and final positions, initial velocity and acceleration).
         */
        void set_points(const Vector& x) override
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

        /**
         * @brief Set polynomial coefficients manually.
         *
         * @param x Coefficients vector.
         */
        void set_coeffs(const Vector& x) override
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

        /**
         * @brief Get the position derivative of the polynomial at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the position derivative at the given time.
         */
        Vector deriv_pos(double t) const override
        {
            const double t2 = t * t;
            const double t3 = t * t2;
            Vector deriv = Vector::Zero(4);

            // initial position derivative
            deriv[0] = 1. - 1. * t3;
            // initial velocity derivative
            deriv[1] = 1. * t - 1. * t3;
            // initial acceleration derivative
            deriv[2] = 0.5 * t2 - 0.5 * t3;
            // final position derivative
            deriv[3] = 1. * t3;

            return deriv;
        }

        /**
         * @brief Get the velocity derivative of the polynomial at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the velocity derivative at the given time.
         */
        Vector deriv_vel(double t) const override
        {
            const double t2 = t * t;
            Vector deriv = Vector::Zero(4);

            // initial position derivative
            deriv[0] = -3. * t2;
            // initial velocity derivative
            deriv[1] = 1. - 3. * t2;
            // initial acceleration derivative
            deriv[2] = t - 1.5 * t2;
            // final position derivative
            deriv[3] = 3. * t2;

            return deriv;
        }

        /**
         * @brief Get the acceleration derivative of the polynomial at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the acceleration derivative at the given time.
         */
        Vector deriv_acc(double t) const override
        {
            Vector deriv = Vector::Zero(4);

            // initial position derivative
            deriv[0] = -6. * t;
            // initial velocity derivative
            deriv[1] = -6. * t;
            // initial acceleration derivative
            deriv[2] = 1. - 3. * t;
            // final position derivative
            deriv[3] = 6. * t;

            return deriv;
        }
    };

    using CubicHermitePolynomialAcc2D = CubicHermitePolynomialAcc<2>;
    using CubicHermitePolynomialAcc3D = CubicHermitePolynomialAcc<3>;
} // namespace cspl
