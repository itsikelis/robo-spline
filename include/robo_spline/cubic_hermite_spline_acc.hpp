#pragma once

#include <robo_spline/cubic_hermite_spline.hpp>

namespace rspl {
    /**
     * @brief Acceleration Cubic Hermite Spline class.
     *
     * @tparam D The dimensionality of the trajectory.
     */
    template <unsigned int D>
    class CubicHermiteSplineAcc : public CubicHermiteSpline<D> {
    public:
        using VecD = typename CubicHermiteSpline<D>::VecD;
        using Vector = typename CubicHermiteSpline<D>::Vector;

        /**
         * @brief Constructs a cubic Hermite spline object given initial position, velocity and acceleration and the final position.
         *
         * @param p0 The initial position vector.
         * @param v0 The initial velocity vector.
         * @param a0 The initial acceleration vector.
         * @param p1 The final position vector.
         * @param T The duration of the polynomial.
         */
        CubicHermiteSplineAcc(const VecD& p0, const VecD& v0, const VecD& a0, const VecD& p1, double T = 1.) : CubicHermiteSpline<D>()
        {
            Vector x(D * 4);
            x << p0, v0, a0, p1;
            this->calc_coeffs_from_points(x, T);
        }

        /**
         * @brief Set spline initial and final points and their derivatives.
         *
         * @param x Parameters vector (initial and final positions, initial velocity and acceleration).
         * @param T The duration of the polynomial.
         */
        void calc_coeffs_from_points(const Vector& x, double T = 1.) override
        {
            // assume x.size() == D*4
            this->_p0 = x.head(D);
            this->_v0 = x.segment(D, D);
            this->_v1 = x.segment(2 * D, D); // _v1 is treated as a0
            this->_p1 = x.tail(D);

            const double o_T = 1. / T;
            const double o_T2 = 1. / (T * T);
            const double o_T3 = 1. / (T * T * T);

            this->_c0 = this->_p0;
            this->_c1 = this->_v0;
            this->_c2 = this->_v1 / 2.; // _v1 is treated as a0
            this->_c3 = this->_p1 * o_T3 - this->_c2 * o_T - this->_c1 * o_T2 - this->_c0 * o_T3;
            this->_T = T;
        }

        /**
         * @brief Get the position derivative of the spline at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the position derivative at the given time.
         */
        Vector deriv_pos(double t) const override
        {
            const double t2 = t * t;
            const double t3 = t * t2;
            Vector deriv = Vector::Zero(4);

            const double o_T = 1. / this->_T;
            const double o_T2 = 1. / (this->_T * this->_T);
            const double o_T3 = 1. / (this->_T * this->_T * this->_T);

            // initial position derivative
            deriv[0] = 1. - t3 * o_T3;
            // initial velocity derivative
            deriv[1] = t - t3 * o_T2;
            // initial acceleration derivative
            deriv[2] = 0.5 * t2 - 0.5 * t3 * o_T;
            // final position derivative
            deriv[3] = t3 * o_T3;

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

            const double o_T = 1. / this->_T;
            const double o_T2 = 1. / (this->_T * this->_T);
            const double o_T3 = 1. / (this->_T * this->_T * this->_T);

            // initial position derivative
            deriv[0] = -3. * t2 * o_T3;
            // initial velocity derivative
            deriv[1] = 1. - 3. * t2 * o_T2;
            // initial acceleration derivative
            deriv[2] = t - 1.5 * t2 * o_T;
            // final position derivative
            deriv[3] = 3. * t2 * o_T3;

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

            const double o_T = 1. / this->_T;
            const double o_T2 = 1. / (this->_T * this->_T);
            const double o_T3 = 1. / (this->_T * this->_T * this->_T);

            // initial position derivative
            deriv[0] = -6. * t * o_T3;
            // initial velocity derivative
            deriv[1] = -6. * t * o_T2;
            // initial acceleration derivative
            deriv[2] = 1. - 3. * t * o_T;
            // final position derivative
            deriv[3] = 6. * t * o_T3;

            return deriv;
        }

        /**
         * @brief Returns a Vector object representing all parameters (initial and final) of the cubic Hermite spline.
         *
         * @return A Vector object containing all parameters (initial and final) of the cubic Hermite spline.
         */
        Vector points_all() const override
        {
            Vector points(D * 4);
            points << this->_p0, this->_v0, this->_v1, this->_p1; // _v1 is treated as a0
            return points;
        }

        /**
         * @brief Returns a Vector object representing the initial parameters (position, velocity, acceleration) of the cubic Hermite spline.
         *
         * @return A Vector object containing the initial parameters (position, velocity, acceleration) of the cubic Hermite spline.
         */
        Vector points_initial() const override
        {
            Vector points(D * 3);
            points << this->_p0, this->_v0, this->_v1; // _v1 is treated as a0
            return points;
        }

        /**
         * @brief Returns a Vector object representing the final parameters (position, velocity, acceleration) of the cubic Hermite spline.
         *
         * @return A Vector object containing the final parameters (position, velocity, acceleration) of the cubic Hermite spline.
         */
        Vector points_target() const override
        {
            Vector points(D);
            points << this->_p1;
            return points;
        }

        /**
         * @brief Set spline coefficients manually.
         *
         * @param x Coefficients vector.
         * @param T The duration of the polynomial.
         */
        void set_coeffs(const Vector& x, double T = 1.) override
        {
            // assume x.size() == D*4
            this->_c0 = x.head(D);
            this->_c1 = x.segment(D, D);
            this->_c2 = x.segment(2 * D, D);
            this->_c3 = x.tail(D);

            this->_p0 = this->position(0.);
            this->_v0 = this->velocity(0.);
            this->_v1 = this->acceleration(0.); // _v1 is treated as a0
            this->_p1 = this->position(T);

            this->_T = T;
        }
    };

    using CubicHermiteSplineAcc2D = CubicHermiteSplineAcc<2>;
    using CubicHermiteSplineAcc3D = CubicHermiteSplineAcc<3>;
} // namespace rspl
