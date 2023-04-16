#pragma once

#include <Eigen/Dense>

namespace cspl {
    // Interface for the polynomial class.
    template <unsigned int D>
    class CubicHermitePolynomial {
    public:
        using VecD = Eigen::Matrix<double, D, 1>; // D dimensional Vector.
        using Vector = Eigen::Matrix<double, -1, 1>; // X dimensional Vector.

        /**
         * Constructs a cubic Hermite polynomial object using the given initial and final position and velocity vectors.
         *
         * @tparam D The dimension of the vector space.
         * @param p0 The initial position vector.
         * @param v0 The initial velocity vector.
         * @param p1 The final position vector.
         * @param v1 The final velocity vector.
         */
        CubicHermitePolynomial(const VecD& p0, const VecD& v0, const VecD& p1, const VecD& v1)
        {
            Vector p(D * 4);
            p << p0, v0, p1, v1;
            set_points(p);
        }

        /**
         * Returns the position of a cubic Hermite polynomial at a normalized time in the range [0,1].
         *
         * @param t The normalized time value in the range [0,1].
         * @return A VecD object representing the position of the cubic Hermite polynomial at the given time.
         */
        VecD position(double t) const
        {
            return _c0 + (_c1 * t) + (_c2 * t * t) + (_c3 * t * t * t);
        }

        /**
         * Returns the velocity of a cubic Hermite polynomial at a normalized time in the range [0,1].
         *
         * @param t The normalized time value in the range [0,1].
         * @return A VecD object representing the velocity of the cubic Hermite polynomial at the given time.
         */
        VecD velocity(double t) const
        {
            return _c1 + (2 * _c2 * t) + (3 * _c3 * t * t);
        }

        /**
         * Returns the acceleration of a cubic Hermite polynomial at a normalized time in the range [0,1].
         *
         * @param t The normalized time value in the range [0,1].
         * @return A VecD object representing the acceleration of the cubic Hermite polynomial at the given time.
         */
        VecD acceleration(double t) const
        {
            return (2 * _c2) + (6 * _c3 * t);
        }

        /**
         * Returns a Vector object representing the coefficients of the cubic Hermite polynomial.
         *
         * @return A Vector object containing the coefficients of the cubic Hermite polynomial.
         */
        Vector coeffs() const
        {
            Vector coeffs(D * 4);
            coeffs << _c0, _c1, _c2, _c3;
            return coeffs;
        }

        /**
         * Returns a Vector object representing all parameters (initial and final) of the cubic Hermite polynomial.
         *
         * @return A Vector object containing all parameters (initial and final) of the cubic Hermite polynomial.
         */
        virtual Vector points_all() const
        {
            Vector points(D * 4);
            points << _p0, _v0, _p1, _v1;
            return points;
        }

        /**
         * Returns a Vector object representing the initial parameters (position and velocity) of the cubic Hermite polynomial.
         *
         * @return A Vector object containing the initial parameters (position and velocity) of the cubic Hermite polynomial.
         */
        virtual Vector points_initial() const
        {
            Vector points(D * 2);
            points << _p0, _v0;
            return points;
        }

        /**
         * Returns a Vector object representing the final parameters (position and velocity) of the cubic Hermite polynomial.
         *
         * @return A Vector object containing the final parameters (position and velocity) of the cubic Hermite polynomial.
         */
        virtual Vector points_target() const
        {
            Vector points(D * 2);
            points << _p1, _v1;
            return points;
        }

        /**
         * Set polynomial initial and final positions and velocities.
         *
         * @param x Parameters vector (initial and final positions and velocities).
         */
        virtual void set_points(const Vector& x)
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

        /**
         * Set polynomial coefficients manually.
         *
         * @param x Coefficients vector.
         */
        virtual void set_coeffs(const Vector& x)
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

        /**
         * Get the position derivative of the polynomial at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the position derivative at the given time.
         */
        virtual Vector deriv_pos(double t) const
        {
            const double t2 = t * t;
            const double t3 = t * t2;
            Vector deriv = Vector::Zero(4);

            // initial position derivative
            deriv[0] = 1. - 3. * t2 + 2. * t3;
            // initial velocity derivative
            deriv[1] = 1. * t - 2. * t2 + 1. * t3;
            // final position derivative
            deriv[2] = 3. * t2 - 2. * t3;
            // final velocity derivative
            deriv[3] = -1. * t2 + 1. * t3;

            return deriv;
        }

        /**
         * Get the velocity derivative of the polynomial at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the velocity derivative at the given time.
         */
        virtual Vector deriv_vel(double t) const
        {
            const double t2 = t * t;
            Vector deriv = Vector::Zero(4);

            // initial position
            deriv[0] = -6. * t + 6. * t2;
            // initial velocity
            deriv[1] = 1. - 4. * t + 3. * t2;
            // final position
            deriv[2] = 6. * t - 6. * t2;
            // final velocity
            deriv[3] = -2. * t + 3. * t2;

            return deriv;
        }

        /**
         * Get the acceleration derivative of the polynomial at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the acceleration derivative at the given time.
         */
        virtual Vector deriv_acc(double t) const
        {
            Vector deriv = Vector::Zero(4);

            // initial position
            deriv[0] = -6. + 12. * t;
            // initial velocity
            deriv[1] = -4. + 6. * t;
            // final position
            deriv[2] = 6. - 12. * t;
            // final velocity
            deriv[3] = -2. + 6. * t;

            return deriv;
        }

    protected:
        VecD _c0, _c1, _c2, _c3; // Polynomial Coefficients.
        VecD _p0, _v0, _p1, _v1; // Points.

        CubicHermitePolynomial() {}
    };

    using CubicHermitePolynomialReg2D = CubicHermitePolynomial<2>;
    using CubicHermitePolynomialReg3D = CubicHermitePolynomial<3>;
} // namespace cspl