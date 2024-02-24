#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace rspl {
    /**
     * @brief Regular Cubic Hermite Polynomial class.
     *
     * @tparam D The dimensionality of the trajectory.
     */
    template <unsigned int D>
    class CubicHermiteSpline {
    public:
        using VecD = Eigen::Matrix<double, D, 1>; // D dimensional Vector.
        using Vector = Eigen::Matrix<double, -1, 1>; // dynamic Vector.
        using Jacobian = Eigen::Matrix<double, -1, -1>; // dynamic Matrix.
        using SparseJacobian = Eigen::SparseMatrix<double, Eigen::RowMajor>; // dynamic Sparse Matrix.

        CubicHermiteSpline() = default;

        /**
         * @brief Constructs a cubic Hermite polynomial object using the given initial and final position and velocity vectors.
         *
         * @param p0 The initial position vector.
         * @param v0 The initial velocity vector.
         * @param p1 The final position vector.
         * @param v1 The final velocity vector.
         * @param T The duration of the polynomial.
         */
        CubicHermiteSpline(const VecD& p0, const VecD& v0, const VecD& p1, const VecD& v1, double T = 1.)
        {
            Vector x(D * 4);
            x << p0, v0, p1, v1;
            calc_coeffs_from_points(x, T);
        }

        /**
         * @brief Returns the duration of the polynomial.
         *
         * @return The duration of the polynomial.
         */
        double duration() const { return _T; }

        /**
         * @brief Set polynomial initial and final positions and velocities.
         *
         * @param x Parameters vector (initial and final positions and velocities).
         * @param T The duration of the polynomial.
         */
        virtual void calc_coeffs_from_points(const Vector& x, double T = 1.)
        {
            // assume x.size() == D*4
            _p0 = x.head(D); // initial position
            _v0 = x.segment(D, D); // initial velocity
            _p1 = x.segment(2 * D, D); // final position
            _v1 = x.tail(D); // final velocity

            const double o_T = 1. / T;
            const double o_T2 = 1. / (T * T);
            const double o_T3 = 1. / (T * T * T);

            _c0 = _p0;
            _c1 = _v0;
            _c2 = 3. * _p1 * o_T2 - 3. * _p0 * o_T2 - 2. * _v0 * o_T - _v1 * o_T;
            _c3 = -2. * _p1 * o_T3 + 2. * _p0 * o_T3 + _v0 * o_T2 + _v1 * o_T2;
            _T = T;
        }

        /**
         * @brief Returns the position of a cubic Hermite polynomial at a normalized time in the range [0,T].
         *
         * @param t The normalized time value in the range [0,T].
         * @return A VecD object representing the position of the cubic Hermite polynomial at the given time.
         */
        VecD position(double t) const
        {
            const double t2 = t * t;
            const double t3 = t * t2;

            return _c0 + (_c1 * t) + (_c2 * t2) + (_c3 * t3);
        }

        /**
         * @brief Returns the velocity of a cubic Hermite polynomial at a normalized time in the range [0,T].
         *
         * @param t The normalized time value in the range [0,T].
         * @return A VecD object representing the velocity of the cubic Hermite polynomial at the given time.
         */
        VecD velocity(double t) const
        {
            return _c1 + (2. * _c2 * t) + (3. * _c3 * t * t);
        }

        /**
         * @brief Returns the acceleration of a cubic Hermite polynomial at a normalized time in the range [0,T].
         *
         * @param t The normalized time value in the range [0,T].
         * @return A VecD object representing the acceleration of the cubic Hermite polynomial at the given time.
         */
        VecD acceleration(double t) const
        {
            return (2. * _c2) + (6. * _c3 * t);
        }

        /**
         * @brief Get the position derivative of the polynomial at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the partial derivatives at the given time.
         */
        virtual Vector deriv_pos(double t) const
        {
            const double t2 = t * t;
            const double t3 = t * t2;
            Vector deriv = Vector::Zero(4);

            const double o_T = 1. / _T;
            const double o_T2 = 1. / (_T * _T);
            const double o_T3 = 1. / (_T * _T * _T);

            // derivative w.r.t x0
            deriv[0] = 1. - 3. * t2 * o_T2 + 2. * t3 * o_T3;
            // derivative w.r.t v0
            deriv[1] = t - 2. * t2 * o_T + t3 * o_T2;
            // derivative w.r.t x1
            deriv[2] = 3. * t2 * o_T2 - 2. * t3 * o_T3;
            // derivative w.r.t v1
            deriv[3] = -t2 * o_T + t3 * o_T2;

            return deriv;
        }

        /**
         * @brief Get the position jacobian of the polynomial at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A matrix containing the jacobian at the given time.
         */
        virtual Jacobian jacobian_dense_pos(double t) const
        {
            Vector deriv = deriv_pos(t);

            Jacobian jac = Jacobian::Zero(D, D * 4);

            for (unsigned int i = 0; i < D; i++) {
                for (unsigned int j = 0; j < 4; j++) {
                    jac(i, i + j * D) = deriv[j];
                }
            }

            return jac;
        }

        /**
         * @brief Get the position jacobian of the polynomial at time t. (SparseMatrix version)
         *
         * @param t The time to evaluate the derivative at.
         * @return A matrix containing the jacobian at the given time.
         */
        virtual SparseJacobian jacobian_pos(double t) const
        {
            Vector deriv = deriv_pos(t);

            SparseJacobian jac(D, D * 4);

            for (unsigned int i = 0; i < D; i++) {
                for (unsigned int j = 0; j < 4; j++) {
                    jac.coeffRef(i, i + j * D) = deriv[j];
                }
            }

            return jac;
        }

        /**
         * @brief Get the velocity derivative of the polynomial at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the velocity derivative at the given time.
         */
        virtual Vector deriv_vel(double t) const
        {
            const double t2 = t * t;
            Vector deriv = Vector::Zero(4);

            const double o_T = 1. / _T;
            const double o_T2 = 1. / (_T * _T);
            const double o_T3 = 1. / (_T * _T * _T);

            // initial position
            deriv[0] = -6. * t * o_T2 + 6. * t2 * o_T3;
            // initial velocity
            deriv[1] = 1. - 4. * t * o_T + 3. * t2 * o_T2;
            // final position
            deriv[2] = 6. * t * o_T2 - 6. * t2 * o_T3;
            // final velocity
            deriv[3] = -2. * t * o_T + 3. * t2 * o_T2;

            return deriv;
        }

        /**
         * @brief Get the velocity jacobian of the polynomial at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A matrix containing the jacobian at the given time.
         */
        virtual Jacobian jacobian_dense_vel(double t) const
        {
            Vector deriv = deriv_vel(t);

            Jacobian jac = Jacobian::Zero(D, D * 4);

            for (unsigned int i = 0; i < D; i++) {
                for (unsigned int j = 0; j < 4; j++) {
                    jac(i, i + j * D) = deriv[j];
                }
            }

            return jac;
        }

        /**
         * @brief Get the velocity jacobian of the polynomial at time t. (SparseMatrix version)
         *
         * @param t The time to evaluate the derivative at.
         * @return A matrix containing the jacobian at the given time.
         */
        virtual SparseJacobian jacobian_vel(double t) const
        {
            Vector deriv = deriv_vel(t);

            SparseJacobian jac(D, D * 4);

            for (unsigned int i = 0; i < D; i++) {
                for (unsigned int j = 0; j < 4; j++) {
                    jac.coeffRef(i, i + j * D) = deriv[j];
                }
            }

            return jac;
        }

        /**
         * @brief Get the acceleration derivative of the polynomial at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the acceleration derivative at the given time.
         */
        virtual Vector deriv_acc(double t) const
        {
            Vector deriv = Vector::Zero(4);

            const double o_T = 1. / _T;
            const double o_T2 = 1. / (_T * _T);
            const double o_T3 = 1. / (_T * _T * _T);

            // initial position
            deriv[0] = -6. * o_T2 + 12. * t * o_T3;
            // initial velocity
            deriv[1] = -4. * o_T + 6. * t * o_T2;
            // final position
            deriv[2] = 6. * o_T2 - 12. * t * o_T3;
            // final velocity
            deriv[3] = -2. * o_T + 6. * t * o_T2;

            return deriv;
        }

        /**
         * @brief Get the acceleration jacobian of the polynomial at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A matrix containing the jacobian at the given time.
         */
        virtual Jacobian jacobian_dense_acc(double t) const
        {
            Vector deriv = deriv_acc(t);

            Jacobian jac = Jacobian::Zero(D, D * 4);

            for (unsigned int i = 0; i < D; i++) {
                for (unsigned int j = 0; j < 4; j++) {
                    jac(i, i + j * D) = deriv[j];
                }
            }

            return jac;
        }

        /**
         * @brief Get the acceleration jacobian of the polynomial at time t. (SparseMatrix version)
         *
         * @param t The time to evaluate the derivative at.
         * @return A matrix containing the jacobian at the given time.
         */
        virtual SparseJacobian jacobian_acc(double t) const
        {
            Vector deriv = deriv_acc(t);

            SparseJacobian jac(D, D * 4);

            for (unsigned int i = 0; i < D; i++) {
                for (unsigned int j = 0; j < 4; j++) {
                    jac.coeffRef(i, i + j * D) = deriv[j];
                }
            }

            return jac;
        }

        /**
         * @brief Returns a Vector object representing the coefficients of the cubic Hermite polynomial.
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
         * @brief Returns a Vector object representing all parameters (initial and final) of the cubic Hermite polynomial.
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
         * @brief Returns a Vector object representing the initial parameters (position and velocity) of the cubic Hermite polynomial.
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
         * @brief Returns a Vector object representing the final parameters (position and velocity) of the cubic Hermite polynomial.
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
         * @brief Set polynomial coefficients manually.
         *
         * @param x Coefficients vector.
         * @param T The duration of the polynomial.
         */
        virtual void set_coeffs(const Vector& x, double T = 1.)
        {
            // assume x.size() == D*4
            _c0 = x.head(D);
            _c1 = x.segment(D, D);
            _c2 = x.segment(2 * D, D);
            _c3 = x.tail(D);

            _p0 = position(0.);
            _v0 = velocity(0.);
            _p1 = position(T);
            _v1 = velocity(T);

            _T = T;
        }

    protected:
        VecD _c0, _c1, _c2, _c3; // Polynomial Coefficients.
        VecD _p0, _v0, _p1, _v1; // Points.
        double _T = 1.;
    };

    using CubicHermitePolynomialReg2D = CubicHermiteSpline<2>;
    using CubicHermitePolynomialReg3D = CubicHermiteSpline<3>;
} // namespace rspl
