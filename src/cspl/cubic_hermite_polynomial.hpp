#include <Eigen/Dense>

namespace cspl {
    template <unsigned int D>
    class CubicHermitePolynomial {
    public:
        using Vec = Eigen::Matrix<double, D, 1>;
        using Vector = Eigen::Matrix<double, -1, 1>;
        using Matrix = Eigen::Matrix<double, -1, -1>;

        CubicHermitePolynomial(const Vec& initial_position, const Vec& initial_velocity, const Vec& p1, const Vec& p2, bool regular = true)
        {
            Vector p(D * 4);
            p << initial_position, initial_velocity, p1, p2;
            set_node_params(p, regular);
        }

        // Get position at normalised time in [0-1].
        Vec position(double t) const
        {
            return _a0 + (_a1 * t) + (_a2 * t * t) + (_a3 * t * t * t);
        }

        // Get velocity at normalised time in [0-1].
        Vec velocity(double t) const
        {
            return _a1 + (2 * _a2 * t) + (3 * _a3 * t * t);
        }

        // Get acceleration at normalised time in [0-1].
        Vec acceleration(double t) const
        {
            return (2 * _a2) + (6 * _a3 * t);
        }

        // Get polynomial parameters.
        Vector coeff_params() const
        {
            Vector params(D * 4);
            params << _a0, _a1, _a2, _a3;

            return params;
        }

        // Set polynomial parameters manually.
        void set_coeff_params(const Vector& x, bool regular = true)
        {
            // assume x.size() == D*4
            _a0 = x.head(D);
            _a1 = x.segment(D, D);
            _a2 = x.segment(2 * D, D);
            _a3 = x.tail(D);

            if (regular) {
                _p0 = position(0.);
                _v0 = velocity(0.);
                _p1 = position(0.);
                _v1 = velocity(1.);
            }
            else {
                _p0 = position(0.);
                _v0 = velocity(0.);
                _p1 = acceleration(0.);
                _v1 = position(1.);
            }
        }

        // Get intial/target points (positions and velocities).
        Vector node_params() const
        {
            Vector params(D * 4);
            params << _p0, _v0, _p1, _v1;

            return params;
        }

        void set_node_params(const Vector& x, bool regular = true)
        {
            // assume x.size() == D*4
            _p0 = x.head(D);
            _v0 = x.segment(D, D);
            _p1 = x.segment(2 * D, D);
            _v1 = x.tail(D);

            Vec initial_position = _p0; // This is for clarity, we could avoid copies here
            Vec initial_velocity = _v0;

            if (regular) {
                // Regular constructor with initial and final positions and velocities
                Vec final_position = _p1; // This is for clarity, we could avoid copies here
                Vec final_velocity = _v1;

                _a0 = initial_position;
                _a1 = initial_velocity;
                _a2 = 3. * final_position - 3. * initial_position - 2. * initial_velocity - final_velocity;
                _a3 = -2. * final_position + 2. * initial_position + initial_velocity + final_velocity;
            }
            else {
                // Constructor with initial position, velocity, acceleration and final position
                Vec initial_acceleration = _p1; // This is for clarity, we could avoid copies here
                Vec final_position = _v1;

                _a0 = initial_position;
                _a1 = initial_velocity;
                _a2 = initial_acceleration / 2.;
                _a3 = final_position - _a2 - _a1 - _a0;
            }
        }

        Matrix jac_pos(double t, bool regular = true) const
        {
            const double t2 = t * t;
            const double t3 = t * t2;
            Matrix jac = Matrix::Zero(D, D * 4);
            if (regular) {
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
            }
            else {
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
            }

            return jac;
        }

        Matrix jac_vel(double t, bool regular = true) const
        {
            const double t2 = t * t;
            const double t3 = t * t2;
            Matrix jac = Matrix::Zero(D, D * 4);
            if (regular) {
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
            }
            else {
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
            }

            return jac;
        }

        Matrix jac_acc(double t, bool regular = true) const
        {
            const double t2 = t * t;
            const double t3 = t * t2;
            Matrix jac = Matrix::Zero(D, D * 4);
            if (regular) {
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
            }
            else {
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
            }

            return jac;
        }

    protected:
        Vec _a0, _a1, _a2, _a3;
        Vec _p0, _v0, _p1, _v1;
    };

    using CubicHermitePolynomial2D = CubicHermitePolynomial<2>;
    using CubicHermitePolynomial3D = CubicHermitePolynomial<3>;
} // namespace cspl
