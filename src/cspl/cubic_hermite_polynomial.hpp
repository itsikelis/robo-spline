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
            if (regular) {
                // Regular constructor with initial and final positions and velocities
                Vec final_position = p1; // This is for clarity, we could avoid copies here
                Vec final_velocity = p2;

                _a0 = initial_position;
                _a1 = initial_velocity;
                _a2 = 3. * final_position - 3. * initial_position - 2. * initial_velocity - final_velocity;
                _a3 = -2. * final_position + 2. * initial_position + initial_velocity + final_velocity;
            }
            else {
                // Constructor with initial position, velocity, acceleration and final position
                Vec initial_acceleration = p1; // This is for clarity, we could avoid copies here
                Vec final_position = p2;

                _a0 = initial_position;
                _a1 = initial_velocity;
                _a2 = initial_acceleration / 2.;
                _a3 = final_position - _a2 - _a1 - _a0;
            }

            _p0 = initial_position;
            _v0 = initial_velocity;
            _p1 = p1;
            _v1 = p2;
        }

        Vec position(double t) const
        {
            return _a0 + (_a1 * t) + (_a2 * t * t) + (_a3 * t * t * t);
        }

        Vec velocity(double t) const
        {
            return _a1 + (2 * _a2 * t) + (3 * _a3 * t * t);
        }

        Vec acceleration(double t) const
        {
            return (2 * _a2) + (6 * _a3 * t);
        }

        Vector params_nodes() const
        {
            Vector params(D * 4);
            params << _p0, _v0, _p1, _v1;

            return params;
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

    protected:
        Vec _a0, _a1, _a2, _a3;
        Vec _p0, _v0, _p1, _v1;
    };

    using CubicHermitePolynomial2D = CubicHermitePolynomial<2>;
    using CubicHermitePolynomial3D = CubicHermitePolynomial<3>;
} // namespace cspl
