#include <Eigen/Dense>

namespace cspl {
    template <unsigned int D>
    class CubicHermitePolynomial {
    public:
        using Vec = Eigen::Matrix<double, D, 1>;

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

    protected:
        Vec _a0, _a1, _a2, _a3;
    };

    using CubicHermitePolynomial2D = CubicHermitePolynomial<2>;
    using CubicHermitePolynomial3D = CubicHermitePolynomial<3>;
} // namespace cspl
