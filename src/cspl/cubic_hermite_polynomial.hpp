#include <Eigen/Dense>

namespace cspl {
    template <unsigned int DIM>
    class CubicHermitePolynomial {
    public:
        using Vec = Eigen::Matrix<double, DIM, 1>;

        CubicHermitePolynomial(const Vec& initial_position, const Vec& initial_velocity, const Vec& final_position, const Vec& final_velocity)
        {
            _a0 = initial_position;
            _a1 = initial_velocity;
            _a2 = 3. * final_position - 3. * initial_position - 2. * initial_velocity - final_velocity;
            _a3 = -2. * final_position + 2. * initial_position + initial_velocity + final_velocity;
        }

        Vec position(double t)
        {
            return _a0 + (_a1 * t) + (_a2 * t * t) + (_a3 * t * t * t);
        }

        Vec velocity(double t)
        {
            return _a1 + (2 * _a2 * t) + (3 * _a3 * t * t);
        }

        Vec acceleration(double t)
        {
            return (2 * _a2) + (6 * _a3 * t);
        }

    protected:
        Vec _a0, _a1, _a2, _a3;
    };

    using CubicHermitePolynomial2D = CubicHermitePolynomial<2>;
    using CubicHermitePolynomial3D = CubicHermitePolynomial<3>;
} // namespace cspl
