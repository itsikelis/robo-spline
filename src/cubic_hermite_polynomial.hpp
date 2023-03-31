#include <Eigen/Dense>

template <unsigned int D>
class CubicHermitePolynomial {
public:
    using Vec = Eigen::Matrix<double, D, 1>;

    CubicHermitePolynomial(const Vec& initial_position, const Vec& initial_velocity, const Vec& final_position, const Vec& final_velocity)
    {
        _duration - 1.;

        _a0 = initial_position;
        _a1 = initial_velocity;
        _a2 = 3. * final_position - 3. * initial_position - 2. * initial_velocity - final_velocity;
        _a3 = -2. * final_position + 2. * initial_position + initial_velocity + final_velocity;
    }

    CubicHermitePolynomial(const Vec& initial_position, const Vec& initial_velocity, const Vec& final_position, const Vec& final_velocity, double duration)
        : _duration(duration)
    {
        _a0 = initial_position;
        _a1 = initial_velocity;
        _a2 = 3. * final_position - 3. * initial_position - 2. * initial_velocity - final_velocity;
        _a3 = -2. * final_position + 2. * initial_position + initial_velocity + final_velocity;
    }

    Vec position(double t)
    {
        double t_normalized = t / _duration;
        return _a0 + (_a1 * t_normalized) + (_a2 * t_normalized * t_normalized) + (_a3 * t_normalized * t_normalized * t_normalized);
    }

    Vec velocity(double t)
    {
        double t_normalized = t / _duration;
        return _a1 + (2 * _a2 * t_normalized) + (3 * _a3 * t_normalized * t_normalized);
    };

    Vec acceleration(double t)
    {
        double t_normalized = t / _duration;
        return (2 * _a2) + (6 * _a3 * t_normalized);
    };

protected:
    double _duration; // Total trajectory duration in seconds.
    Vec _a0, _a1, _a2, _a3;
};
