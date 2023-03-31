#include <Eigen/Dense>

#include <iostream>

template <int D>
class CubicHermitePolynomial {
public:
    using Vec = Eigen::Matrix<double, D, 1>;

    CubicHermitePolynomial(const Vec& initial_position, const Vec& initial_velocity, const Vec& final_position, const Vec& final_velocity)
    {
        _a0 = initial_position;
        _a1 = initial_velocity;
        _a2 = 3. * final_position - 3. * initial_position - 2. * initial_velocity - final_velocity;
        _a3 = -2. * final_position + 2. * initial_position + initial_velocity + final_velocity;
    }

    Eigen::Vector2d position(double t)
    {
        return _a0 + (_a1 * t) + (_a2 * t * t) + (_a3 * t * t * t);
    };

    Eigen::Vector2d velocity(double t)
    {
        return _a1 + (2 * _a2 * t) + (3 * _a3 * t * t);
    };

protected:
    Vec _a0, _a1, _a2, _a3;
};

int main()
{
    Eigen::Vector2d x0 = {0., 0.};
    Eigen::Vector2d v0 = {1., 0.};

    Eigen::Vector2d x1 = {1., 1.};
    Eigen::Vector2d v1 = {1., 0.};

    auto pol = CubicHermitePolynomial<2>(x0, x1, v0, v1);

    std::cout << "Positions:" << std::endl;
    for (double t = 0.; t <= 1.; t += 0.1) {
        std::cout << pol.position(t).transpose() << std::endl;
    }

    std::cout << "Velocities:" << std::endl;
    for (double t = 0.; t <= 1.; t += 0.1) {
        std::cout << pol.velocity(t).transpose() << std::endl;
    }

    return 0;
}