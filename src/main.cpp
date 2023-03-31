#include <Eigen/Dense>

#include <iostream>

class PolynomialPair {
public:
    PolynomialPair(const Eigen::Vector2d& pos0, const Eigen::Vector2d& posFin, const Eigen::Vector2d& vel0, const Eigen::Vector2d& velFin, double t0, double tFin)
    // : _startPos(pos0), _finalPos(posFin), _startVel(vel0), _finalVel(velFin), _startTime(t0), _finalTime(tFin)
    {
        _a0 = pos0;
        _a1 = vel0;
        _a2 = (3 / std::pow(tFin, 2)) * (posFin - pos0) - (2 / tFin) * vel0 - (1 / tFin) * velFin;
        _a3 = -(2 / std::pow(tFin, 3)) * (posFin - pos0) + (1 / std::pow(tFin, 2)) * (velFin - vel0);
    };

    Eigen::Vector2d position(double t)
    {
        return _a0 + (_a1 * t) + (_a2 * std::pow(t, 2)) + (_a3 * std::pow(t, 3));
    };

    Eigen::Vector2d velocity(double t)
    {
        return _a1 + (2 * _a2 * t) + (3 * _a3 * std::pow(t, 2));
    };

private:
    // Eigen::Vector2d _startPos, _finalPos, _startVel, _finalVel;
    // double _startTime, _finalTime;
    Eigen::Vector2d _a0, _a1, _a2, _a3;
};

int main()
{
    Eigen::Vector2d x0 = {0., 0.};
    Eigen::Vector2d v0 = {1., 0.};

    Eigen::Vector2d x1 = {1., 1.};
    Eigen::Vector2d v1 = {1., 0.};

    auto pair = new PolynomialPair(x0, x1, v0, v1, 0., 1.);

    std::cout << "Positions:" << std::endl;
    for (double t = 0.; t <= 1.; t += 0.1) {
        std::cout << pair->position(t).transpose() << std::endl;
    }

    std::cout << "Velocities:" << std::endl;
    for (double t = 0.; t <= 1.; t += 0.1) {
        std::cout << pair->velocity(t).transpose() << std::endl;
    }

    return 0;
}