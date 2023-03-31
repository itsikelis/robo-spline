#include "cubic_hermite_polynomial.hpp"

#include <iostream>

using CubicHermitePolynomial2D = CubicHermitePolynomial<2>;
using CubicHermitePolynomial3D = CubicHermitePolynomial<3>;

int main()
{
    Eigen::Vector3d x0 = {0., 0., 0.};
    Eigen::Vector3d v0 = {1., 0., 1.};

    Eigen::Vector3d x1 = {1., 1., 0.};
    Eigen::Vector3d v1 = {1., 1., 1.};

    double duration = 5.;

    auto pol = CubicHermitePolynomial3D(x0, v0, x1, v1, duration);

    std::cout << "Positions:" << std::endl;
    for (double t = 0.; t <= duration; t += 0.1 * duration) {
        std::cout << pol.position(t).transpose() << std::endl;
    }

    std::cout << "Velocities:" << std::endl;
    for (double t = 0.; t <= duration; t += 0.1 * duration) {
        std::cout << pol.velocity(t).transpose() << std::endl;
    }

    std::cout << "Accelerations:" << std::endl;
    for (double t = 0.; t <= duration; t += 0.1 * duration) {
        std::cout << pol.acceleration(t).transpose() << std::endl;
    }

    return 0;
}