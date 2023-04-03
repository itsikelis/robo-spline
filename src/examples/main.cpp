#include <iostream>

#include <cspl/cubic_hermite_polynomial.hpp>
#include <cspl/trajectory.hpp>

int main()
{
    Eigen::Vector3d x0 = {0., 0., 0.};
    Eigen::Vector3d v0 = {1., 0., 1.};

    Eigen::Vector3d x1 = {1., 1., 0.};
    Eigen::Vector3d v1 = {1., 1., 1.};

    Eigen::Vector3d x2 = {2., 2., 0.};
    Eigen::Vector3d v2 = {2., 2., 1.};

    cspl::Trajectory<3> traj;
    traj.add_point(x0, v0);
    traj.add_point(x1, v1, 1.);
    traj.add_point(x2, v2, 0.5);

    double total_duration = traj.total_duration();

    std::cout << "Positions:" << std::endl;
    for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
        std::cout << t << ": " << traj.position(t).transpose() << std::endl;
    }

    std::cout << "Velocities:" << std::endl;
    for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
        std::cout << t << ": " << traj.velocity(t).transpose() << std::endl;
    }

    std::cout << "Accelerations:" << std::endl;
    for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
        std::cout << t << ": " << traj.acceleration(t).transpose() << std::endl;
    }

    std::cout << "Compare Accelerations:" << std::endl;
    auto pols = traj.polynomials();
    for (size_t i = 0; i < pols.size() - 1; i++) {
        std::cout << i << ": " << pols[i].polynomial.acceleration(pols[i].duration).transpose() << std::endl;
        std::cout << (i + 1) << ": " << pols[i + 1].polynomial.acceleration(0.).transpose() << std::endl;
    }

    return 0;
}