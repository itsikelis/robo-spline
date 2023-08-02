#include <iostream>

#include <robo_spline/trajectory.hpp>

int main()
{
    rspl::Trajectory2D traj;
    traj.add_point({0., 0.});
    traj.add_point({1., 1.}, 1.);
    traj.add_point({2., 2.}, 0.5);
    traj.add_point({0., 0.}, 1.);

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

    auto pols = traj.splines();

    std::cout << "Compare Velocities:" << std::endl;
    for (size_t i = 0; i < pols.size() - 1; i++) {
        std::cout << i << ": " << pols[i]->velocity(pols[i]->duration()).transpose() << std::endl;
        std::cout << (i + 1) << ": " << pols[i + 1]->velocity(0.).transpose() << std::endl;
    }

    std::cout << "Compare Accelerations:" << std::endl;
    for (size_t i = 0; i < pols.size() - 1; i++) {
        std::cout << i << ": " << pols[i]->acceleration(pols[i]->duration()).transpose() << std::endl;
        std::cout << (i + 1) << ": " << pols[i + 1]->acceleration(0.).transpose() << std::endl;
    }

    std::cout << "Derivatives:" << std::endl;
    size_t idx = 0;
    std::cout << pols[idx]->deriv_pos(0.125).transpose() << std::endl
              << std::endl;
    std::cout << pols[idx]->deriv_vel(0.125).transpose() << std::endl
              << std::endl;
    std::cout << pols[idx]->deriv_acc(0.125).transpose() << std::endl
              << std::endl;

    return 0;
}
