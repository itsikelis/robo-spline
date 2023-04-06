#include <iostream>

// #include <cspl/cubic_hermite_polynomial.hpp>
#include <cspl/cubic_hermite_polynomial_acc.hpp>
#include <cspl/cubic_hermite_polynomial_reg.hpp>
// #include <cspl/trajectory.hpp>

int main()
{
    // cspl::Trajectory3D traj;
    // traj.add_point({0., 0., 0.}, {1., 0., 1.});
    // traj.add_point({1., 1., 0.}, {1., 1., 1.}, 1.);
    // traj.add_point({2., 2., 0.}, {2., 2., 1.}, 0.5);
    // traj.add_point({0., 0., 0.}, {0., 0., 0.}, 1.);

    // double total_duration = traj.total_duration();

    // std::cout << "Positions:" << std::endl;
    // for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
    //     std::cout << t << ": " << traj.position(t).transpose() << std::endl;
    // }

    // std::cout << "Velocities:" << std::endl;
    // for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
    //     std::cout << t << ": " << traj.velocity(t).transpose() << std::endl;
    // }

    // std::cout << "Accelerations:" << std::endl;
    // for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
    //     std::cout << t << ": " << traj.acceleration(t).transpose() << std::endl;
    // }

    // std::cout << "Compare Accelerations:" << std::endl;
    // auto pols = traj.polynomials();
    // for (size_t i = 0; i < pols.size() - 1; i++) {
    //     std::cout << i << ": " << pols[i].polynomial.acceleration(pols[i].duration).transpose() << std::endl;
    //     std::cout << (i + 1) << ": " << pols[i + 1].polynomial.acceleration(0.).transpose() << std::endl;
    // }

    // std::cout << "Jacobian/Derivative" << std::endl;
    // std::cout << pols[0].polynomial.jac_pos(0.125) << std::endl;

    // return 0;

    cspl::CubicHermitePolynomialReg2D pol1({0., 0.}, {0., 0.}, {1., 1.}, {1., 1.});
    cspl::CubicHermitePolynomialAcc2D pol2({0., 0.}, {0., 0.}, {4., 4.}, {1., 1.});

    // auto acc1 = pol1.acceleration(0.5);

    // std::cout << "Positions:" << std::endl;
    // for (double t = 0.; t <= 1.; t += 0.1) {
    //     std::cout << t << ": " << pol1.position(t).transpose() << std::endl;
    // }

    // std::cout << "Velocities:" << std::endl;
    // for (double t = 0.; t <= 1.; t += 0.1) {
    //     std::cout << t << ": " << pol1.velocity(t).transpose() << std::endl;
    // }

    std::cout << "Accelerations:" << std::endl;
    for (double t = 0.; t <= 1.; t += 0.1) {
        std::cout << t << ": " << pol1.acceleration(t).transpose() << "\t" << pol2.acceleration(t).transpose() << std::endl;
    }
}