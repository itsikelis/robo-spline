#include "cubic_hermite_spline.hpp"
#include "hermite_spline.hpp"
#include "trajectory.hpp"
#include <iostream>

int main()
{
    // Eigen::VectorXd vec(8);
    // vec << 0., 0., 0., 0., 1., 1., 0., 0.;
    // rspl::CubicHermiteSpline<2> spline(vec, 1.);

    // const size_t num_knot_points = 3;
    // const size_t dim = 2;
    // Eigen::VectorXd points(num_knot_points * 2 * dim);
    // Eigen::VectorXd times(num_knot_points - 1);

    // points << 0., 0., 0., 0., 1., 1., 1., 1., 1., 2., 0., 0.;
    // times << 1., 0.5;

    // rspl::Trajectory<rspl::CubicHermiteSpline<dim>, num_knot_points> traj(points, times);

    // double total_duration
    //     = traj.total_duration();

    // std::cout << "Positions:" << std::endl;
    // for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
    //     std::cout << t << ": " << traj.position(t).transpose() << std::endl;
    // }3

    // std::cout << "Velocities:" << std::endl;
    // for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
    //     std::cout << t << ": " << traj.velocity(t).transpose() << std::endl;
    // }

    // std::cout << "Accelerations:" << std::endl;
    // for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
    //     std::cout << t << ": " << traj.acceleration(t).transpose() << std::endl;
    // }

    // auto splines = traj.splines();

    // // These should be equal.
    // std::cout << "Compare Positions:" << std::endl;
    // for (size_t i = 0; i < splines.size() - 1; i++) {
    //     std::cout << i << ": " << splines[i]->positions(1.).transpose() << std::endl;
    //     std::cout << (i + 1) << ": " << splines[i + 1]->positions(0.).transpose() << std::endl;
    // }

    // // These should be equal.
    // std::cout << "Compare Velocities:" << std::endl;
    // for (size_t i = 0; i < splines.size() - 1; i++) {
    //     std::cout << i << ": " << splines[i]->velocity(1.).transpose() << std::endl;
    //     std::cout << (i + 1) << ": " << splines[i + 1]->velocity(0.).transpose() << std::endl;
    // }

    // // These should be inequal.
    // std::cout << "Compare Accelerations:" << std::endl;
    // for (size_t i = 0; i < splines.size() - 1; i++) {
    //     std::cout << i << ": " << splines[i]->acceleration(1.).transpose() << std::endl;
    //     std::cout << (i + 1) << ": " << splines[i + 1]->acceleration(0.).transpose() << std::endl;
    // }

    // std::cout << "Derivatives:" << std::endl;
    // size_t idx = 0;
    // std::cout << pols[idx].spline->deriv_pos(0.125).transpose() << std::endl
    //           << std::endl;
    // std::cout << pols[idx].spline->deriv_vel(0.125).transpose() << std::endl
    //           << std::endl;
    // std::cout << pols[idx].spline->deriv_acc(0.125).transpose() << std::endl
    //           << std::endl;

    // // Test derivatives
    // double t = 1.25;
    // size_t idx;
    // rspl::Jacobian jac_pos;
    // std::tie(idx, jac_pos) = traj.jacobian_position(t);

    // std::cout << "Spline Index:" << idx << std::endl;
    // std::cout << Eigen::MatrixXd(jac_pos) << std::endl;

    // std::cout << pols[idx].spline->deriv_pos(0.125).transpose();

    const size_t num_knot_points = 2;
    const size_t dim = 2;

    Eigen::Vector2d p0(0., 0.);
    Eigen::Vector2d p1(1., 2.);

    Eigen::Vector2d v0(0., 0.);
    Eigen::Vector2d v1(0., 0.);

    double t0 = 3.;

    Eigen::VectorXd points(num_knot_points * 2 * dim);
    Eigen::VectorXd times(num_knot_points - 1);

    points << p0, v0, p1, v1;
    times << t0;

    rspl::Trajectory<rspl::CubicHermiteSpline<dim>, num_knot_points> traj(points, times);

    // Test derivatives
    double t = 1.2;
    size_t idx;
    rspl::Jacobian jac_pos;
    std::tie(idx, jac_pos) = traj.jacobian_position(t);

    Eigen::MatrixXd jac_pos_est = Eigen::MatrixXd(dim, points.rows());

    double eps = 1e-6;
    for (size_t i = 0.; i < static_cast<size_t>(points.rows()); ++i) {
        auto points_p = points;
        auto points_m = points;
        points_p[i] += eps;
        points_m[i] -= eps;

        auto traj_p = rspl::Trajectory<rspl::CubicHermiteSpline<dim>, num_knot_points>(points_p, times);
        auto traj_m = rspl::Trajectory<rspl::CubicHermiteSpline<dim>, num_knot_points>(points_m, times);

        auto pos_p = traj_p.position(t);
        auto pos_m = traj_m.position(t);

        jac_pos_est.col(i) = (pos_p - pos_m) / (2 * eps);
    }

    std::cout << "Absolute norm of difference: (eps: " << eps << ")" << std::endl;
    std::cout << abs((jac_pos - jac_pos_est).norm()) << std::endl;

    return 0;
}