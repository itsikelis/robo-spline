#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "cubic_hermite_spline.hpp"
#include "hermite_spline.hpp"
#include "trajectory.hpp"

Eigen::VectorXd random_uniform_vector(size_t rows, double lower, double upper)
{
    std::uniform_real_distribution<double> unif(lower, upper);
    std::default_random_engine re;
    re.seed(time(0));

    Eigen::VectorXd rand(rows);
    for (size_t i = 0; i < rows; ++i) {
        rand[i] = unif(re);
    }

    return rand;
}

bool test_duration()
{
    static constexpr size_t NumKnotPoints = 5;
    static constexpr size_t Dim = 3;
    std::cout << "Trajectory Duration Test" << std::endl;
    // for (size_t i = 0; i < static_cast<size_t>(times.rows()); ++i) {

    std::uniform_real_distribution<double> unif_knots(-10., 10.);
    std::uniform_real_distribution<double> unif_time(0., 3.);
    std::default_random_engine re;
    re.seed(time(0));

    // const size_t dim = _Traj.dim();

    Eigen::VectorXd knots = random_uniform_vector(NumKnotPoints * 2 * Dim, -10., 10.);
    Eigen::VectorXd times = random_uniform_vector(NumKnotPoints - 1, 0., 3.);

    rspl::Trajectory<rspl::CubicHermiteSpline<Dim>, NumKnotPoints> traj(knots, times);

    std::cout << "Duration 1, Duration 2 " << std::endl;
    std::cout << times.sum() << ", " << traj.total_duration() << " (these should be approx. equal)" << std::endl;

    double eps = 1e-12;
    return (std::abs(times.sum() - traj.total_duration()) < eps);
}

int main()
{
    srand(static_cast<unsigned int>(time(0)));

    rspl::Trajectory<rspl::CubicHermiteSpline<3>, 2> traj1;

    rspl::Trajectory<rspl::CubicHermiteSpline<2>, 2> traj2;

    std::cout << traj1.dim() << ", " << traj2.dim() << std::endl;

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

    // const size_t num_knot_points = 3;
    // const size_t dim = 2;

    // Eigen::Vector2d p0(0., 0.);
    // Eigen::Vector2d v0(0., 0.);

    // Eigen::Vector2d p1(1., 2.);
    // Eigen::Vector2d v1(1., 1.);

    // Eigen::Vector2d p2(2., 3.);
    // Eigen::Vector2d v2(0., 0.);

    // double t0 = 3.;
    // double t1 = 2.;

    // Eigen::VectorXd points(num_knot_points * 2 * dim);
    // Eigen::VectorXd times(num_knot_points - 1);

    // points << p0, v0, p1, v1, p2, v2;
    // times << t0, t1;

    // rspl::Trajectory<rspl::CubicHermiteSpline<dim>, num_knot_points> traj(points, times);

    // // Test derivatives
    // double t = 3.5;
    // size_t idx;
    // rspl::Jacobian jac_pos;
    // std::tie(idx, jac_pos) = traj.jacobian_velocity(t);

    // Eigen::MatrixXd jac_pos_est = Eigen::MatrixXd(dim, points.rows());

    // double eps = 1e-6;
    // for (size_t i = 0.; i < static_cast<size_t>(points.rows()); ++i) {
    //     auto points_p = points;
    //     auto points_m = points;
    //     points_p[i] += eps;
    //     points_m[i] -= eps;

    //     auto traj_p = rspl::Trajectory<rspl::CubicHermiteSpline<dim>, num_knot_points>(points_p, times);
    //     auto traj_m = rspl::Trajectory<rspl::CubicHermiteSpline<dim>, num_knot_points>(points_m, times);

    //     auto pos_p = traj_p.velocity(t);
    //     auto pos_m = traj_m.velocity(t);

    //     jac_pos_est.col(i) = (pos_p - pos_m) / (2 * eps);
    // }

    // std::cout << "Jacobian Block:" << std::endl;
    // std::cout << Eigen::MatrixXd(jac_pos) << std::endl;

    // std::cout << "Jacobian Block Approximation:" << std::endl;
    // std::cout << jac_pos_est << std::endl;

    // std::cout << "Absolute norm of difference: (eps: " << eps << ")" << std::endl;
    // std::cout << abs((jac_pos - jac_pos_est).norm()) << std::endl;
    test_duration();

    return 0;
}