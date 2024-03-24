#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "cubic_hermite_spline.hpp"
#include "phased_trajectory.hpp"

static constexpr size_t Dim = 2;
static constexpr double eps = 1e-6;
static constexpr size_t NumStancePhases = 2;
static constexpr size_t NumSwingPhases = 1;

static constexpr size_t NumKnotsPerSwing = 3;

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

// Test calculated.
template <typename _Traj>
bool test_duration(const _Traj& traj, const Eigen::VectorXd& times)
{
    // std::cout << "Trajectory Duration Test" << std::endl;

    // std::cout << "Duration 1, Duration 2 " << std::endl;
    // std::cout << times.sum() << ", " << traj.total_duration() << " (these should be approx. equal)" << std::endl;
    return (std::abs(times.sum() - traj.total_duration()) < eps);
}

// Test pos and vel at knot points are approximately equal.
template <typename _Traj>
bool test_splines(const _Traj& traj)
{
    bool flag = true;

    for (size_t deriv_order = 0; deriv_order < 2; ++deriv_order) {
        for (size_t i = 0; i < traj.num_splines() - 1; i++) {
            double T = traj.spline(i)->duration();
            auto pos1 = traj.spline(i)->evaluate(T, deriv_order);
            auto pos2 = traj.spline(i + 1)->evaluate(0., deriv_order);
            auto res = std::abs((pos1 - pos2).norm());
            if (res > eps)
                flag = false;
        }
    }
    return flag;
}

// Compare calculated with estimated Jacobians
template <typename _Traj>
bool test_jacobians(const _Traj& traj, const Eigen::VectorXd& knots, const Eigen::VectorXd& times, const std::vector<rspl::Phase>& phase_seq)
{
    bool flag = true;
    Eigen::VectorXd test_times(times.rows() - 1);

    // for (size_t i = 0; i < static_cast<size_t>(times.rows()) - 1; ++i) {
    //     if (i == 0) {
    //         test_times[i] = (times[i + 1] - times[i]) / 2.;
    //         continue;
    //     }
    //     test_times[i] = times[i] + ((times[i + 1] - times[i]) / 2.);
    // }
    // Evaluate all Jacobian orders.
    for (size_t order = 0; order < 3; ++order) {
        for (size_t i = 0; i < static_cast<size_t>(test_times.rows()); ++i) {
            // double t = test_times[i];

            double t = 0.2;

            size_t idx;
            rspl::Jacobian jac;
            std::tie(idx, jac) = traj.jac_block(t, order);

            Eigen::MatrixXd jac_est = Eigen::MatrixXd(Dim, knots.rows());

            for (size_t j = 0.; j < static_cast<size_t>(knots.rows()); ++j) {
                auto knots_p = knots;
                auto knots_m = knots;
                knots_p[j] += eps;
                knots_m[j] -= eps;

                auto traj_p = rspl::PhasedTrajectory<Dim>(knots_p, times, phase_seq, NumKnotsPerSwing);
                auto traj_m = rspl::PhasedTrajectory<Dim>(knots_m, times, phase_seq, NumKnotsPerSwing);

                auto pos_p = traj_p.evaluate(t, order);
                auto pos_m = traj_m.evaluate(t, order);

                jac_est.col(j) = (pos_p - pos_m) / (2 * eps);
            }

            std::cout << "Jacobian Block:" << std::endl;
            std::cout << Eigen::MatrixXd(jac) << std::endl;

            std::cout << "Jacobian Block Approximation:" << std::endl;
            std::cout << jac_est << std::endl;

            std::cout << "Absolute norm of difference: (eps: " << eps << ")" << std::endl;
            std::cout << abs((jac - jac_est).norm()) << std::endl;
            double err = abs((jac - jac_est).norm());

            if (err > 2 * eps) {
                // std::cerr << "Absolute norm of difference: (eps: " << eps << ")" << std::endl;
                // std::cerr << abs((jac - jac_est).norm()) << std::endl;
                flag = false;
                break;
            }
        }
    }

    return flag;
}

int main()
{
    Eigen::VectorXd knots((NumStancePhases * Dim) + (NumSwingPhases * NumKnotsPerSwing * 2 * Dim));
    Eigen::VectorXd phase_times(NumStancePhases + NumSwingPhases);

    knots << 0., 0., 0.25, 0.25, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.75, 0.25, 0.0, 0.0, 1.0, 0.0;

    // knots = random_uniform_vector(knots.rows(), -10., 10.);
    // std::cout << phase_times.rows() << std::endl;
    phase_times << 0.5, 1., 0.5;

    std::vector<rspl::Phase> phase_seq = {rspl::Phase::Stance, rspl::Phase::Swing, rspl::Phase::Stance};

    rspl::PhasedTrajectory<Dim> traj(knots, phase_times, phase_seq, NumKnotsPerSwing);

    // std::cout << "Phase Times: " << phase_times.transpose() << std::endl;
    // std::cout << "Duration: " << traj.total_duration() << std::endl;

    std::cout << "Positions:" << std::endl;
    for (double t = 0.; t <= (traj.total_duration() + 1e-6); t += 0.1) {
        std::cout << traj.pos(t)[1] << ", ";
    }
    std::cout << std::endl;

    // std::cout << "Velocities:" << std::endl;
    // for (double t = 0.; t <= (traj.total_duration() + 1e-6); t += 0.1) {
    //     std::cout << t << ": " << traj.vel(t).transpose() << std::endl;
    // }

    // std::cout << "Accelerations:" << std::endl;
    // for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
    //     std::cout << t << ": " << traj.acceleration(t).transpose() << std::endl;
    // }

    // if (!test_duration(traj, phase_times)) {
    //     std::cerr << "Test Failed: test_duration" << std::endl;
    //     return -1;
    // }

    // if (!test_splines(traj)) {
    //     std::cerr << "Test Failed: test_splines" << std::endl;
    //     return -1;
    // }

    // if (!test_jacobians(traj, knots, phase_times, phase_seq)) {
    //     std::cerr << "Test Failed: test_jacobians" << std::endl;
    //     return -1;
    // }

    return 0;
}