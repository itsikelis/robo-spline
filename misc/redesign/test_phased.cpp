#include "cubic_hermite_spline.hpp"
#include "phased_trajectory.hpp"

static constexpr size_t NumKnotPoints = 5;
static constexpr size_t Dim = 2;
static constexpr double eps = 1e-6;

static constexpr size_t NumKnotsPerSwing = 2;

int main()
{
    Eigen::VectorXd knots(NumKnotPoints * 2 * Dim);
    Eigen::VectorXd times(NumKnotPoints - 1);

    std::vector<rspl::Phase> phase_seq = {rspl::Phase::Stance, rspl::Phase::Swing, rspl::Phase::Stance, rspl::Phase::Swing};

    rspl::PhasedTrajectory<rspl::CubicHermiteSpline<Dim>, NumKnotPoints> traj(knots, times, phase_seq, NumKnotsPerSwing);

    std::cout << "Positions:" << std::endl;
    for (double t = 0.; t <= (traj.total_duration() + 1e-6); t += 0.1) {
        std::cout << t << ": " << traj.pos(t).transpose() << std::endl;
    }

    // std::cout << "Velocities:" << std::endl;
    // for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
    //     std::cout << t << ": " << traj.velocity(t).transpose() << std::endl;
    // }

    // std::cout << "Accelerations:" << std::endl;
    // for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
    //     std::cout << t << ": " << traj.acceleration(t).transpose() << std::endl;
    // }
    return 0;
}