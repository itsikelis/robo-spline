#include "cubic_hermite_spline.hpp"
#include "phased_trajectory.hpp"

static constexpr size_t Dim = 2;
static constexpr double eps = 1e-6;
static constexpr size_t NumStancePhases = 2;
static constexpr size_t NumSwingPhases = 1;

static constexpr size_t NumKnotsPerSwing = 3;

int main()
{
    Eigen::VectorXd knots((NumStancePhases * Dim) + (NumSwingPhases * NumKnotsPerSwing * 2 * Dim));
    Eigen::VectorXd phase_times(NumStancePhases + NumSwingPhases);

    knots << 0., 0., 0.25, 0.25, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.75, 0.75, 0.0, 0.0, 0., 0.;
    std::cout << phase_times.rows() << std::endl;
    phase_times << 0.5, 1., 0.5;

    std::vector<rspl::Phase> phase_seq = {rspl::Phase::Stance, rspl::Phase::Swing, rspl::Phase::Stance};

    rspl::PhasedTrajectory<Dim> traj(knots, phase_times, phase_seq, NumKnotsPerSwing);

    std::cout << "Phase Times: " << phase_times.transpose() << std::endl;
    std::cout << "Duration: " << traj.total_duration() << std::endl;

    // std::cout << "Positions:" << std::endl;
    // for (double t = 0.; t <= (traj.total_duration() + 1e-6); t += 0.1) {
    //     std::cout << traj.pos(t)[0] << ", ";
    // }
    // std::cout << std::endl;

    std::cout << "Velocities:" << std::endl;
    for (double t = 0.; t <= (traj.total_duration() + 1e-6); t += 0.1) {
        std::cout << t << ": " << traj.vel(t).transpose() << std::endl;
    }

    // std::cout << "Accelerations:" << std::endl;
    // for (double t = 0.; t <= (total_duration + 1e-6); t += 0.1) {
    //     std::cout << t << ": " << traj.acceleration(t).transpose() << std::endl;
    // }
    return 0;
}