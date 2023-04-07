#include <Eigen/Dense>

#include <iostream>
#include <vector>

#include <cspl/cubic_hermite_polynomial.hpp>
#include <cspl/cubic_hermite_polynomial_acc.hpp>

namespace cspl {
    template <unsigned int D>
    class Trajectory {
    public:
        using VectorD = typename CubicHermitePolynomial<D>::VectorD;
        using VectorX = typename CubicHermitePolynomial<D>::VectorX;
        using MatrixX = typename CubicHermitePolynomial<D>::MatrixX;
        using VectorTimePair = Eigen::VectorXd;

        // Struct to store each polynomial trajectory and its corresponding duration.
        struct PolynomialTimePair {
            CubicHermitePolynomial<D> polynomial;
            double duration;
        };

        // Default constructor.
        Trajectory() : _total_duration(-1.) {}

        // Initialize a trajectory given the polynomial time durations.
        Trajectory(VectorX durations, bool all_regular = false) : _total_duration(0.)
        {
            // Fill _polynomial_pairs vector with dummy polynomials (all-zeros).
            for (int i = 0; i < durations.size(); i++) {
                if (all_regular || i == 0 || i == (durations.size() - 1))
                    _polynomial_pairs.push_back({CubicHermitePolynomial<D>{}, durations[i]});
                else
                    _polynomial_pairs.push_back({CubicHermitePolynomialAcc<D>{}, durations[i]});
                _total_duration += durations[i];
            }
        }

        // Get total duration of trajectory.
        double total_duration() const { return _total_duration; }

        // Adds a new target position and velocity and calculates a position-velocity polynomial to get there.
        void add_point(const VectorD& next_pos, const VectorD& next_vel, double duration = 0.)
        {
            // Check if it's the first point added to trajectory.
            if (_total_duration < 0.) {
                _total_duration = 0.;
                _last_pos = next_pos;
                _last_vel = next_vel;

                // Initial point begins at time 0.
                if (duration > 0.) {
                    std::cerr << "You cannot have a duration > 0. for the initial point!" << std::endl;
                }
                return;
            }

            _polynomial_pairs.push_back({CubicHermitePolynomial<D>{_last_pos, _last_vel, next_pos, next_vel}, duration});
            _total_duration += duration;

            _last_pos = next_pos;
            _last_vel = next_vel;
        }

        // Adds a new target position and velocity and calculates a position-velocity polynomial to get there.
        void add_point(const VectorD& next_pos, double duration = 1.)
        {
            // Check if it's the first point added to trajectory.
            if (_total_duration < 0.) {
                std::cerr << "You cannot use acceleration constructor for first point!" << std::endl;
                return;
            }

            VectorD a0 = _polynomial_pairs.back().polynomial.acceleration(1.);

            _polynomial_pairs.push_back({CubicHermitePolynomialAcc<D>{_last_pos, _last_vel, a0, next_pos}, duration});

            _total_duration += duration;

            _last_pos = next_pos;
            _last_vel = _polynomial_pairs.back().polynomial.velocity(1.);
        }

        // Get position at time t.
        VectorD position(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _polynomial_pairs.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _polynomial_pairs[i].duration;
                if (t <= (sum + _epsilon)) {
                    double t_norm = (t - prev_sum) / _polynomial_pairs[i].duration;
                    return _polynomial_pairs[i].polynomial.position(t_norm);
                }
                prev_sum = sum;
            }

            // If t > trajectory duration, return final trajectory point.
            return _polynomial_pairs.back().polynomial.position(1.);
        }

        // Get velocity at time t.
        VectorD velocity(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _polynomial_pairs.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _polynomial_pairs[i].duration;
                if (t <= (sum + _epsilon)) {
                    double t_norm = (t - prev_sum) / _polynomial_pairs[i].duration;
                    return _polynomial_pairs[i].polynomial.velocity(t_norm);
                }
                prev_sum = sum;
            }

            // If t > trajectory duration, return final trajectory point.
            return _polynomial_pairs.back().polynomial.velocity(1.);
        }

        // Get acceleration at time t.
        VectorD acceleration(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _polynomial_pairs.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _polynomial_pairs[i].duration;
                if (t <= (sum + _epsilon)) {
                    double t_norm = (t - prev_sum) / _polynomial_pairs[i].duration;
                    return _polynomial_pairs[i].polynomial.acceleration(t_norm);
                }
                prev_sum = sum;
            }

            // If t > trajectory duration, return final trajectory point.
            return _polynomial_pairs.back().polynomial.acceleration(1.);
        }

        // Get polynomials vector (const ref).
        const std::vector<PolynomialTimePair>& polynomials() const { return _polynomial_pairs; }
        // Get polynomials vector (pass-by-reference to modify outside class and avoid copies).
        std::vector<PolynomialTimePair>& polynomials() { return _polynomial_pairs; }

    protected:
        static constexpr double _epsilon = 1e-12;
        double _total_duration; // Total duration of trajectory.
        VectorD _last_pos, _last_vel; // Target position and velocity of last point entered.
        std::vector<PolynomialTimePair> _polynomial_pairs; // Polynomials and their time durations stored in and std::vector.
    };

    // Aliases.
    using Trajectory2D = Trajectory<2>;
    using Trajectory3D = Trajectory<3>;
} // namespace cspl
