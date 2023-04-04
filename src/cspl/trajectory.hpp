#include <Eigen/Dense>

#include <iostream>
#include <vector>

namespace cspl {
    template <unsigned int D>
    class Trajectory {
    public:
        using Vec = typename Eigen::Matrix<double, D, 1>;

        // Struct to store each polynomial in chain and its corresponding duration.
        struct PolynomialTimePair {
            CubicHermitePolynomial<D> polynomial;
            double duration;
        };

        // Default constructor.
        Trajectory() : _total_duration(-1.) {}

        // Initialize a trajectory given the polynomial time durations.
        Trajectory(Eigen::VectorXd durations) : _total_duration(-1.)
        {
            // Fill _polynomials vector with dummy polynomials (all-zeros).
            for (int i = 0; i < durations.size(); i++) {
                _polynomials.push_back({{0, 0, 0, 0}, durations(i)});
            }
        }

        // Get total duration of trajectory.
        double total_duration() const { return _total_duration; }

        // Adds a new target position and velocity and calculates a position-velocity polynomial to get there.
        void add_point(const Vec& next_pos, const Vec& next_vel, double duration = 0.)
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

            _polynomials.push_back({{_last_pos, _last_vel, next_pos, next_vel}, duration});
            _total_duration += duration;

            _last_pos = next_pos;
            _last_vel = next_vel;
        }

        // Adds a new target position and calculate a position-velocity-acceleration polynomial to get there.
        void add_point(const Vec& next_pos, double duration = 0.)
        {
            if (_total_duration < 0.) {
                std::cerr << "You cannot use the acceleration constructor for the first point!" << std::endl;
                return;
            }

            // use the acceleration constructor
            _polynomials.push_back({{_last_pos, _last_vel, _polynomials.back().polynomial.acceleration(1.), next_pos, false}, duration});
            _total_duration += duration;

            _last_pos = next_pos;
            _last_vel = _polynomials.back().polynomial.velocity(1.);
        }

        // Set polynomial parameters
        // full for position-velocity at each point, !full for position-velocity at start-end and positions in-between.
        // regular for position-velocity at start-end, !regular for position-velocity-acceleration at start, position at end.
        void set_params(Eigen::VectorXd params, bool full = true)
        {
            if (full) {
                for (int i = 0; i < _polynomials.size(); i += 2) {
                    double x0 = params.segment(i * 2 * D, D);
                    double v0 = params.segment(i * 2 * D + D, D);
                    double x1 = params.segment((i + 1) * 2 * D, D);
                    double v1 = params.segment((i + 1) * 2 * D + D, D);

                    _polynomials.at(i).polynomial.set_node_node_params({x0, v0, x1, v1});
                }
            }
            else {
                auto x0 = params.segment(0, D);
                auto v0 = params.segment(D, D);
                auto x1 = params.segment(2 * D, D);
                auto acc = _polynomials.at(0).acceleration(1.);
                _polynomials.at(0).polynomial.set_node_node_params({x0, v0, acc, x1}, false);

                // Fill in-between points.
                for (int i = 1; i < _polynomials.size() - 1; i += 2) {
                    auto x0 = params.segment(i * 2 * D, D);
                    auto x1 = params.segment(i * 2 * D + D, D);

                    auto v0 = _polynomials.at(i - 1).polynomial.velocity(1.);
                    auto acc = _polynomials.at(i - 1).polynomial.acceleration(1.);

                    _polynomials.at(i).polynomial.set_node_node_params({x0, v0, acc, x1}, false);
                }

                // Configure last polynomial.

                // second-to-last position begins 3D positions from end of params.
                auto x0 = params.segment(params.size() - 3 * D, D);
                // get second-to-last velocity from second-to-last polynomial.
                auto v0 = _polynomials.at(_polynomials.size() - 2).polynomial.velocity(1.);

                auto x1 = params.segment(params.size() - 2 * D, D);
                auto v1 = params.segment(params.size() - D, D);

                _polynomials.back().polynomial.set_node_node_params({x0, v0, x1, v1});
            }
        }

        // TODO: Get params and Jacobian.

        // Get position at time t.
        Vec position(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _polynomials.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _polynomials[i].duration;
                if (t <= (sum + _epsilon)) {
                    double t_norm = (t - prev_sum) / _polynomials[i].duration;
                    return _polynomials[i].polynomial.position(t_norm);
                }
                prev_sum = sum;
            }

            // If t > trajectory duration, return final trajectory point.
            return _polynomials.back().polynomial.position(1.);
        }

        // Get velocity at time t.
        Vec velocity(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _polynomials.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _polynomials[i].duration;
                if (t <= (sum + _epsilon)) {
                    double t_norm = (t - prev_sum) / _polynomials[i].duration;
                    return _polynomials[i].polynomial.velocity(t_norm);
                }
                prev_sum = sum;
            }

            // If t > trajectory duration, return final trajectory point.
            return _polynomials.back().polynomial.velocity(1.);
        }

        // Get acceleration at time t.
        Vec acceleration(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _polynomials.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _polynomials[i].duration;
                if (t <= (sum + _epsilon)) {
                    double t_norm = (t - prev_sum) / _polynomials[i].duration;
                    return _polynomials[i].polynomial.acceleration(t_norm);
                }
                prev_sum = sum;
            }

            // If t > trajectory duration, return final trajectory point.
            return _polynomials.back().polynomial.acceleration(1.);
        }

        // Get polynomials vector (const ref).
        const std::vector<PolynomialTimePair>& polynomials() const { return _polynomials; }
        // Get polynomials vector (pass-by-reference to modify outside class and avoid copies).
        std::vector<PolynomialTimePair>& polynomials() { return _polynomials; }

    protected:
        static constexpr double _epsilon = 1e-12;
        double _total_duration; // Total duration of trajectory.
        Vec _last_pos, _last_vel; // Target position and velocity of last point entered.
        std::vector<PolynomialTimePair> _polynomials; // Polynomials and their time durations stored in and std::vector.
    };

    // Aliases.
    using Trajectory2D = Trajectory<2>;
    using Trajectory3D = Trajectory<3>;
} // namespace cspl
