#include <Eigen/Dense>
#include <iostream>

#include "i_cubic_hermite_polynomial.hpp"

namespace cspl {
    template <unsigned int D>
    class Trajectory {
    public:
        using VectorD = typename ICubicHermitePolynomial<D, CubicHermitePolynomialReg<D>>::VectorD;
        using VectorX = typename ICubicHermitePolynomial<D, CubicHermitePolynomialReg<D>>::VectorX;
        using MatrixX = typename ICubicHermitePolynomial<D, CubicHermitePolynomialReg<D>>::MatrixX;
        using VectorTimePair = Eigen::VectorXd;

        // Struct to store each polynomial trajectory and its corresponding duration.
        struct PolynomialTimePair {
            ICubicHermitePolynomial<D, CubicHermitePolynomialReg<D>> polynomial;
            double duration;
        };

        // Default constructor.
        Trajectory() : _total_duration(-1.) {}

        // Initialize a trajectory given the polynomial time durations.
        Trajectory(VectorX durations) : _total_duration(-1.)
        {
            // Fill _polynomial_pairs vector with dummy polynomials (all-zeros).
            for (int i = 0; i < durations.size(); i++) {
                _polynomial_pairs.push_back({{0, 0, 0, 0}, durations(i)});
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

            _polynomial_pairs.push_back({{_last_pos, _last_vel, next_pos, next_vel}, duration});
            _total_duration += duration;

            _last_pos = next_pos;
            _last_vel = next_vel;
        }

        // Adds a new target position and velocity and calculates a position-velocity polynomial to get there.
        void add_point(const VectorD& next_position, double duration = 0.)
        {
            // Check if it's the first point added to trajectory.
            if (_total_duration < 0.) {
                std::cerr << "You cannot use acceleration constructor for first point!" << std::endl;
            }
            return;

            VectorD a0 = _polynomial_pairs.back().acceleration(1.);

            CubicHermitePolynomialAcc<D> chpa(_last_pos, _last_vel, a0, next_position);

            VectorD next_velocity = chpa.velocity(1.);

            _polynomial_pairs.push_back({{_last_pos, _last_vel, next_position, next_velocity}, duration});

            _total_duration += duration;

            _last_pos = next_position;
            _last_vel = next_velocity;
        }

        // Get position at time t.
        Vec position(double t) const
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

    protected:
        static constexpr double _epsilon
            = 1e-12;
        double _total_duration; // Total duration of trajectory.
        VectorD _last_pos, _last_vel; // Target position and velocity of last point entered.
        std::vector<PolynomialTimePair> _polynomial_pairs; // Polynomials and their time durations stored in and std::vector.
    };

    // Aliases.
    using Trajectory2D = Trajectory<2>;
    using Trajectory3D = Trajectory<3>;
} // namespace cspl

////////

//         // Set polynomial parameters
//         // full for position-velocity at each point, !full for position-velocity at start-end and positions in-between.
//         // regular for position-velocity at start-end, !regular for position-velocity-acceleration at start, position at end.
//         void set_params(Eigen::VectorXd params, bool full = true)
//         {
//             if (full) {
//                 for (int i = 0; i < _polynomial_pairs.size(); i += 2) {
//                     Vec x0 = params.segment(i * 2 * D, D);
//                     Vec v0 = params.segment(i * 2 * D + D, D);
//                     Vec x1 = params.segment((i + 1) * 2 * D, D);
//                     Vec v1 = params.segment((i + 1) * 2 * D + D, D);

//                     _polynomial_pairs.at(i).polynomial.set_node_node_params({x0, v0, x1, v1});
//                 }
//             }
//             else {
//                 Vec x0 = params.segment(0, D);
//                 Vec v0 = params.segment(D, D);
//                 Vec x1 = params.segment(2 * D, D);
//                 Vec acc = _polynomial_pairs.at(0).acceleration(1.);
//                 _polynomial_pairs.at(0).polynomial.set_node_node_params({x0, v0, acc, x1}, false);

//                 // Fill in-between points.
//                 for (int i = 1; i < _polynomial_pairs.size() - 1; i += 2) {
//                     Vec x0 = params.segment(i * 2 * D, D);
//                     Vec x1 = params.segment(i * 2 * D + D, D);
//                     Vec v0 = _polynomial_pairs.at(i - 1).polynomial.velocity(1.);
//                     Vec acc = _polynomial_pairs.at(i - 1).polynomial.acceleration(1.);

//                     _polynomial_pairs.at(i).polynomial.set_node_node_params({x0, v0, acc, x1}, false);
//                 }

//                 // Configure last polynomial.

//                 // second-to-last position begins 3D positions from end of params.
//                 auto x0 = params.segment(params.size() - 3 * D, D);
//                 // get second-to-last velocity from second-to-last polynomial.
//                 auto v0 = _polynomial_pairs.at(_polynomial_pairs.size() - 2).polynomial.velocity(1.);

//                 auto x1 = params.segment(params.size() - 2 * D, D);
//                 auto v1 = params.segment(params.size() - D, D);

//                 _polynomial_pairs.back().polynomial.set_node_node_params({x0, v0, x1, v1});
//             }
//         }

//         // Eigen::VectorXd params(bool full = 1) const
//         // {
//         //     Eigen::VectorXd vec;
//         //     if (full) {
//         //         for (int i = 0; i < _polynomial_pairs.size(); i++) {
//         //             vec << _polynomial_pairs.at(i).polynomial.node_params("initial");
//         //         }
//         //         vec << _polynomial_pairs.back().polynomial.node_params("target");
//         //     }
//         //     else {
//         //         vec << _polynomial_pairs.front().polynomial.node_params("initial");
//         //         for (int i = 1; i < _polynomial_pairs.size(); i++) {
//         //             vec << _polynomial_pairs.at(i).polynomial.node_params("p0");
//         //         }
//         //         vec << _polynomial_pairs.back().polynomial.node_params("target");
//         //     }
//         //     return vec;
//         // }

//         // TODO: Get Jacobian.

//             // If t > trajectory duration, return final trajectory point.
//             return _polynomial_pairs.back().polynomial.position(1.);
//         }

//         // Get velocity at time t.
//         Vec velocity(double t) const
//         {
//             double sum = 0;
//             double prev_sum = 0;
//             for (unsigned int i = 0; i < _polynomial_pairs.size(); i++) {
//                 // Add current polynomial's duration to sum and check if t belongs to this interval.
//                 sum += _polynomial_pairs[i].duration;
//                 if (t <= (sum + _epsilon)) {
//                     double t_norm = (t - prev_sum) / _polynomial_pairs[i].duration;
//                     return _polynomial_pairs[i].polynomial.velocity(t_norm);
//                 }
//                 prev_sum = sum;
//             }

//             // If t > trajectory duration, return final trajectory point.
//             return _polynomial_pairs.back().polynomial.velocity(1.);
//         }

//         // Get acceleration at time t.
//         Vec acceleration(double t) const
//         {
//             double sum = 0;
//             double prev_sum = 0;
//             for (unsigned int i = 0; i < _polynomial_pairs.size(); i++) {
//                 // Add current polynomial's duration to sum and check if t belongs to this interval.
//                 sum += _polynomial_pairs[i].duration;
//                 if (t <= (sum + _epsilon)) {
//                     double t_norm = (t - prev_sum) / _polynomial_pairs[i].duration;
//                     return _polynomial_pairs[i].polynomial.acceleration(t_norm);
//                 }
//                 prev_sum = sum;
//             }

//             // If t > trajectory duration, return final trajectory point.
//             return _polynomial_pairs.back().polynomial.acceleration(1.);
//         }

//         // Get polynomials vector (const ref).
//         const std::vector<PolynomialTimePair>& polynomials() const { return _polynomial_pairs; }
//         // Get polynomials vector (pass-by-reference to modify outside class and avoid copies).
//         std::vector<PolynomialTimePair>& polynomials() { return _polynomial_pairs; }
