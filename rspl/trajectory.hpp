#pragma once

#include <Eigen/Dense>

#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "cubic_hermite_spline.hpp"
#include "cubic_hermite_spline_acc.hpp"

namespace rspl {
    /**
     * @brief Trajectory class.
     *
     * @tparam D The dimensionality of the trajectory.
     */
    template <unsigned int D>
    class Trajectory {
    public:
        using VecD = typename CubicHermiteSpline<D>::VecD;
        using Vector = typename CubicHermiteSpline<D>::Vector;
        using Vec2 = Eigen::Vector2d;

        struct SplineDurationPair {
            std::shared_ptr<CubicHermiteSpline<D>> spline;
            double duration;
        };

        Trajectory() : _total_duration(-1.) {}

        /**
         * @brief Get the total duration of the trajectory.
         *
         * @return double The total duration of the trajectory.
         */
        double total_duration() const { return _total_duration; }

        /**
         * @brief Adds a new target position and velocity and calculates a position-velocity polynomial to get there.
         *
         * @param next_pos The next target position.
         * @param next_vel The next target velocity.
         * @param duration The duration to reach the next target. Default is 0.
         * @note If it's the first point added to trajectory, the initial point begins at time 0.
         * @warning You cannot have a duration > 0. for the initial point!
         */
        void add_point(const VecD& next_pos, const VecD& next_vel, double duration = 0.)
        {
            // Check if it's the first point added to trajectory.
            if (_total_duration < 0.) {

                if (duration != 0.) {
                    std::cerr << "You cannot have a duration > 0. for the initial point! Defaulting to 0." << std::endl;
                }

                _last_pos = next_pos;
                _last_vel = next_vel;
                _total_duration = 0.;

                return;
            }

            _spline_duration_pairs.push_back({std::make_shared<CubicHermiteSpline<D>>(_last_pos, _last_vel, next_pos, next_vel), duration});
            _total_duration += duration;

            _last_pos = next_pos;
            _last_vel = next_vel;
        }

        /**
         * @brief Adds a new target position and velocity and calculates a position-velocity polynomial to get there.
         *
         * @param next_pos Target position to reach.
         * @param duration Duration of the trajectory to reach the target position.
         * @warning You cannot use acceleration contructor for the initial point!
         */
        void add_point(const VecD& next_pos, double duration = 0.)
        {
            // Check if it's the first point added to trajectory.
            if (_total_duration < 0.) {
                std::cerr << "No initial velocity and acceleration specified for first point in Trajectory! Assuming 0." << std::endl;
                _last_vel = VecD::Zero();

                if (duration != 0.) {
                    std::cerr << "You cannot have a duration > 0. for the initial point! Defaulting to 0." << std::endl;
                }

                _last_pos = next_pos;
                _last_vel = VecD::Zero();
                _total_duration = 0.;

                return;
            }

            VecD a0;
            if (_total_duration == 0) // We assume zero initial acceleration!
                a0 = VecD::Zero();
            else
                a0 = _spline_duration_pairs.back().spline->acceleration(1.);

            _spline_duration_pairs.push_back({std::make_shared<CubicHermiteSplineAcc<D>>(_last_pos, _last_vel, a0, next_pos), duration});

            _total_duration += duration;

            _last_pos = next_pos;
            _last_vel = _spline_duration_pairs.back().spline->velocity(1.);
        }

        /**
         * @brief Get position at time t.
         *
         * @param t Time.
         * @return VecD Position vector.
         */
        VecD position(double t) const
        {
            if (t <= _total_duration) {
                IndexTimePair pair = t_normalised_and_spline_index(t);
                return _spline_duration_pairs[pair.spline_idx].spline->position(pair.t_norm);
            }
            else {
                return _spline_duration_pairs.back().spline->position(1.);
            }
        }

        /**
         * @brief Get velocity at time t.
         *
         * @param t Time.
         * @return VecD Velocity vector.
         */
        VecD velocity(double t) const
        {
            if (t <= _total_duration) {
                IndexTimePair pair = t_normalised_and_spline_index(t);
                return _spline_duration_pairs[pair.spline_idx].spline->velocity(pair.t_norm);
            }
            else {
                return _spline_duration_pairs.back().spline->velocity(1.);
            }
        }

        /**
         * @brief Get velocity at time t.
         *
         * @param t Time.
         * @return VecD Velocity vector.
         */
        VecD acceleration(double t) const
        {
            if (t <= _total_duration) {
                IndexTimePair pair = t_normalised_and_spline_index(t);
                return _spline_duration_pairs[pair.spline_idx].spline->acceleration(pair.t_norm);
            }
            else {
                return _spline_duration_pairs.back().spline->acceleration(1.);
            }
        }

        /**
         * @brief Get the position derivative of the trajectory at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A pair containing a 4D vector containing the position derivative at the given time amd the current polynomial index.
         */
        std::pair<Vector, int> deriv_pos(double t) const
        {
            if (t <= _total_duration) {
                IndexTimePair pair = t_normalised_and_spline_index(t);
                return std::make_pair(_spline_duration_pairs[pair.spline_idx].spline->deriv_pos(pair.t_norm), pair.spline_idx);
            }
            else {
                return std::make_pair(_spline_duration_pairs.back().spline->deriv_pos(1.), _spline_duration_pairs.size() - 1);
            }
        }

        /**
         * @brief Get the velocity derivative of the trajectory at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the velocity derivative at the given time.
         */
        std::pair<Vector, int> deriv_vel(double t) const
        {
            if (t <= _total_duration) {
                IndexTimePair pair = t_normalised_and_spline_index(t);
                return std::make_pair(_spline_duration_pairs[pair.spline_idx].spline->deriv_vel(pair.t_norm), pair.spline_idx);
            }
            else {
                return std::make_pair(_spline_duration_pairs.back().spline->deriv_vel(1.), _spline_duration_pairs.size() - 1);
            }
        }

        /**
         * @brief Get the acceleration derivative of the trajectory at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the acceleration derivative at the given time.
         */
        std::pair<Vector, int> deriv_acc(double t) const
        {
            if (t <= _total_duration) {
                IndexTimePair pair = t_normalised_and_spline_index(t);
                return std::make_pair(_spline_duration_pairs[pair.spline_idx].spline->deriv_acc(pair.t_norm), pair.spline_idx);
            }
            else {
                return std::make_pair(_spline_duration_pairs.back().spline->deriv_acc(1.), _spline_duration_pairs.size() - 1);
            }
        }

        /**
         * @brief Get the polynomials vector of the trajectory (unmodifiable).
         *
         * @return const std::vector<PolynomialTimePair>& Reference to the vector of polynomial-time pairs.
         */
        const std::vector<SplineDurationPair>& splines() const { return _spline_duration_pairs; }

        /**
         * @brief Get the polynomials vector of the trajectory (modifiable).
         *
         * @note pass-by-reference to modify outside class and avoid copies.
         *
         * @return std::vector<PolynomialTimePair>& Reference to the vector of polynomial-time pairs.
         */
        std::vector<SplineDurationPair>& splines() { return _spline_duration_pairs; }

    protected:
        struct IndexTimePair {
            int spline_idx;
            double t_norm;
        };

        IndexTimePair t_normalised_and_spline_index(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (int i = 0; i < static_cast<int>(_spline_duration_pairs.size()); i++) {
                sum += _spline_duration_pairs[i].duration;

                if (t <= sum - _epsilon) {
                    double t_norm = (t - prev_sum) / (sum - prev_sum);

                    return {i, t_norm};
                }

                prev_sum = sum;
            }
            // If time is not found, return last spline index and total time.
            return {static_cast<int>(_spline_duration_pairs.size() - 1), _total_duration};
        }

    protected:
        static constexpr double _epsilon = 1e-12;
        double _total_duration; // Total duration of trajectory.
        VecD _last_pos, _last_vel; // Target position and velocity of last point entered.
        std::vector<SplineDurationPair> _spline_duration_pairs; // Polynomials and their time durations stored in and std::vector.
    };

    // Aliases.
    using Trajectory2D = Trajectory<2>;
    using Trajectory3D = Trajectory<3>;
} // namespace rspl
