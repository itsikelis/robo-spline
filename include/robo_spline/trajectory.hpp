#pragma once

#include <Eigen/Dense>

#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include <robo_spline/cubic_hermite_spline.hpp>
#include <robo_spline/cubic_hermite_spline_acc.hpp>

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
        using PointIndex = typename CubicHermiteSpline<D>::PointIndex;
        using Vec2 = Eigen::Vector2d;
        using Time = double;
        using SplineIndex = unsigned int;

        struct SplineDurationPair {
            std::shared_ptr<CubicHermiteSpline<D>> spline;
            double duration;
        };

    public:
        Trajectory() : _total_duration(-1.) {}

        void clear()
        {
            _spline_duration_pairs.clear();
            _total_duration = -1.;
        }

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
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            return _spline_duration_pairs[idx].spline->position(t_norm);
        }

        /**
         * @brief Get velocity at time t.
         *
         * @param t Time.
         * @return VecD Velocity vector.
         */
        VecD velocity(double t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            return _spline_duration_pairs[idx].spline->velocity(t_norm);
        }

        /**
         * @brief Get velocity at time t.
         *
         * @param t Time.
         * @return VecD Velocity vector.
         */
        VecD acceleration(double t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            return _spline_duration_pairs[idx].spline->acceleration(t_norm);
        }

        /**
         * @brief Get the position derivative of the trajectory at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A pair containing a 4D vector containing the position derivative at the given time amd the current polynomial index.
         */
        std::pair<SplineIndex, Vector> deriv_pos(double t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Vector deriv = _spline_duration_pairs[idx].spline->deriv_pos(t_norm);

            return std::make_pair(idx, deriv);
        }

        std::pair<SplineIndex, double> deriv_pos(double t, PointIndex knot_index) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            double deriv = _spline_duration_pairs[idx].spline->deriv_pos(t_norm, knot_index);

            return std::make_pair(idx, deriv);
        }

        /**
         * @brief Get the velocity derivative of the trajectory at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the velocity derivative at the given time.
         */
        std::pair<SplineIndex, Vector> deriv_vel(double t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Vector deriv = _spline_duration_pairs[idx].spline->deriv_vel(t_norm);

            return std::make_pair(idx, deriv);
        }

        std::pair<SplineIndex, double> deriv_vel(double t, PointIndex knot_index) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            double deriv = _spline_duration_pairs[idx].spline->deriv_vel(t_norm, knot_index);

            return std::make_pair(idx, deriv);
        }

        /**
         * @brief Get the acceleration derivative of the trajectory at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A 4D vector containing the acceleration derivative at the given time.
         */
        std::pair<SplineIndex, Vector> deriv_acc(double t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Vector deriv = _spline_duration_pairs[idx].spline->deriv_acc(t_norm);

            return std::make_pair(idx, deriv);
        }

        std::pair<SplineIndex, double> deriv_acc(double t, PointIndex knot_index) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            double deriv = _spline_duration_pairs[idx].spline->deriv_acc(t_norm, knot_index);

            return std::make_pair(idx, deriv);
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
        std::pair<SplineIndex, Time> normalise_time(double t) const
        {
            if (t < _total_duration) {
                double sum = 0;
                double prev_sum = 0;
                for (int i = 0; i < static_cast<int>(_spline_duration_pairs.size()); i++) {
                    sum += _spline_duration_pairs[i].duration;

                    if (t <= sum - _epsilon) {
                        Time t_norm = (t - prev_sum) / (sum - prev_sum);

                        return std::make_pair(i, t_norm);
                    }

                    prev_sum = sum;
                }
            }
            // Keep this outside an else statement, because for loop may fail to return due to floating point error.
            return std::make_pair(static_cast<int>(_spline_duration_pairs.size() - 1), 1.);
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
