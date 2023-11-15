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
        using Jacobian = typename CubicHermiteSpline<D>::Jacobian;
        using Time = double;
        using SplineIndex = unsigned int;
        using SplinePtr = std::shared_ptr<CubicHermiteSpline<D>>;

    public:
        Trajectory() : _total_duration(-1.) {}

        void clear()
        {
            _splines.clear();
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

            _splines.push_back(std::make_shared<CubicHermiteSpline<D>>(_last_pos, _last_vel, next_pos, next_vel, duration));
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
                // std::cerr << "No initial velocity and acceleration specified for first point in Trajectory! Assuming 0." << std::endl;
                _last_vel = VecD::Zero();

                if (duration != 0.) {
                    // std::cerr << "You cannot have a duration > 0. for the initial point! Defaulting to 0." << std::endl;
                }

                _last_pos = next_pos;
                _last_vel = VecD::Zero();
                _total_duration = 0.;

                return;
            }

            VecD a0;
            if (_total_duration < 1e-12) // We assume zero initial acceleration!
                a0 = VecD::Zero();
            else
                a0 = _splines.back()->acceleration(_splines.back()->duration());

            _splines.push_back(std::make_shared<CubicHermiteSplineAcc<D>>(_last_pos, _last_vel, a0, next_pos, duration));

            _total_duration += duration;

            _last_pos = next_pos;
            _last_vel = _splines.back()->velocity(_splines.back()->duration());
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

            return _splines[idx]->position(t_norm);
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

            return _splines[idx]->velocity(t_norm);
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

            return _splines[idx]->acceleration(t_norm);
        }

        /**
         * @brief Get the position derivative of the trajectory at time t.
         *
         * @param t The time to evaluate the derivative at.
         * @return A pair containing a 4D vector containing the position derivative at the given time and the current polynomial index.
         */
        std::pair<SplineIndex, Vector> deriv_pos(double t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Vector deriv = _splines[idx]->deriv_pos(t_norm);

            return std::make_pair(idx, deriv);
        }

        /**
         * @brief Get the position jacobian of the trajectory at time t.
         *
         * @param t The time to evaluate the derjacobianivative at.
         * @return A pair containing the position jacobian at the given time and the current polynomial index.
         */
        std::pair<SplineIndex, Jacobian> jacobian_pos(double t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Jacobian deriv = _splines[idx]->jacobian_pos(t_norm);

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

            Vector deriv = _splines[idx]->deriv_vel(t_norm);

            return std::make_pair(idx, deriv);
        }

        /**
         * @brief Get the velocity jacobian of the trajectory at time t.
         *
         * @param t The time to evaluate the derjacobianivative at.
         * @return A pair containing the velocity jacobian at the given time and the current polynomial index.
         */
        std::pair<SplineIndex, Jacobian> jacobian_vel(double t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Jacobian deriv = _splines[idx]->jacobian_vel(t_norm);

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

            Vector deriv = _splines[idx]->deriv_acc(t_norm);

            return std::make_pair(idx, deriv);
        }

        /**
         * @brief Get the acceleration jacobian of the trajectory at time t.
         *
         * @param t The time to evaluate the derjacobianivative at.
         * @return A pair containing the acceleration jacobian at the given time and the current polynomial index.
         */
        std::pair<SplineIndex, Jacobian> jacobian_acc(double t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Jacobian deriv = _splines[idx]->jacobian_acc(t_norm);

            return std::make_pair(idx, deriv);
        }

        /**
         * @brief Get the polynomials vector of the trajectory (unmodifiable).
         *
         * @return const std::vector<PolynomialTimePair>& Reference to the vector of polynomial-time pairs.
         */
        const std::vector<SplinePtr>& splines() const { return _splines; }

        /**
         * @brief Get the polynomials vector of the trajectory (modifiable).
         *
         * @note pass-by-reference to modify outside class and avoid copies.
         *
         * @return std::vector<PolynomialTimePair>& Reference to the vector of polynomial-time pairs.
         */
        std::vector<SplinePtr>& splines() { return _splines; }

    protected:
        std::pair<SplineIndex, Time> normalise_time(double t) const
        {
            if (t < _total_duration) {
                double sum = 0;
                double prev_sum = 0;
                for (int i = 0; i < static_cast<int>(_splines.size()); i++) {
                    sum += _splines[i]->duration();

                    if (t <= sum - _epsilon) {
                        Time t_norm = (t - prev_sum);

                        return std::make_pair(i, t_norm);
                    }

                    prev_sum = sum;
                }
            }
            // Keep this outside an else statement, because for loop may fail to return due to floating point error.
            return std::make_pair(static_cast<int>(_splines.size() - 1), _splines.back()->duration());
        }

    protected:
        static constexpr double _epsilon = 1e-12;
        double _total_duration; // Total duration of trajectory.
        VecD _last_pos, _last_vel; // Target position and velocity of last point entered.
        std::vector<SplinePtr> _splines; // Polynomials stored in and std::vector.
    };

    // Aliases.
    using Trajectory2D = Trajectory<2>;
    using Trajectory3D = Trajectory<3>;
} // namespace rspl
