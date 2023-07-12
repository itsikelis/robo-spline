#include <Eigen/Dense>

#include <iostream>
#include <memory>
#include <vector>

#include <cspl/cubic_hermite_polynomial.hpp>
#include <cspl/cubic_hermite_polynomial_acc.hpp>

namespace cspl {
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

        /**
         * @brief Struct to store each polynomial trajectory and its corresponding duration.
         */
        struct PolynomialTimePair {
            std::shared_ptr<CubicHermiteSpline<D>> polynomial;
            double duration;
        };

        /**
         * @brief Default constructor for Trajectory class.
         */
        Trajectory() : _total_duration(-1.) {}

        /**
         * @brief Constructor for Trajectory class that initializes the trajectory given polynomial time durations.
         *
         * @param durations A vector containing the duration of each polynomial.
         * @param all_regular A flag indicating if all polynomials should be regular.
         */
        Trajectory(const std::vector<double>& durations, bool all_regular = false) : _total_duration(0.)
        {
            // Fill _polynomial_pairs vector with dummy polynomials (all-zeros).
            for (size_t i = 0; i < durations.size(); i++) {
                // Last polynomial cannot be acceleration only.
                if (all_regular || i == (durations.size() - 1))
                    _polynomial_pairs.push_back({std::make_shared<CubicHermiteSpline<D>>(VecD::Zero(), VecD::Zero(), VecD::Zero(), VecD::Zero()), durations[i]});
                else
                    _polynomial_pairs.push_back({std::make_shared<CubicHermiteSplineAcc<D>>(VecD::Zero(), VecD::Zero(), VecD::Zero(), VecD::Zero()), durations[i]});
                _total_duration += durations[i];
            }
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
                _total_duration = 0.;
                _last_pos = next_pos;
                _last_vel = next_vel;

                // Initial point begins at time 0.
                if (duration > 0.) {
                    std::cerr << "You cannot have a duration > 0. for the initial point!" << std::endl;
                }
                return;
            }

            _polynomial_pairs.push_back({std::make_shared<CubicHermiteSpline<D>>(_last_pos, _last_vel, next_pos, next_vel), duration});
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
        void add_point(const VecD& next_pos, double duration = 1.)
        {
            // Check if it's the first point added to trajectory.
            if (_total_duration < 0.) {
                std::cerr << "You cannot use acceleration constructor for first point!" << std::endl;
                return;
            }

            VecD a0;
            if (_polynomial_pairs.size() == 0) // We assume zero initial acceleration!
                a0 = VecD::Zero();
            else
                a0 = _polynomial_pairs.back().polynomial->acceleration(1.);

            _polynomial_pairs.push_back({std::make_shared<CubicHermiteSplineAcc<D>>(_last_pos, _last_vel, a0, next_pos), duration});

            _total_duration += duration;

            _last_pos = next_pos;
            _last_vel = _polynomial_pairs.back().polynomial->velocity(1.);
        }

        /**
         * @brief Get position at time t.
         *
         * @param t Time.
         * @return VecD Position vector.
         */
        VecD position(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _polynomial_pairs.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _polynomial_pairs[i].duration;
                if (t <= (sum + _epsilon)) {
                    double t_norm = (t - prev_sum) / _polynomial_pairs[i].duration;
                    return _polynomial_pairs[i].polynomial->position(t_norm);
                }
                prev_sum = sum;
            }

            // If t > trajectory duration, return final trajectory point.
            return _polynomial_pairs.back().polynomial->position(1.);
        }

        /**
         * @brief Get velocity at time t.
         *
         * @param t Time at which to get the velocity.
         * @return Velocity vector at time t.
         */
        VecD velocity(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _polynomial_pairs.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _polynomial_pairs[i].duration;
                if (t <= (sum + _epsilon)) {
                    double t_norm = (t - prev_sum) / _polynomial_pairs[i].duration;
                    return _polynomial_pairs[i].polynomial->velocity(t_norm);
                }
                prev_sum = sum;
            }

            // If t > trajectory duration, return final trajectory point.
            return _polynomial_pairs.back().polynomial->velocity(1.);
        }

        /**
         * @brief Get acceleration at time t.
         *
         * @param t Time at which to get the acceleration.
         * @return Acceleration vector at time t.
         */
        VecD acceleration(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _polynomial_pairs.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _polynomial_pairs[i].duration;
                if (t <= (sum + _epsilon)) {
                    double t_norm = (t - prev_sum) / _polynomial_pairs[i].duration;
                    return _polynomial_pairs[i].polynomial->acceleration(t_norm);
                }
                prev_sum = sum;
            }

            // If t > trajectory duration, return final trajectory point.
            return _polynomial_pairs.back().polynomial->acceleration(1.);
        }

        /**
         * @brief Get the polynomials vector of the trajectory.
         *
         * @return const std::vector<PolynomialTimePair>& Reference to the vector of polynomial-time pairs.
         */
        const std::vector<PolynomialTimePair>& polynomials() const { return _polynomial_pairs; }

        /**
         * @brief Get the polynomials vector of the trajectory (modifiable).
         *
         * @note pass-by-reference to modify outside class and avoid copies.
         *
         * @return std::vector<PolynomialTimePair>& Reference to the vector of polynomial-time pairs.
         */
        std::vector<PolynomialTimePair>& polynomials() { return _polynomial_pairs; }

    protected:
        static constexpr double _epsilon = 1e-12;
        double _total_duration; // Total duration of trajectory.
        VecD _last_pos, _last_vel; // Target position and velocity of last point entered.
        std::vector<PolynomialTimePair> _polynomial_pairs; // Polynomials and their time durations stored in and std::vector.
    };

    // Aliases.
    using Trajectory2D = Trajectory<2>;
    using Trajectory3D = Trajectory<3>;
} // namespace cspl
