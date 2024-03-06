#include <cstdint>
#include <iostream>
#include <memory>

#include "hermite_spline.hpp"
#include "types.hpp"

namespace rspl {
    template <typename _Spline, size_t _NumKnots>
    class Trajectory {
    public:
        using VecD = typename _Spline::VecD;
        using SplinePtr = std::shared_ptr<_Spline>;
        // Number of required knot points for spline (4 for cubic, 6 for quintic etc).
        // using NumSplineReqKnots = typename _Spline<_Dim>::NumReqKnots;

    public:
        Trajectory() : _total_duration(-1.) {}

        Trajectory(const Vector& knot_points, const Vector& times) : _num_vars_total(knot_points.rows())
        {
            // TODO:
            // assert knot_points and times size are correct.
            // assert NumKnots - 1 == times.size()

            // Keep num of total vars

            // Create a temp spline object to get info needed.
            _Spline temp;
            _spline_type = temp.type();
            _dim = temp.dimension();

            generate_trajectory_from_points(knot_points, times, _spline_type);
        }

        inline VecD evaluate(Time t, size_t order)
        {
            switch (order) {
            case 0:
                return position(t);
                break;
            case 1:
                return velocity(t);
                break;
            case 2:
                return acceleration(t);
                break;
            default:
                std::cerr << "Invalid derivative order!" << std::endl;
                return VecD::Zero();
                break;
            }
        }

        // @brief Get position at time t.
        VecD position(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            return _splines[idx]->position(t_norm);
        }

        // @brief Get velocity at time t.
        VecD velocity(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            return _splines[idx]->velocity(t_norm);
        }

        // @brief Get acceleration at time t.
        VecD acceleration(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            return _splines[idx]->acceleration(t_norm);
        }

        std::pair<SplineIndex, Jacobian> jacobian_position(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Jacobian spline_jac = _splines[idx]->jacobian_position(t_norm);
            std::cout << "Jac Rows, Cols:" << spline_jac.rows() << ", " << spline_jac.cols() << std::endl;

            size_t col_offset = jac_column_offset(_spline_type);
            Jacobian jac(_num_vars_total, _dim);
            jac.middleRows(idx * _dim * col_offset, spline_jac.cols()) = spline_jac.transpose();

            return std::make_pair(idx, jac.transpose());
        }

        std::pair<SplineIndex, Jacobian> jacobian_velocity(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Jacobian spline_jac = _splines[idx]->jacobian_velocity(t_norm);

            Jacobian jac(_num_vars_total, _dim);
            jac.middleRows(idx * 2 * _dim, spline_jac.cols()) = spline_jac.transpose();

            return std::make_pair(idx, jac.transpose());
        }

        std::pair<SplineIndex, Jacobian> jacobian_acceleration(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Jacobian spline_jac = _splines[idx]->jacobian_acceleration(t_norm);

            Jacobian jac(_num_vars_total, _dim);
            jac.middleRows(idx * 2 * _dim, spline_jac.cols()) = spline_jac.transpose();

            return std::make_pair(idx, jac.transpose());
        }

        void clear()
        {
            _splines.clear();
            _total_duration = -1.;
        }

        std::vector<SplinePtr>& splines() { return _splines; }

        double total_duration() const { return _total_duration; }

        ~Trajectory() { _splines.clear(); };

    protected:
        std::pair<SplineIndex, Time> normalise_time(Time t) const
        {
            double sum = 0;
            double prev_sum = 0;
            size_t iters = _splines.size();
            for (size_t i = 0; i < iters; ++i) {
                sum += _splines[i]->duration();

                if (t <= sum - _epsilon) {
                    Time t_norm = (t - prev_sum);
                    return std::make_pair(i, t_norm);
                }

                prev_sum = sum;
            }

            // If t > _total_duration, return final time of last spline.
            return std::make_pair(static_cast<int>(_splines.size() - 1), _splines.back()->duration());
        }

        inline void generate_trajectory_from_points(const Vector& knot_points, const Vector& times, SplineType type)
        {
            const size_t num_splines = _NumKnots - 1;
            switch (type) {
            case SplineType::CubicHermite:
                // Break up knot_points and .
                for (size_t i = 0; i < num_splines; ++i) {
                    Eigen::VectorXd spline_knots = knot_points.segment(i * 2 * _dim, 4 * _dim);
                    Time dt = times(i);
                    _splines.push_back(std::make_shared<_Spline>(spline_knots, dt));
                    _total_duration += dt;
                }
                break;
            case SplineType::QuinticHermite:
                // Add points in quintic hermite spline.
                break;
            default:
                std::cerr << "Error, unsuported spline type!" << std::endl;
                break;
            }
        }

        inline size_t jac_column_offset(SplineType type) const
        {
            if (type == SplineType::CubicHermite) {
                return 2;
            }
            else if (type == SplineType::QuinticHermite) {
                return 3;
            }
            else {
                std::cerr << "Wrong Spline Type" << std::endl;
                return 2;
            }
        }

    protected:
        static constexpr double _epsilon = 1e-12;
        const size_t _num_vars_total{0};
        SplineType _spline_type;
        size_t _dim;

        Time _total_duration;
        std::vector<SplinePtr> _splines;
    };
} // namespace rspl
