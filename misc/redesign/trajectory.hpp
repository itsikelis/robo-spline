#pragma once

#include <cstdint>
#include <iostream>
#include <memory>

#include "hermite_spline.hpp"
#include "types.hpp"

namespace rspl {

    template <size_t _Dim>
    class Trajectory {
    public:
        using VecD = typename CubicHermiteSpline<_Dim>::VecD;
        using Spline = CubicHermiteSpline<_Dim>;

    public:
        Trajectory(const Vector& knot_points, const Vector& times) : _total_duration(times.sum()), _num_vars_total(knot_points.rows()), _num_knot_points(_num_vars_total / (2 * _Dim))
        {
            const size_t num_splines = _num_knot_points - 1;

            for (size_t i = 0; i < num_splines; ++i) {
                Eigen::VectorXd spline_knots = knot_points.segment(i * 2 * _Dim, 4 * _Dim);
                Time dt = times[i];
                _splines.push_back(std::make_shared<Spline>(spline_knots, dt));
            }
        }

        ~Trajectory() { _splines.clear(); };

        void clear()
        {
            _splines.clear();
            _total_duration = 0.;
        }

        inline size_t dim() const { return _Dim; }
        inline std::shared_ptr<Spline> spline(size_t idx) const { return _splines[idx]; }
        inline size_t num_knot_points() const { return _num_knot_points; }
        inline double total_duration() const { return _total_duration; }

        inline VecD evaluate(double t, size_t order) const
        {
            switch (order) {
            case 0:
                return pos(t);
            case 1:
                return vel(t);
            case 2:
                return acc(t);
            default:
                std::cerr << "Invalid derivative order!" << std::endl;
                return VecD::Zero();
            }
        }

        // @brief Get position at time t.
        inline VecD pos(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            return _splines[idx]->pos(t_norm);
        }

        // @brief Get velocity at time t.
        inline VecD vel(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            return _splines[idx]->vel(t_norm);
        }

        // @brief Get acceleration at time t.
        inline VecD acc(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            return _splines[idx]->acc(t_norm);
        }

        std::pair<SplineIndex, Jacobian> jac_block(Time t, size_t order) const
        {
            switch (order) {
            case 0:
                return jac_block_pos(t);
            case 1:
                return jac_block_vel(t);
            case 2:
                return jac_block_acc(t);
            default:
                std::cerr << "Invalid derivative order!" << std::endl;
                return std::make_pair(0, Jacobian());
            }
        }

        std::pair<SplineIndex, Jacobian> jac_block_pos(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Jacobian spline_jac = _splines[idx]->jac_block_pos(t_norm);

            Jacobian jac(_num_vars_total, _Dim);
            jac.middleRows(idx * 2 * _Dim, spline_jac.cols()) = spline_jac.transpose();

            return std::make_pair(idx, jac.transpose());
        }

        std::pair<SplineIndex, Jacobian> jac_block_vel(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Jacobian spline_jac = _splines[idx]->jac_block_vel(t_norm);

            Jacobian jac(_num_vars_total, _Dim);
            jac.middleRows(idx * 2 * _Dim, spline_jac.cols()) = spline_jac.transpose();

            return std::make_pair(idx, jac.transpose());
        }

        std::pair<SplineIndex, Jacobian> jac_block_acc(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Jacobian spline_jac = _splines[idx]->jac_block_acc(t_norm);

            Jacobian jac(_num_vars_total, _Dim);
            jac.middleRows(idx * 2 * _Dim, spline_jac.cols()) = spline_jac.transpose();

            return std::make_pair(idx, jac.transpose());
        }

    protected:
        std::pair<SplineIndex, Time> normalise_time(Time t) const
        {
            double sum = 0.;
            double prev_sum = 0.;
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
            return std::make_pair(static_cast<size_t>(_splines.size() - 1), _splines.back()->duration());
        }

    protected:
        static constexpr double _epsilon = 1e-12;

        const double _total_duration;
        const size_t _num_vars_total;
        const size_t _num_knot_points;
        std::vector<std::shared_ptr<Spline>> _splines;
    };
} // namespace rspl