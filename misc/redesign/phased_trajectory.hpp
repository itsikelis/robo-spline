#pragma once

#include <cstdint>
#include <iostream>
#include <memory>

#include "hermite_spline.hpp"
#include "types.hpp"

namespace rspl {

    enum class Phase {
        Stance,
        Swing,
    };

    template <size_t _Dim>
    class PhasedTrajectory {
    public:
        using VecD = typename CubicHermiteSpline<_Dim>::VecD;
        using Spline = CubicHermiteSpline<_Dim>;

    public:
        PhasedTrajectory(const Vector& knot_points, const Vector& phase_times, const std::vector<Phase>& phase_sequence, size_t knots_per_swing)
            : _phase_sequence(phase_sequence), _num_knots_per_swing(knots_per_swing), _total_duration(phase_times.sum())
        {
            size_t k_idx = 0; // knot_points index
            size_t t_idx = 0; // times index
            for (auto& phase : phase_sequence) {
                if (phase == Phase::Stance) {
                    VecD pos = knot_points.segment(k_idx, _Dim);
                    VecD vel = VecD::Zero();
                    Vector spline_knots(4 * _Dim);
                    spline_knots << pos, vel, pos, vel;

                    Time t = phase_times[t_idx];

                    _splines.push_back(std::make_shared<Spline>(spline_knots, t));

                    // k_idx += _Dim;
                    t_idx += 1;
                }
                else {
                    // If curr_phase == Phase::Swing
                    Time t = phase_times[t_idx] / static_cast<double>(_num_knots_per_swing + 1);

                    // Add initial point (from stance to swing)
                    VecD pos0 = knot_points.segment(k_idx, _Dim);
                    std::cout << "pos0: " << pos0.transpose() << "#";
                    VecD vel0 = VecD::Zero();
                    VecD pos1 = knot_points.segment(k_idx + _Dim, _Dim);
                    VecD vel1 = knot_points.segment(k_idx + 2 * _Dim, _Dim);
                    Vector spline_knots(4 * _Dim);
                    spline_knots << pos0, vel0, pos1, vel1;

                    _splines.push_back(std::make_shared<Spline>(spline_knots, t));
                    k_idx += _Dim;

                    // Add in-between spline points
                    for (size_t i = 0; i < _num_knots_per_swing - 1; ++i) {

                        pos0 = knot_points.segment(k_idx, _Dim);
                        vel0 = knot_points.segment(k_idx + _Dim, _Dim);
                        pos1 = knot_points.segment(k_idx + 2 * _Dim, _Dim);
                        vel1 = knot_points.segment(k_idx + 3 * _Dim, _Dim);
                        Vector spline_knots(4 * _Dim);
                        spline_knots << pos0, vel0, pos1, vel1;

                        _splines.push_back(std::make_shared<Spline>(spline_knots, t));
                        k_idx += 2 * _Dim;
                    }
                    // Add last swing knot (swing to stance)
                    pos0 = knot_points.segment(k_idx, _Dim);
                    vel0 = knot_points.segment(k_idx + _Dim, _Dim);
                    pos1 = knot_points.segment(k_idx + 2 * _Dim, _Dim);
                    vel1 = VecD::Zero();
                    // Vector spline_knots(4 * _Dim);
                    spline_knots << pos0, vel0, pos1, vel1;

                    k_idx += 2 * _Dim;
                    t_idx++;
                }
            }
        }

        ~PhasedTrajectory() { _splines.clear(); };

        void clear()
        {
            _splines.clear();
            _total_duration = 0.;
        }

        // inline size_t dim() const { return _Dim; }
        // inline std::shared_ptr<Spline> spline(size_t idx) const { return _splines[idx]; }
        // inline size_t num_knot_points() const { return _num_knot_points; }
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

        Jacobian jac_block_pos(Time t) const
        {
            std::pair<SplineIndex, Time> pair = normalise_time(t);
            SplineIndex idx = pair.first;
            Time t_norm = pair.second;

            Jacobian spline_jac = _splines[idx]->jac_block_pos(t_norm);

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
        static constexpr double _epsilon{1e-12};
        const std::vector<Phase> _phase_sequence;
        const size_t _num_knots_per_swing;

        double _total_duration{0.};
        std::vector<std::shared_ptr<Spline>> _splines;
    };
} // namespace rspl