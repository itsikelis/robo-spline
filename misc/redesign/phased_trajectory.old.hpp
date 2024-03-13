#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

#include "hermite_spline.hpp"
#include "types.hpp"

namespace rspl {

    enum class Phase {
        Stance,
        Swing,
    };

    struct PhasedTuple {
        SplineIndex spline_idx;
        Phase phase;
        Time time_norm;
    };

    template <typename _Spline, size_t _NumKnots>
    class PhasedTrajectory {
    public:
        using VecD = typename _Spline::VecD;
        using SplinePtr = std::shared_ptr<_Spline>;

        PhasedTrajectory(const Vector& knot_points, const Vector& times, const std::vector<Phase>& phase_sequence, size_t knots_per_swing)
            : _spline_type(_Spline::Type),
              _dim(_Spline::Dim),
              _num_vars_total(knot_points.rows()),
              _total_duration(times.sum()),
              _phase_sequence(phase_sequence),
              _num_knots_per_swing(knots_per_swing)
        {
            // TODO:
            // assert knot_points and times size are correct.
            // assert NumKnots - 1 == times.size()

            generate_trajectory_from_points(knot_points, times, phase_sequence);
        }

        void generate_trajectory_from_points(const Vector& knot_points, const Vector& times, const std::vector<Phase>& phase_sequence)
        {
            size_t k_idx = 0; // knot_points index
            size_t t_idx = 0; // times index
            for (size_t phase_idx = 0; phase_idx < _phase_sequence.size(); ++phase_idx) {
                Phase curr_phase = _phase_sequence[phase_idx];
                if (curr_phase == Phase::Stance) {

                    add_stance_knots(knot_points, times, k_idx, t_idx);
                    k_idx += _Spline::Dim;
                    t_idx += 1;
                }
                else {
                    // If curr_phase == Phase::Swing
                    // Add initial point (from stance to swing)
                    add_init_swing_knots(knot_points, times, k_idx, t_idx);
                    k_idx += _Spline::Dim;

                    // Add in-between spline points
                    for (size_t i = 1; i < _num_knots_per_swing - 1; ++i) {

                        add_intermediate_swing_knots(knot_points, times, k_idx, t_idx);
                        k_idx += 2 * _Spline::Dim;
                    }

                    // Add last swing knot with stance.
                    add_final_swing_knots(knot_points, times, k_idx, t_idx);

                    k_idx += 2 * _Spline::Dim;
                    t_idx++;
                }

                curr_phase = phase_sequence[phase_idx];
            }
        }

        void reset(const Vector& knot_points, const Vector& times)
        {
            _splines.clear();
            _total_duration = 0.;
        }

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

        inline Time total_duration() const { return _total_duration; };

        // @brief Get position at time t.
        inline VecD pos(Time t) const
        {
            auto tuple = normalise_time(t);
            SplineIndex idx = tuple.spline_idx;
            Time t_norm = tuple.time_norm;

            return _splines[idx]->pos(t_norm);
        }

        // @brief Get velocity at time t.
        inline VecD vel(Time t) const
        {
            auto tuple = normalise_time(t);
            SplineIndex idx = tuple.spline_idx;
            Time t_norm = tuple.time_norm;

            return _splines[idx]->vel(t_norm);
        }

        // @brief Get acceleration at time t.
        inline VecD acc(Time t) const
        {
            auto tuple = normalise_time(t);
            SplineIndex idx = tuple.spline_idx;
            Time t_norm = tuple.time_norm;

            return _splines[idx]->acc(t_norm);
        }

        // std::pair<SplineIndex, Jacobian> jac_block(Time t, size_t order) const
        // {
        //     switch (order) {
        //     case 0:
        //         return jac_block_pos(t);
        //     case 1:
        //         return jac_block_vel(t);
        //     case 2:
        //         return jac_block_acc(t);
        //     default:
        //         std::cerr << "Invalid derivative order!" << std::endl;
        //         return std::make_pair(0, Jacobian());
        //     }
        // }

        // std::pair<SplineIndex, Jacobian> jac_block_pos(Time t) const
        // {
        //     std::pair<SplineIndex, Time> pair = normalise_time(t);
        //     SplineIndex idx = pair.first;
        //     Time t_norm = pair.second;

        //     Jacobian spline_jac = _splines[idx]->jac_block_pos(t_norm);

        //     size_t col_offset = jac_column_offset(_spline_type);
        //     Jacobian jac(_num_vars_total, _dim);
        //     jac.middleRows(idx * _dim * col_offset, spline_jac.cols()) = spline_jac.transpose();

        //     return std::make_pair(idx, jac.transpose());
        // }

        // Phase phase_at_time(Time t) const
        // {
        //     double sum = 0.;
        //     double prev_sum = 0.;
        //     size_t iters = _phase_sequence.size();
        //     for (size_t i = 0; i < iters; ++i) {
        //         sum += phase_times[i];

        //         if (t <= sum - _epsilon) {
        //             Time t_norm = (t - prev_sum);
        //             return _phase_sequence[i];
        //         }

        //         prev_sum = sum;
        //     }

        //     // If t > _total_duration, return final phase.
        //     return _phase_sequence.back();
        // }

        PhasedTuple normalise_time(Time t) const
        {
            double sum = 0.;
            double prev_sum = 0.;
            size_t iters = _splines.size();
            for (size_t i = 0; i < iters; ++i) {
                sum += _splines[i]->duration();

                if (t <= sum - _epsilon) {
                    Time t_norm = (t - prev_sum);
                    return {i, _phase_sequence[i], t_norm};
                }

                prev_sum = sum;
            }

            // If t > _total_duration, return final time of last spline.
            return {static_cast<size_t>(_splines.size() - 1), _phase_sequence.back(), _splines.back()->duration()};
        }

    protected:
        void add_stance_knots(const Vector& knot_points, const Vector& times, size_t k_idx, size_t t_idx)
        {
            VecD pos = knot_points.segment(k_idx, _dim);
            VecD vel = VecD::Zero();

            VecD spline_knots((_Spline::Order + 1) * _dim);
            spline_knots << pos, vel, pos, vel;

            Time t = times[t_idx];

            _splines.push_back(std::make_shared<_Spline>(spline_knots, t));
        }

        void add_init_swing_knots(const Vector& knot_points, const Vector& times, size_t k_idx, size_t t_idx)
        {
            Time t = times[t_idx] / static_cast<double>(_num_knots_per_swing + 1);

            VecD pos1 = knot_points.segment(k_idx, _Spline::Dim);
            VecD vel1 = VecD::Zero();
            VecD pos2 = knot_points.segment(k_idx + _Spline::Dim, _Spline::Dim);
            VecD vel2 = knot_points.segment(k_idx + 2 * _Spline::Dim, _Spline::Dim);

            Vector spline_knots((_Spline::Order + 1) * _dim);
            spline_knots << pos1, vel1, pos2, vel2;

            _splines.push_back(std::make_shared<_Spline>(spline_knots, t));
        }

        void add_intermediate_swing_knots(const Vector& knot_points, const Vector& times, size_t k_idx, size_t t_idx)
        {
            Time t = times[t_idx] / static_cast<double>(_num_knots_per_swing + 1);

            VecD pos1 = knot_points.segment(k_idx, _Spline::Dim);
            VecD vel1 = knot_points.segment(k_idx + 1 * _Spline::Dim, _Spline::Dim);
            VecD pos2 = knot_points.segment(k_idx + 2 * _Spline::Dim, _Spline::Dim);
            VecD vel2 = knot_points.segment(k_idx + 3 * _Spline::Dim, _Spline::Dim);

            Vector spline_knots((_Spline::Order + 1) * _dim);
            spline_knots << pos1, vel1, pos2, vel2;

            _splines.push_back(std::make_shared<_Spline>(spline_knots, t));
        }

        void add_final_swing_knots(const Vector& knot_points, const Vector& times, size_t k_idx, size_t t_idx)
        {
            Time t = times[t_idx] / static_cast<double>(_num_knots_per_swing + 1);

            VecD pos1 = knot_points.segment(k_idx, _Spline::Dim);
            VecD vel1 = knot_points.segment(k_idx + _Spline::Dim, _Spline::Dim);
            VecD pos2 = knot_points.segment(k_idx + 2 * _Spline::Dim, _Spline::Dim);
            VecD vel2 = VecD::Zero();

            Vector spline_knots((_Spline::Order + 1) * _dim);
            spline_knots << pos1, vel1, pos2, vel2;

            _splines.push_back(std::make_shared<_Spline>(spline_knots, t));
        }

    protected:
        static constexpr double _epsilon = 1e-12;
        const SplineType _spline_type;
        const size_t _dim;
        const size_t _num_vars_total{0};
        double _total_duration;

        const std::vector<Phase> _phase_sequence;
        const size_t _num_knots_per_swing;

        std::vector<SplinePtr> _splines;
    };

} // namespace rspl
