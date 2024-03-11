#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

#include "hermite_spline.hpp"
#include "types.hpp"

namespace rspl {

    enum class Phase {
        Step,
        Swing,
    };

    template <typename _Spline, size_t _NumKnots>
    class PhasedTrajectory {
    public:
        using VecD = typename _Spline::VecD;
        using SplinePtr = std::shared_ptr<_Spline>;

        PhasedTrajectory() : _spline_type(_Spline::Type), _dim(_Spline::Dim), _total_duration(-1.) {}

        PhasedTrajectory(const Vector& knot_points, const Vector& times) : _spline_type(_Spline::Type), _dim(_Spline::Dim), _num_vars_total(knot_points.rows()), _total_duration(times.sum())
        {
            // TODO:
            // assert knot_points and times size are correct.
            // assert NumKnots - 1 == times.size()

            generate_trajectory_from_points(knot_points, times, _spline_type);
        }

        void generate_trajectory_from_points(const Vector& knot_points, const Vector& times, const std::vector<Phase>& phase_sequence, const Vector& phase_times)
        {
        }

        Phase phase_at_time(Time t) const
        {
            double sum = 0.;
            double prev_sum = 0.;
            size_t iters = _phase_sequence.size();
            for (size_t i = 0; i < iters; ++i) {
                sum += _phase_times[i];

                if (t <= sum - _epsilon) {
                    Time t_norm = (t - prev_sum);
                    return _phase_sequence[i];
                }

                prev_sum = sum;
            }

            // If t > _total_duration, return final phase.
            return _phase_sequence.back();
        }

    protected:
        static constexpr double _epsilon = 1e-12;
        const SplineType _spline_type;
        const size_t _dim;
        const size_t _num_vars_total{0};

        const Vector _phase_times;
        const std::vector<Phase> _phase_sequence;

        double _total_duration;
        std::vector<SplinePtr> _splines;
    };

} // namespace rspl
