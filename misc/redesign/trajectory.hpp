#include <cstdint>
#include <iostream>
#include <memory>

#include "hermite_spline.hpp"
#include "types.hpp"

namespace rspl {
    template <typename _SplineType, size_t _NumKnots>
    class Trajectory {
    public:
        using VecD = typename HermiteSpline<_Dim>::VecD;
        using SplinePtr = std::shared_ptr<_SplineType>;
        // Number of required knot points for spline (4 for cubic, 6 for quintic etc).
        // using NumSplineReqKnots = typename _SplineType<_Dim>::NumReqKnots;

    public:
        Trajectory() : _total_duration(-1.) {}

        Trajectory(const Vector& knot_points, const Vector& times)
        {
            // Todo:
            // assert knot_points and times size are correct.
            // assert NumKnots - 1 == times.size()

            // Create a temp spline object to get num of knots needed (4 for cubic, 6 for quintic etc).
            _SplineType temp;
            const size_t R = temp.knots_required();
            const size_t D = temp.dimension();

            const size_t iters = _NumKnots - 1;
            for (size_t i = 0; i < iters; ++i) {
                Eigen::VectorXd spline_knots = knot_points.segment(i * R * D, R * D);
                Time dt = times(i);

                _splines.push_back(std::make_shared<SplinePtr>(spline_knots, dt));
                _total_duration += dt;
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
            if (t > _total_duration) {
                // If t > _total_duration, return final time of last spline.
                return std::make_pair(static_cast<int>(_splines.size() - 1), _splines.back()->duration());
            }

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
        }

    protected:
        static constexpr double _epsilon = 1e-12;

        Time _total_duration;
        std::vector<SplinePtr> _splines;
    };

    // Aliases.
    using Trajectory2D = Trajectory<2>;
    using Trajectory3D = Trajectory<3>;

} // namespace rspl
