#include <cstdint>

#include "hermite_spline.hpp"
#include "types.hpp"

namespace rspl {
    template <size_t _D>
    class Trajectory {
    public:
        using VecD = typename HermiteSpline<_D>::VecD;
        using SplineIndex = size_t;
        using SplinePtr = std::shared_ptr<CubicHermiteSpline<D>>;

    public:
        Trajectory() : _total_duration(-1.) {}

        Trajectory(const Vector& knot_points, const Vector& times)
        {
            // Todo: assert knot_points and times size are correct.
        }

        ~Trajectory();

        void clear()
        {
            _splines.clear();
            _total_duration = -1.;
        }

        double total_duration() const { return _total_duration; }

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

    protected:
    };

} // namespace rspl