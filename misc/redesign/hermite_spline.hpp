#pragma once

#include "types.hpp"
#include <cstdint>

namespace rspl {
    template <size_t _Dim>
    class HermiteSpline {
    public:
        using VecD = Eigen::Matrix<double, _Dim, 1>;

        virtual SplineType type() const = 0;
        virtual size_t knots_required() const = 0;

        virtual size_t dimension() const = 0;
        virtual Time duration() const = 0;

        // virtual Vector knot_points() const = 0;
        // virtual Vector knot_points_start() const = 0;
        // virtual Vector knot_points_end() const = 0;

        virtual VecD position(Time t) const = 0;
        virtual VecD velocity(Time t) const = 0;
        virtual VecD acceleration(Time t) const = 0;
        // virtual VecD jerk(Time t) const = 0;
        // virtual VecD snap(Time t) const = 0;

        virtual Jacobian jacobian_position(Time t) const = 0;
        virtual Jacobian jacobian_velocity(Time t) const = 0;
        virtual Jacobian jacobian_acceleration(Time t) const = 0;
        // virtual Jacobian jacobian_jerk(Time t) const = 0;
        // virtual Jacobian jacobian_snap(Time t) const = 0;
    };

} // namespace rspl