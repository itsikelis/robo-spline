#pragma once

#include "types.hpp"
#include <cstdint>

namespace rspl {
    template <size_t _Dim>
    class HermiteSpline {
    public:
        static const size_t Dim = _Dim;
        using VecD = Eigen::Matrix<double, _Dim, 1>;

        virtual size_t dim() const = 0;
        virtual Time duration() const = 0;

        // virtual Vector knot_points() const = 0;
        // virtual Vector knot_points_start() const = 0;
        // virtual Vector knot_points_end() const = 0;

        virtual VecD evaluate(Time t, size_t order) const = 0;

        virtual VecD pos(Time t) const = 0;
        virtual VecD vel(Time t) const = 0;
        virtual VecD acc(Time t) const = 0;
        // virtual VecD jerk(Time t) const = 0;
        // virtual VecD snap(Time t) const = 0;

        virtual Jacobian jac_block(Time t, size_t order) const = 0;

        virtual Jacobian jac_block_pos(Time t) const = 0;
        virtual Jacobian jac_block_vel(Time t) const = 0;
        virtual Jacobian jac_block_acc(Time t) const = 0;
        // virtual Jacobian jac_block_jerk(Time t) const = 0;
        // virtual Jacobian jac_block_snap(Time t) const = 0;
    };

} // namespace rspl