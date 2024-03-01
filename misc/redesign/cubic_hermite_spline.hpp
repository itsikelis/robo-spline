#include "hermite_spline.hpp"
#include <cstdint>

namespace rspl {

    template <size_t _D>
    class CubicHermiteSpline : public HermiteSpline<_D> {
    public:
        using VecD = typename rspl::HermiteSpline<_D>::VecD;

        CubicHermiteSpline(const Vector& knot_points, Time duration) : _T(duration)
        {
            // To-Do Assert knot points is 4 * _D in size.
            _p0 = knot_points.head(_D); // p0
            _v0 = knot_points.segment(_D, _D); // v0
            _p1 = knot_points.segment(2 * _D, _D); // p1
            _v1 = knot_points.tail(_D); // v1

            calc_coeffs_from_knot_points();
        }

        inline void calc_coeffs_from_knot_points()
        {
            const double o_T = 1. / _T;
            const double o_T2 = 1. / (_T * _T);
            const double o_T3 = 1. / (_T * _T * _T);

            _c0 = _p0;
            _c1 = _v0;
            _c2 = 3. * _p1 * o_T2 - 3. * _p0 * o_T2 - 2. * _v0 * o_T - _v1 * o_T;
            _c3 = -2. * _p1 * o_T3 + 2. * _p0 * o_T3 + _v0 * o_T2 + _v1 * o_T2;
        }

        inline SplineType type() const override { return _type; }

        inline size_t dimension() const override { return _D; }

        inline Time duration() const override { return _T; }

        inline VecD position(Time t) const override
        {
            const double t2 = t * t;
            const double t3 = t * t2;

            return _c0 + (_c1 * t) + (_c2 * t2) + (_c3 * t3);
        }

        inline VecD velocity(Time t) const override
        {
            return _c1 + (2. * _c2 * t) + (3. * _c3 * t * t);
        }

        inline VecD acceleration(Time t) const override
        {
            return (2. * _c2) + (6. * _c3 * t);
        }

        Jacobian jacobian_position(Time t) const override
        {
            Vector deriv = deriv_pos(t);

            Jacobian jac(_D, _D * 4);

            for (size_t i = 0; i < _D; ++i) {
                for (size_t j = 0; j < 4; ++j) {
                    jac.coeffRef(i, i + j * _D) = deriv[j];
                }
            }

            return jac;
        }

        Jacobian jacobian_velocity(Time t) const override
        {
            Vector deriv = deriv_vel(t);

            Jacobian jac(_D, _D * 4);

            for (size_t i = 0; i < _D; ++i) {
                for (size_t j = 0; j < 4; ++j) {
                    jac.coeffRef(i, i + j * _D) = deriv[j];
                }
            }

            return jac;
        }

        Jacobian jacobian_acceleration(Time t) const override
        {
            Vector deriv = deriv_acc(t);

            Jacobian jac(_D, _D * 4);

            for (size_t i = 0; i < _D; ++i) {
                for (size_t j = 0; j < 4; ++j) {
                    jac.coeffRef(i, i + j * _D) = deriv[j];
                }
            }

            return jac;
        }

    protected:
        inline Vector deriv_pos(Time t) const
        {
            Vector deriv = Vector::Zero(4);

            const double t2 = t * t;
            const double t3 = t * t2;

            const double o_T = 1. / _T;
            const double o_T2 = 1. / (_T * _T);
            const double o_T3 = 1. / (_T * _T * _T);

            // Position spline derivatives w.r.t. x0, v0, x1, v1.
            deriv[0] = 1. - 3. * t2 * o_T2 + 2. * t3 * o_T3;
            deriv[1] = t - 2. * t2 * o_T + t3 * o_T2;
            deriv[2] = 3. * t2 * o_T2 - 2. * t3 * o_T3;
            deriv[3] = -t2 * o_T + t3 * o_T2;

            return deriv;
        }

        inline Vector deriv_vel(Time t) const
        {
            const double t2 = t * t;
            Vector deriv = Vector::Zero(4);

            const double o_T = 1. / _T;
            const double o_T2 = 1. / (_T * _T);
            const double o_T3 = 1. / (_T * _T * _T);

            // Velocity spline derivatives w.r.t. x0, v0, x1, v1.
            deriv[0] = -6. * t * o_T2 + 6. * t2 * o_T3;
            deriv[1] = 1. - 4. * t * o_T + 3. * t2 * o_T2;
            deriv[2] = 6. * t * o_T2 - 6. * t2 * o_T3;
            deriv[3] = -2. * t * o_T + 3. * t2 * o_T2;

            return deriv;
        }

        inline Vector deriv_acc(Time t) const
        {
            Vector deriv = Vector::Zero(4);

            const double o_T = 1. / _T;
            const double o_T2 = 1. / (_T * _T);
            const double o_T3 = 1. / (_T * _T * _T);

            // Acceleration spline derivatives w.r.t. x0, v0, x1, v1.
            deriv[0] = -6. * o_T2 + 12. * t * o_T3;
            deriv[1] = -4. * o_T + 6. * t * o_T2;
            deriv[2] = 6. * o_T2 - 12. * t * o_T3;
            deriv[3] = -2. * o_T + 6. * t * o_T2;

            return deriv;
        }

    private:
        const SplineType _type{SplineType::CubicHermite};

        // Initial and final knot points.
        VecD _p0;
        VecD _v0;
        VecD _p1;
        VecD _v1;

        // Spline duration.
        Time _T;

        // Spline polynomial coefficients.
        VecD _c0;
        VecD _c1;
        VecD _c2;
        VecD _c3;
    };

} // namespace rspl