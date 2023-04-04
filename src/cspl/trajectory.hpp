#include <Eigen/Dense>

#include <iostream>
#include <vector>

namespace cspl {
    template <unsigned int D>
    class Trajectory {
    public:
        using Vec = typename Eigen::Matrix<double, D, 1>;

        struct PolynomialTimePair {
            CubicHermitePolynomial<D> polynomial;
            double duration;
        };

        Trajectory() : _total_duration(-1.) {}

        double total_duration() const { return _total_duration; }

        void add_point(const Vec& next_pos, const Vec& next_vel, double duration = 0.)
        {
            if (_total_duration < 0.) {
                _total_duration = 0.;
                _last_pos = next_pos;
                _last_vel = next_vel;

                if (duration > 0.) {
                    std::cerr << "You cannot have a duration > 0. for the initial point!" << std::endl;
                }
                return;
            }

            _polynomials.push_back({{_last_pos, _last_vel, next_pos, next_vel}, duration});
            _total_duration += duration;

            _last_pos = next_pos;
            _last_vel = next_vel;
        }

        void add_point(const Vec& next_pos, double duration = 0.)
        {
            if (_total_duration < 0.) {
                std::cerr << "You cannot use the acceleration constructor for the first point!" << std::endl;
                return;
            }

            // use the acceleration constructor
            _polynomials.push_back({{_last_pos, _last_vel, _polynomials.back().polynomial.acceleration(1.), next_pos, false}, duration});
            _total_duration += duration;

            _last_pos = next_pos;
            _last_vel = _polynomials.back().polynomial.velocity(1.);
        }

        Vec position(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _polynomials.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _polynomials[i].duration;
                if (t <= (sum + _epsilon)) {
                    double t_norm = (t - prev_sum) / _polynomials[i].duration;
                    return _polynomials[i].polynomial.position(t_norm);
                }
                prev_sum = sum;
            }

            // If t > trajectory duration, return final trajectory point.
            return _polynomials.back().polynomial.position(1.);
        }

        Vec velocity(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _polynomials.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _polynomials[i].duration;
                if (t <= (sum + _epsilon)) {
                    double t_norm = (t - prev_sum) / _polynomials[i].duration;
                    return _polynomials[i].polynomial.velocity(t_norm);
                }
                prev_sum = sum;
            }

            // If t > trajectory duration, return final trajectory point.
            return _polynomials.back().polynomial.velocity(1.);
        }

        Vec acceleration(double t) const
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _polynomials.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _polynomials[i].duration;
                if (t <= (sum + _epsilon)) {
                    double t_norm = (t - prev_sum) / _polynomials[i].duration;
                    return _polynomials[i].polynomial.acceleration(t_norm);
                }
                prev_sum = sum;
            }

            // If t > trajectory duration, return final trajectory point.
            return _polynomials.back().polynomial.acceleration(1.);
        }

        const std::vector<PolynomialTimePair>& polynomials() const { return _polynomials; }
        std::vector<PolynomialTimePair>& polynomials() { return _polynomials; }

    protected:
        static constexpr double _epsilon = 1e-12;
        double _total_duration;
        Vec _last_pos, _last_vel;
        std::vector<PolynomialTimePair> _polynomials;
    };

    using Trajectory2D = Trajectory<2>;
    using Trajectory3D = Trajectory<3>;
} // namespace cspl
