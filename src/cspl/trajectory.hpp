#include <Eigen/Dense>
#include <vector>

namespace cspl {
    template <unsigned int D>
    class Trajectory {
    public:
        using Vec = typename Eigen::Matrix<double, D, 1>;

        Trajectory(const CubicHermitePolynomial<D>& initial_pol, double duration)
        {
            PolynomialTimePair pair{initial_pol, duration};
            _vec.push_back(pair);
            _total_duration = duration;
        }

        void add_point(const Vec& next_pos, const Vec& next_vel, double duration)
        {
            auto prev_pos = _vec.back().pol.position(1); // REVIEW: Maybe store initial and last positions in polynomial object to avoid calling position(1)?
            auto prev_vel = _vec.back().pol.velocity(1);

            auto pol = CubicHermitePolynomial<D>(prev_pos, prev_vel, next_pos, next_vel);
            PolynomialTimePair pair{pol, duration};

            _vec.push_back(pair);
            _total_duration += duration;
        }

        Vec position(double t)
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _vec.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _vec[i].duration;
                if (t <= sum) {
                    double t_norm = (t - prev_sum) / _vec[i].duration;
                    return _vec[i].pol.position(t_norm);
                }
                prev_sum = sum;
            }
            // If t > trajectory duration, return final trajectory point.
            return _vec.back().pol.position(1);
        }

        Vec velocity(double t)
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _vec.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _vec[i].duration;
                if (t <= sum) {
                    double t_norm = (t - prev_sum) / _vec[i].duration;
                    return _vec[i].pol.velocity(t_norm);
                }
                prev_sum = sum;
            }
            // If t > trajectory duration, return final trajectory point.
            return _vec.back().pol.velocity(1);
        }

        Vec acceleration(double t)
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _vec.size(); i++) {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _vec[i].duration;
                if (t <= sum) {
                    double t_norm = (t - prev_sum) / _vec[i].duration;
                    return _vec[i].pol.acceleration(t_norm);
                }
                prev_sum = sum;
            }
            // If t > trajectory duration, return final trajectory point.
            return _vec.back().pol.acceleration(1);
        }

    protected:
        struct PolynomialTimePair {
            CubicHermitePolynomial<D> pol;
            double duration;
        };

        std::vector<PolynomialTimePair> _vec;
        double _total_duration;
    };
} // namespace cspl
