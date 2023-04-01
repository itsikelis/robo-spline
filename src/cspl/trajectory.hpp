#include <vector>
#include <Eigen/Dense>

namespace cspl
{
    template <unsigned int DIM>
    class Trajectory
    {
    public:
        using Vec = Eigen::Matrix<double, DIM, 1>;

        Trajectory(const CubicHermitePolynomial<DIM> &initial_pol, double duration)
        {
            PolynomialTimePair pair{initial_pol, duration};
            _vec.push_back(pair);
            _total_duration = duration;
        }

        void add_point(const Vec &next_pos, const Vec &next_vel, double duration)
        {
            auto prev_pos = _vec.back().pol.position(1); // REVIEW: Maybe store initial and last positions in polynomial object to avoid calling position(1)?
            auto prev_vel = _vec.back().pol.velocity(1);

            auto pol = CubicHermitePolynomial3D(prev_pos, prev_vel, next_pos, next_vel);
            PolynomialTimePair pair{pol, duration};

            _vec.push_back(pair);
            _total_duration += duration;
        }

        Vec position(double t)
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _vec.size(); i++)
            {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _vec[i].duration;
                if (t <= sum)
                {
                    double t_norm = (t - prev_sum) / _vec[i].duration;
                    return _vec[i].pol.position(t_norm);
                }
                prev_sum = sum;
            }
        }

        Vec velocity(double t)
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _vec.size(); i++)
            {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _vec[i].duration;
                if (t <= sum)
                {
                    double t_norm = (t - prev_sum) / _vec[i].duration;
                    return _vec[i].pol.velocity(t_norm);
                }
                prev_sum = sum;
            }
        }

        Vec acceleration(double t)
        {
            double sum = 0;
            double prev_sum = 0;
            for (unsigned int i = 0; i < _vec.size(); i++)
            {
                // Add current polynomial's duration to sum and check if t belongs to this interval.
                sum += _vec[i].duration;
                if (t <= sum)
                {
                    double t_norm = (t - prev_sum) / _vec[i].duration;
                    return _vec[i].pol.acceleration(t_norm);
                }
                prev_sum = sum;
            }
        }

    protected:
        struct PolynomialTimePair
        {
            CubicHermitePolynomial<DIM> pol;
            double duration;
        };

        std::vector<PolynomialTimePair> _vec;
        double _total_duration;
    };
} // namespace cspl
