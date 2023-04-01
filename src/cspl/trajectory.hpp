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
        }

        void add_point(const Vec &next_pos, const Vec &next_vel, double duration)
        {
            auto prev_pos = _vec.back().pol.position(1); // REVIEW: Maybe store initial and last positions in polynomial object to avoid calling position(1)?
            auto prev_vel = _vec.back().pol.velocity(1);

            auto pol = CubicHermitePolynomial3D(prev_pos, prev_vel, next_pos, next_vel);
            PolynomialTimePair pair{pol, duration};

            _vec.push_back(pair);
        }

    protected:
        struct PolynomialTimePair
        {
            CubicHermitePolynomial<DIM> pol;
            double duration;
        };

        std::vector<PolynomialTimePair> _vec;
    };
} // namespace cspl
