#include <vector>

namespace cspl
{
    template <unsigned int DIM>
    class Trajectory
    {
    public:
        Trajectory(){};

        void push_back(const CubicHermitePolynomial<DIM> &pol, double duration)
        {
            PolynomialTimePair pair{pol, duration};
            _vec.push_back(pair);
        };

    protected:
        struct PolynomialTimePair
        {
            CubicHermitePolynomial<DIM> pol; // REVIEW : Should these be _pol and _duration?
            double duration;
        };

        std::vector<PolynomialTimePair> _vec;
    };
} // namespace cspl
