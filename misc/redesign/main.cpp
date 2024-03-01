#include "cubic_hermite_spline.hpp"
#include "hermite_spline.hpp"
#include <iostream>

int main()
{
    Eigen::VectorXd vec(8);
    vec << 0., 0., 0., 0., 1., 1., 0., 0.;
    rspl::CubicHermiteSpline<2> spline(vec, 1.);
    std::cout << "Hello, World!" << std::endl;
    return 0;
}