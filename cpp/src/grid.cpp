#include "grid.h"

namespace grid {
std::vector<Complex> boundary_points(const Curve& c)
{
    const auto boundaries{c.boundaries()};
    std::vector<Complex> points;
    points.reserve(boundaries.size());
    for (const auto& b: boundaries)
        points.emplace_back(c.curve_func(b));
    return points;
}

Knots generate_knots(double start, double end, std::size_t points)
{
    const gsl::GaussLegendre g{points};
    Knots path(points);
    for (std::size_t i{0}; i<points; ++i)
        path[i] = g.point(start,end,i);
    return path;
}
} // grid
