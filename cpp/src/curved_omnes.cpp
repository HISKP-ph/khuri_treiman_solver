#include "curved_omnes.h"

namespace curved_omnes {
std::vector<Complex> first_points(const grid::Curve& curve, std::size_t size)
{
    const std::vector<double> boundaries{curve.boundaries()};
    if (boundaries.size() < size) {
        std::stringstream message;
        message << "Tried to retrieve " << size << " elements, but curve has "
            << "only " << boundaries.size() << " boundary points.";
        throw std::invalid_argument{message.str()};
    }
    std::vector<Complex> points;
    points.reserve(size);
    const auto start{boundaries.cbegin()};
    std::transform(start, start + size, std::back_inserter(points),
            [&curve](double x){return curve.curve_func(x);});
    return points;
}

std::vector<Complex> all_points(const grid::Curve& curve)
{
    return first_points(curve, curve.boundaries().size());
}

bool on_second_sheet(const std::vector<Complex>& points,
        const Complex& mandelstam_s)
{
    const auto size{points.size()};
    if (size == 2)
        return false;
    if (size > 3)
        return points[0].real() < mandelstam_s.real()
            && mandelstam_s.real() < points[3].real()
            && points[1].imag() < mandelstam_s.imag()
            && mandelstam_s.imag() < 0.0;
    throw std::runtime_error{"Don't know how to handle this case"};
}
} // curved_omnes
