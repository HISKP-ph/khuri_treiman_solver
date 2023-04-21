#include "piecewise.h"
#include <iterator>

namespace piecewise {
Piecewise::Piecewise(const std::vector<Complex>& knots,
        const std::vector<Para>& parametrisations)
    : parametrisations{parametrisations}
{
    const auto s{parametrisations.size()};
    if (s+1 != knots.size())
        throw std::invalid_argument{
            "Each curve segment needs one parametrisation."};

    std::vector<Complex> differences(knots.size());
    std::adjacent_difference(knots.cbegin(),knots.cend(),differences.begin());

    pieces.resize(s);
    adjacent.resize(s);

    for (std::size_t i{0}; i<s; ++i) {
        pieces[i] = std::make_pair(differences[i+1],knots[i]);
        adjacent[i] = std::make_pair(knots[i],knots[i+1]);
    }
}

Piecewise::Piecewise(const std::vector<Complex>& knots)
{
    const auto size{knots.size()};
    if (size < 2)
        throw std::invalid_argument{"Curve needs to have at least two knots."};
    *this = Piecewise{knots, Piecewise::all_linear(size - 1)};
}

std::size_t Piecewise::piece_index(double x) const
{
    if (x<lower() || x>upper())
        throw std::out_of_range{
            "Tried to evaluate piecewise curve outside domain of definition."};
    const auto index{static_cast<std::size_t>(x)};
    // class invariant assures that upper()>=1.
    return index==upper() ? index-1 : index;
}

Complex Piecewise::curve_func(double x) const
{
    const auto k{piece_index(x)};
    const auto& p{pieces.at(k)};
    switch (parametrisations[k]) {
        case linear:
            return p.first*(x-k) + p.second;
        case quadratic:
            return p.first*square(x-k) + p.second;
        default:
            throw UnknownPara{};
    }
}

Complex Piecewise::derivative_func(double x) const
{
    const auto k{piece_index(x)};
    const auto& p{pieces.at(k)};
    switch (parametrisations[k]) {
        case linear:
            return p.first;
        case quadratic:
            return 2.0*p.first*(x-k);
        default:
            throw UnknownPara{};
    }
}

bool in_between(const Complex& x, const Complex& a, const Complex& b)
    // Determine whether `x` is on the connection line of `a` and `b`
    // in between `a` and `b`.
{
    constexpr double minimal_distance{1e-10};
    const double difference{std::abs(x-a) + std::abs(x-b) - std::abs(a-b)};
    return std::abs(difference) < minimal_distance;
}

Piecewise::Segment Piecewise::hits(const Complex& s) const
{
    auto position{std::find_if(adjacent.cbegin(),adjacent.cend(),
            [s](const auto& p){return in_between(s,p.first,p.second);})};
    if (position==adjacent.cend())
        return std::nullopt;
    double lower = std::distance(adjacent.cbegin(),position);
    double upper{lower+1.0};
    return std::make_pair(lower,upper);
}

std::vector<double> Piecewise::boundaries() const
{
    std::vector<double> result(pieces.size()+1);
    std::iota(result.begin(),result.end(),lower());
    return result;
}

std::vector<Complex> real_points(double threshold, double cut,
                                const std::vector<double>& intermediate)
{
    std::vector<Complex> result;
    result.reserve(intermediate.size() + 2);
    result.push_back(threshold);
    std::copy(intermediate.cbegin(), intermediate.cend(),
              std::back_inserter(result));
    result.push_back(cut);
    return result;
}

Real::Real(double threshold, double cut,
           const std::vector<double>& intermediate)
    : Piecewise{real_points(threshold, cut, intermediate)}
{
}

std::vector<Complex> vector_decay_points(double pion_mass, double virtuality,
        double cut)
{
    const double m2{square(pion_mass)};
    const double a{virtuality-2.5*m2};
    const double b{-7.0*m2};

    const double x1{4.0*m2};
    const Complex x2{5.0*m2,b};
    const Complex x3{a,b};
    const double x4{a};
    const double x5{mandelstam::s_greater(pion_mass,virtuality)};

    return {x1,x2,x3,x4,x5,cut};
}

VectorDecay::VectorDecay(double pion_mass, double virtuality, double cut)
    : Piecewise{vector_decay_points(pion_mass,virtuality,cut)}
{
}

std::vector<Complex> adaptive_points(double pion_mass, double virtuality,
        double cut)
{
    const auto m2{square(pion_mass)};
    const mandelstam::Critical critical{pion_mass,virtuality};
    const auto lower{-critical.imaginary_radius()};
    const auto right{critical.right()+m2};
    const double threshold{4.0*m2};

    std::vector<Complex> result;
    result.reserve(4);
    result.emplace_back(threshold);
    result.emplace_back(threshold,lower);
    result.emplace_back(right,lower);
    result.emplace_back(right);
    if (cut <= result.back().real())
        return result;
    result.emplace_back(mandelstam::s_greater(pion_mass,virtuality));
    if (cut <= result.back().real())
        return result;
    result.emplace_back(cut);
    return result;
}

Adaptive::Adaptive(double pion_mass, double virtuality, double cut)
    : Piecewise{adaptive_points(pion_mass,virtuality,cut)}
{
}
} // piecewise
