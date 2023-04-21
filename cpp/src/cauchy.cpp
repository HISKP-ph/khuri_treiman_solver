#include "cauchy.h"

namespace cauchy {
// -- Basic facilities --------------------------------------------------------

// Template argument deduction in calls like `std::transform(.,.,func)`
// does not work if e.g. `func==std::real<Complex::value_type>`, for this
// purpose `real_specified`  and `imag_specified` are provided.

Complex::value_type real_specified(const Complex& z)
{
    return z.real();
}

Complex::value_type imag_specified(const Complex& z)
{
    return z.imag();
}

std::vector<double> real(const std::vector<Complex>& vec)
{
    std::vector<double> real_parts(vec.size());
    std::transform(vec.begin(),vec.end(),real_parts.begin(),
            real_specified);
    return real_parts;
}

std::vector<double> imag(const std::vector<Complex>& vec)
{
    std::vector<double> imaginary_parts(vec.size());
    std::transform(vec.begin(),vec.end(),imaginary_parts.begin(),
            imag_specified);
    return imaginary_parts;
}

// -- Integration -------------------------------------------------------------

std::tuple<Complex,double,double> c_integrate(const Curve& c,
        double lower, double upper, const gsl::Integration& integrate)
{
    gsl::Value real_part{
            integrate(facilities::compose(real_specified,c),lower,upper)};
    gsl::Value imaginary_part{
            integrate(facilities::compose(imag_specified,c),lower,upper)};
    Complex result{real_part.first,imaginary_part.first};
    return std::make_tuple(result,real_part.second,imaginary_part.second);
}

std::tuple<Complex,double,double> c_integrate(const Complex_function& f,
        const Curve& c, const Curve& c_derivative, double lower, double upper,
        const gsl::Integration& integrate)
{
    return c_integrate(
            [&](double x)
            {
                return f(c(x))*c_derivative(x);
            },
            lower,upper,integrate);
}

    
Complex complex_integration(const Curve& f, double lower, double upper, bool adaptive, const std::size_t integration_nodes){
    if (adaptive==false) {
        const std::size_t size{integration_nodes};
        gsl::GaussLegendre g{size};
        Complex value{0.+0.i};
        for (std::size_t i=0; i<size; ++i) {
            auto point{g.point(lower,upper,i)};
            value=value+point.second*f(point.first);
        }
        return value;
    }
    const auto integrate{gsl::Cquad{}};
    const auto result{c_integrate(f, lower, upper, integrate)};
    const auto value{std::get<0>(result)};
    return value;
}
    
// -- Interpolation -----------------------------------------------------------

Interpolate::Interpolate(const Interval& x,
        const std::vector<Complex>& y, gsl::InterpolationMethod m)
    : real_part{x,real(y),m},
    imaginary_part{x,imag(y),m}
{}

Complex Interpolate::operator()(double x) const
{
    return Complex{real_part(x),imaginary_part(x)};
}

Interpolate sample(const Complex_function& f, const Curve& c,
        const Interval& i, gsl::InterpolationMethod m)
{
    std::vector<Complex> y_values(i.size());
    std::transform(i.cbegin(),i.cend(),y_values.begin(),
            facilities::compose(f,c));

    return Interpolate{i,y_values,m};
}

Interpolate sample(const Curve& c, const Interval& i,
        gsl::InterpolationMethod m)
{
    return sample(facilities::identity<Complex>,c,i,m);
}

// -- Differentiation  --------------------------------------------------------

Complex derivative(const Curve& c, double value, double step_size,
        gsl::DerivativeMethod method)
{
    const auto real{gsl::derivative([&c](double x){return c(x).real();},
                                      value, step_size, method)};
    const auto imag{gsl::derivative([&c](double x){return c(x).imag();},
                                      value, step_size, method)};
    return Complex{std::get<0>(real), std::get<0>(imag)};
}

} // cauchy
