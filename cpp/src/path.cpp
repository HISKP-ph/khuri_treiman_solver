#include "path.h"

namespace path {
using facilities::square;

Complex Path::derivative(double curve_parameter, double step_size) const
{
    return cauchy::derivative([this](double x){return (*this)(x);},
                              curve_parameter, step_size);
}

Pinocchio::Pinocchio(double decay_mass, double pion_mass, double cut,
        double epsilon)
    : decay_mass{decay_mass}, pion_mass{pion_mass}, cut{cut},
    epsilon{epsilon}, egg{PolarEgg{decay_mass, pion_mass}}
{
}


Complex Pinocchio::operator()(double curve_parameter) const
{
    if (!domain().contains(curve_parameter)) {
        throw ParametrizationError{"attempt to evaluate curve outside domain"
                                   " of definition"};
    }
    const auto boundaries{parameter_boundaries()};
    if (curve_parameter < boundaries[1]) {
        return mandelstam::t_photon_pion_max(-curve_parameter, pion_mass,
                                             square(decay_mass));
    } else if (curve_parameter < boundaries[2]) {
        const auto phi{egg.phi(-curve_parameter)};
        return egg(phi.lower);
    } else if (curve_parameter < boundaries[3]) {
        return mandelstam::t_photon_pion_max(-curve_parameter, pion_mass,
                                             square(decay_mass))
            - Complex{0.0, epsilon};
    } else if (curve_parameter < boundaries[4]) {
        return mandelstam::t_photon_pion_max(-curve_parameter, pion_mass,
                                             square(decay_mass))
            + Complex{0.0, epsilon};
    } else if (curve_parameter < boundaries[5]) {
        return mandelstam::t_photon_pion_min(curve_parameter + shift(),
                                             pion_mass, square(decay_mass))
            + Complex{0.0, epsilon};
    } else if (curve_parameter < boundaries[6]) {
        const auto phi{egg.phi(curve_parameter + shift())};
        return egg(phi.upper);
    } else {
        return mandelstam::t_photon_pion_min(curve_parameter + shift(),
                                             pion_mass, square(decay_mass));
    }
}

Interval Pinocchio::curve_parameter(double mandelstam_s) const
{
    return {lower(mandelstam_s), upper(mandelstam_s)};
}

std::vector<double> Pinocchio::parameter_boundaries() const
{
    return {-cut, -four(), -three(), -two(), -one(), three() - shift(),
            four() - shift(), cut};
}

std::vector<double> Pinocchio::mandelstam_boundaries() const
{
    return {cut, four(), three(), two(), one(), three(), four(), cut};
}

Interval Pinocchio::domain() const
{
    return {-cut, cut};
}

double Pinocchio::lower(double mandelstam_s) const
{
    if (mandelstam_s > cut) {
        throw ParametrizationError{"attempt to evaluate curve above cutoff"};
    }
    return -mandelstam_s;
}

double Pinocchio::upper(double mandelstam_s) const
{
    if (mandelstam_s > cut) {
        throw ParametrizationError{"attempt to evaluate curve above cutoff"};
    }
    return mandelstam_s - shift();
}

PolarEgg::PolarEgg(double decay_mass, double pion_mass)
    : decay_mass{decay_mass}, pion_mass{pion_mass}
{
    if (decay_mass < 3.0 * pion_mass) {
        throw ParametrizationError{"parametrization cannot be used outside"
                                   " of the decay region"};
    }
}

Complex PolarEgg::operator()(double phi) const
{
    return radius(std::cos(phi)) * std::exp(Complex{0.0, phi});
}

double PolarEgg::radius(double cosine) const
{
    if (cosine < -1) {
        throw ParametrizationError{"cosine cannot be smaller than -1"};
    }
    else if (cosine < 0) {
        return radius_helper(cosine, 0);
    }
    else if (cosine == 0) {
        return squared_diff() * pion_mass / std::sqrt(3.0 * sum());
    }
    else if (cosine <= 1) {
        return radius_helper(cosine, 1);
    }
    else {
        throw ParametrizationError{"cosine cannot be bigger than 1"};
    }
}

Interval PolarEgg::phi(double mandelstam_variable) const
{
    const double x{mandelstam_variable};
    const double b{(square(decay_mass) - square(pion_mass)) * pion_mass};
    const double cosine{(3*sum()-x)*std::sqrt(x)/(2*b)};
    const double lower{-1*std::acos(cosine)};
    const double upper{std::acos(cosine)};
    return Interval{2.*constants::pi()+lower, upper};
}

double PolarEgg::sum() const
{
    return square(pion_mass) + square(decay_mass) / 3.0;
}

double PolarEgg::squared_diff() const
{
    return square(decay_mass) - square(pion_mass);
}

double PolarEgg::radius_helper(double cosine, int flag) const
{
    const int sign{
        [flag]() {
            switch (flag) {
                case 0:
                    return -1;
                case 1:
                    return 1;
                default:
                    throw ParametrizationError{"unknown flag"};
            }
        }()
    };
    const double cos{std::cos(angle(cosine, flag))};
    return sum() / 2.0 / cosine * (1.0 + sign * 2.0 * cos);
}

double PolarEgg::angle(double cosine, int flag) const
{
    const double arg{2.0
        * square(cosine * pion_mass * squared_diff()) / std::pow(sum(), 3)
        - 1};
    return (std::acos(arg) + flag * constants::pi()) / 3.0;
}
} // path
