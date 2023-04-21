#include "splined_path.h"

using namespace splined_path;

//-- Constructors -------------------------------------------------------------
SplinedPath::SplinedPath(double mass_1, double mass_2, double mass_3,
            double mass_4, double cutoff)
:
s_I{std::pow(mass_1+mass_2,2.)},
s_III{std::pow(mass_4-mass_3,2.)}, s_IV{std::pow(mass_4+mass_3,2.)},
s_0{(std::pow(mass_1,2)+std::pow(mass_2,2)+std::pow(mass_3,2)+std::pow(mass_4,2))/3.},
egg(mass_4,mass_1)
{}

SplinedPathEta3Pi::SplinedPathEta3Pi(double mass_1, double mass_4, double cutoff)
:
SplinedPath{mass_1,mass_1,mass_1,mass_4,cutoff},
s_list(mass_1,mass_4, cutoff)
{build_egg(); spline_path();}

SplinedPathEtap3Pi::SplinedPathEtap3Pi(double mass_1, double mass_4, double cutoff)
:
SplinedPath{mass_1,mass_1,mass_1,mass_4,cutoff},
s_list(mass_1,mass_4, cutoff)
{build_egg(); spline_path();}

SplinedPathV3Pi::SplinedPathV3Pi(double mass_1, double mass_4, double cutoff)
:
SplinedPath{mass_1,mass_1,mass_1,mass_4,cutoff},
s_list(mass_1,mass_4, cutoff)
{build_egg(); spline_path();}

SplinedPathX3Pi::SplinedPathX3Pi(double mass_1, double mass_4, double cutoff)
:
SplinedPath{mass_1,mass_1,mass_1,mass_4,cutoff},
s_list(mass_1,mass_4, cutoff)
{build_egg(); spline_path();}

SplinedPathT3Pi::SplinedPathT3Pi(double mass_1, double mass_4, double cutoff)
:
SplinedPath{mass_1,mass_1,mass_1,mass_4,cutoff},
s_list(mass_1,mass_4, cutoff)
{build_egg(); spline_path();}

SplinedPathDerivativeEta3Pi::SplinedPathDerivativeEta3Pi(double mass_1, double mass_4, double cutoff)
:
s_list(mass_1,mass_4, cutoff),
egg(mass_1,mass_4, cutoff)
{build_derivative(); spline_derivative();}

SplinedPathDerivativeEtap3Pi::SplinedPathDerivativeEtap3Pi(double mass_1, double mass_4, double cutoff)
:
s_list(mass_1,mass_4, cutoff),
egg(mass_1,mass_4, cutoff)
{build_derivative(); spline_derivative();}

SplinedPathDerivativeV3Pi::SplinedPathDerivativeV3Pi(double mass_1, double mass_4, double cutoff)
:
s_list(mass_1,mass_4, cutoff),
egg(mass_1,mass_4, cutoff)
{build_derivative(); spline_derivative();}

SplinedPathDerivativeX3Pi::SplinedPathDerivativeX3Pi(double mass_1, double mass_4, double cutoff)
:
s_list(mass_1,mass_4, cutoff),
egg(mass_1,mass_4, cutoff)
{build_derivative(); spline_derivative();}

SplinedPathDerivativeT3Pi::SplinedPathDerivativeT3Pi(double mass_1, double mass_4, double cutoff)
:
s_list(mass_1,mass_4, cutoff),
egg(mass_1,mass_4, cutoff)
{build_derivative(); spline_derivative();}

//------------------------------------------------------------------------------

void SplinedPathEta3Pi::build_egg(){
    egg_list.reserve(s_list.interval_phi.size());
    for (std::size_t i=0; i<s_list.interval_phi.size(); ++i) {
        egg_list.push_back(egg(s_list.interval_phi[i]));
    }
    phi_lower_list.reserve(s_list.interval_egg.size());
    phi_upper_list.reserve(s_list.interval_egg.size());
    for (std::size_t i=0; i<s_list.interval_egg.size(); ++i) {
        phi_lower_list.push_back(egg.phi(s_list.interval_egg[i]).lower);
        phi_upper_list.push_back(egg.phi(s_list.interval_egg[i]).upper);
    }
    return;
}

void SplinedPathEtap3Pi::build_egg(){
    egg_list.reserve(s_list.interval_phi.size());
    for (std::size_t i=0; i<s_list.interval_phi.size(); ++i) {
        egg_list.push_back(egg(s_list.interval_phi[i]));
    }
    phi_lower_list.reserve(s_list.interval_egg.size());
    phi_upper_list.reserve(s_list.interval_egg.size());
    for (std::size_t i=0; i<s_list.interval_egg.size(); ++i) {
        phi_lower_list.push_back(egg.phi(s_list.interval_egg[i]).lower);
        phi_upper_list.push_back(egg.phi(s_list.interval_egg[i]).upper);
    }
    return;
}

void SplinedPathV3Pi::build_egg(){
    egg_list.reserve(s_list.interval_phi.size());
    for (std::size_t i=0; i<s_list.interval_phi.size(); ++i) {
        egg_list.push_back(egg(s_list.interval_phi[i]));
    }
    phi_lower_list.reserve(s_list.interval_egg.size());
    phi_upper_list.reserve(s_list.interval_egg.size());
    for (std::size_t i=0; i<s_list.interval_egg.size(); ++i) {
        phi_lower_list.push_back(egg.phi(s_list.interval_egg[i]).lower);
        phi_upper_list.push_back(egg.phi(s_list.interval_egg[i]).upper);
    }
    return;
}

void SplinedPathX3Pi::build_egg(){
    egg_list.reserve(s_list.interval_phi.size());
    for (std::size_t i=0; i<s_list.interval_phi.size(); ++i) {
        egg_list.push_back(egg(s_list.interval_phi[i]));
    }
    phi_lower_list.reserve(s_list.interval_egg.size());
    phi_upper_list.reserve(s_list.interval_egg.size());
    for (std::size_t i=0; i<s_list.interval_egg.size(); ++i) {
        phi_lower_list.push_back(egg.phi(s_list.interval_egg[i]).lower);
        phi_upper_list.push_back(egg.phi(s_list.interval_egg[i]).upper);
    }
    return;
}

void SplinedPathT3Pi::build_egg(){
    egg_list.reserve(s_list.interval_phi.size());
    for (std::size_t i=0; i<s_list.interval_phi.size(); ++i) {
        egg_list.push_back(egg(s_list.interval_phi[i]));
    }
    phi_lower_list.reserve(s_list.interval_egg.size());
    phi_upper_list.reserve(s_list.interval_egg.size());
    for (std::size_t i=0; i<s_list.interval_egg.size(); ++i) {
        phi_lower_list.push_back(egg.phi(s_list.interval_egg[i]).lower);
        phi_upper_list.push_back(egg.phi(s_list.interval_egg[i]).upper);
    }
    return;
}

void SplinedPathEta3Pi::spline_path(){
    spline_egg = cauchy::Interpolate(s_list.interval_phi, egg_list, gsl::InterpolationMethod::cubic);
    spline_phi_lower = gsl::Interpolate(s_list.interval_egg, phi_lower_list, gsl::InterpolationMethod::cubic);
    spline_phi_upper = gsl::Interpolate(s_list.interval_egg, phi_upper_list, gsl::InterpolationMethod::cubic);
    return;
}

void SplinedPathEtap3Pi::spline_path(){
    spline_egg = cauchy::Interpolate(s_list.interval_phi, egg_list, gsl::InterpolationMethod::cubic);
    spline_phi_lower = gsl::Interpolate(s_list.interval_egg, phi_lower_list, gsl::InterpolationMethod::cubic);
    spline_phi_upper = gsl::Interpolate(s_list.interval_egg, phi_upper_list, gsl::InterpolationMethod::cubic);
    return;
}

void SplinedPathV3Pi::spline_path(){
    spline_egg = cauchy::Interpolate(s_list.interval_phi, egg_list, gsl::InterpolationMethod::cubic);
    spline_phi_lower = gsl::Interpolate(s_list.interval_egg, phi_lower_list, gsl::InterpolationMethod::cubic);
    spline_phi_upper = gsl::Interpolate(s_list.interval_egg, phi_upper_list, gsl::InterpolationMethod::cubic);
    return;
}

void SplinedPathX3Pi::spline_path(){
    spline_egg = cauchy::Interpolate(s_list.interval_phi, egg_list, gsl::InterpolationMethod::cubic);
    spline_phi_lower = gsl::Interpolate(s_list.interval_egg, phi_lower_list, gsl::InterpolationMethod::cubic);
    spline_phi_upper = gsl::Interpolate(s_list.interval_egg, phi_upper_list, gsl::InterpolationMethod::cubic);
    return;
}

void SplinedPathT3Pi::spline_path(){
    spline_egg = cauchy::Interpolate(s_list.interval_phi, egg_list, gsl::InterpolationMethod::cubic);
    spline_phi_lower = gsl::Interpolate(s_list.interval_egg, phi_lower_list, gsl::InterpolationMethod::cubic);
    spline_phi_upper = gsl::Interpolate(s_list.interval_egg, phi_upper_list, gsl::InterpolationMethod::cubic);
    return;
}

void SplinedPathDerivativeEta3Pi::build_derivative(){
    egg_derivative.reserve(s_list.interval_phi.size());
    for (std::size_t i=0; i<s_list.interval_phi.size(); ++i) {
        egg_derivative.push_back(cauchy::derivative(egg.spline_egg, s_list.interval_phi[i],1e-8));
    }
    return;
}

void SplinedPathDerivativeEtap3Pi::build_derivative(){
    egg_derivative.reserve(s_list.interval_phi.size());
    for (std::size_t i=0; i<s_list.interval_phi.size(); ++i) {
        egg_derivative.push_back(cauchy::derivative(egg.spline_egg, s_list.interval_phi[i],1e-8));
    }
    return;
}

void SplinedPathDerivativeV3Pi::build_derivative(){
    egg_derivative.reserve(s_list.interval_phi.size());
    for (std::size_t i=0; i<s_list.interval_phi.size(); ++i) {
        egg_derivative.push_back(cauchy::derivative(egg.spline_egg, s_list.interval_phi[i],1e-8));
    }
    return;
}

void SplinedPathDerivativeX3Pi::build_derivative(){
    egg_derivative.reserve(s_list.interval_phi.size());
    for (std::size_t i=0; i<s_list.interval_phi.size(); ++i) {
        egg_derivative.push_back(cauchy::derivative(egg.spline_egg, s_list.interval_phi[i],1e-8));
    }
    return;
}

void SplinedPathDerivativeT3Pi::build_derivative(){
    egg_derivative.reserve(s_list.interval_phi.size());
    for (std::size_t i=0; i<s_list.interval_phi.size(); ++i) {
        egg_derivative.push_back(cauchy::derivative(egg.spline_egg, s_list.interval_phi[i],1e-8));
    }
    return;
}

void SplinedPathDerivativeEta3Pi::spline_derivative(){
    spline_egg_derivative = cauchy::Interpolate(s_list.interval_phi,egg_derivative,gsl::InterpolationMethod::cubic);
}

void SplinedPathDerivativeEtap3Pi::spline_derivative(){
    spline_egg_derivative = cauchy::Interpolate(s_list.interval_phi,egg_derivative,gsl::InterpolationMethod::cubic);
}

void SplinedPathDerivativeV3Pi::spline_derivative(){
    spline_egg_derivative = cauchy::Interpolate(s_list.interval_phi,egg_derivative,gsl::InterpolationMethod::cubic);
}

void SplinedPathDerivativeX3Pi::spline_derivative(){
    spline_egg_derivative = cauchy::Interpolate(s_list.interval_phi,egg_derivative,gsl::InterpolationMethod::cubic);
}

void SplinedPathDerivativeT3Pi::spline_derivative(){
    spline_egg_derivative = cauchy::Interpolate(s_list.interval_phi,egg_derivative,gsl::InterpolationMethod::cubic);
}


Complex SplinedPath::kappa(double s)const{
    return
     std::sqrt(Complex(1.-s_I/s))
    *std::sqrt(Complex(s_III-s))
    *std::sqrt(Complex(s_IV-s));
}

double SplinedPath::s_upper(double s)const{
    if (s_I<=s && s<s_III) {
        return 0.5*(3.*s_0-s+std::abs(kappa(s)));
    }
    if (s_III<=s && s<s_IV) {
        throw std::domain_error{"When integrating along the complex contour,"
            " the integration limits from path::PolarEgg should be used"};
    }
    if (s_IV<=s) {
        return 0.5*(3.*s_0-s-std::abs(kappa(s)));
    }
    throw std::domain_error{"Value for Mandelstam s is not allowed"};
}

double SplinedPath::s_lower(double s)const{
    if (s_I<=s && s<s_III) {
        return 0.5*(3*s_0-s-std::abs(kappa(s)));
    }
    if (s_III<=s && s<s_IV) {
        throw std::domain_error{"When integrating along the complex contour,"
            " the integration limits from path::PolarEgg should be used"};
    }
    if (s_IV<=s) {
        return 0.5*(3*s_0-s+std::abs(kappa(s)));
    }
    throw std::domain_error{"Value for Mandelstam s is not allowed"};
}
