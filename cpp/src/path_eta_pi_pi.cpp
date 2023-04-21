#include "path_eta_pi_pi.h"

using namespace enums;

namespace path_eta_pi_pi {


//-- Constructors ---------------------------------------------------------------------------------

Path::Path(double mass_1, double mass_2, double mass_3, double mass_4)
:
delta{(std::pow(mass_4,2.)-std::pow(mass_2,2.)) * (std::pow(mass_3,2.)-std::pow(mass_1,2.))},
mass_1{mass_1},mass_2{mass_2},mass_3{mass_3},mass_4{mass_4},
s_I{std::pow(mass_1+mass_2,2)},
s_III{std::pow(mass_4-mass_3,2)}, s_IV{std::pow(mass_4+mass_3,2)},
s_0{(std::pow(mass_1,2)+std::pow(mass_2,2)+std::pow(mass_3,2)+std::pow(mass_4,2))/3.},
t_I{std::pow(mass_3+mass_1,2)},
t_III{std::pow(mass_4-mass_2,2)}, t_IV{std::pow(mass_4+mass_2,2)},
t_0{(std::pow(mass_1,2)+std::pow(mass_2,2)+std::pow(mass_3,2)+std::pow(mass_4,2))/3.}
{}

//---------------------------------------------------------------------------------------------------


Complex Path::kallen(double s, double mass1, double mass2)const{
    return std::pow(s,2.) + std::pow(mass1,4.) + std::pow(mass2,4.)
    -2.*(s*std::pow(mass1,2.) + s*std::pow(mass2,2.) + std::pow(mass1,2.)*std::pow(mass2,2.));
}

double Path::kappa_pi(double s)const{
    return std::sqrt(std::abs(kallen(s, mass_1, mass_2)*kallen(s, mass_3, mass_4)))/s;
}

double Path::kappa_eta(double t)const{
    return std::sqrt(std::abs(kallen(t, mass_1, mass_3)*kallen(t, mass_2, mass_4)))/t;
}



double Path::integration_limit(double x, TypeOfAverage average,
                                PinocchioEndpoints endpoints)const{
    switch (average) {
        case TypeOfAverage::zero:
            if (s_I<=x && x<=s_III) {
                switch (endpoints) {
                    case PinocchioEndpoints::upper: return 0.5*(3.*s_0-x+kappa_pi(x)); break;
                    case PinocchioEndpoints::lower: return 0.5*(3.*s_0-x-kappa_pi(x)); break;
                }
            }
            if (s_III<x && x<s_IV) {
                throw std::domain_error{"When integrating along the complex contour,"
                                    " the integration limits from path_eta_pi_pi::contour should be used"};
            }
            if (s_IV<=x) {
                switch (endpoints) {
                    case PinocchioEndpoints::upper: return 0.5*(3.*s_0-x-kappa_pi(x)); break;
                    case PinocchioEndpoints::lower: return 0.5*(3.*s_0-x+kappa_pi(x)); break;
                }
            }
            throw std::domain_error{"Value for Mandelstam s is not allowed"};
            
        case TypeOfAverage::minus:
            if (t_I<=x && x<=t_III) {
                switch (endpoints) {
                    case PinocchioEndpoints::upper: return 0.5*(3.*t_0-x-delta/x +kappa_eta(x)); break;
                    case PinocchioEndpoints::lower: return 0.5*(3.*t_0-x-delta/x -kappa_eta(x)); break;
                }
            }
            if (t_III<x && x<t_IV) {
                throw std::domain_error{"When integrating along the complex contour,"
                                    " the integration limits from path_eta_pi_pi::contour should be used"};
            }
            if (t_IV<=x) {
                switch (endpoints) {
                    case PinocchioEndpoints::upper: return 0.5*(3.*t_0-x-delta/x -kappa_eta(x)); break;
                    case PinocchioEndpoints::lower: return 0.5*(3.*t_0-x-delta/x +kappa_eta(x)); break;
                }
            }
            throw std::domain_error{"Value for Mandelstam t is not allowed"};
            

        case TypeOfAverage::plus:
            if (t_I<=x && x<=t_III) {
                switch (endpoints) {
                    case PinocchioEndpoints::upper: return 0.5*(3.*t_0-x+delta/x +kappa_eta(x)); break;
                    case PinocchioEndpoints::lower: return 0.5*(3.*t_0-x+delta/x -kappa_eta(x)); break;
                }
            }
            if (t_III<x && x<t_IV) {
                throw std::domain_error{"When integrating along the complex contour,"
                                    " the integration limits from path_eta_pi_pi::contour should be used"};
            }
            if (t_IV<=x) {
                switch (endpoints) {
                    case PinocchioEndpoints::upper: return 0.5*(3.*t_0-x+delta/x -kappa_eta(x)); break;
                    case PinocchioEndpoints::lower: return 0.5*(3.*t_0-x+delta/x +kappa_eta(x)); break;
                }
            }
            throw std::domain_error{"Value for Mandelstam t is not allowed"};

        default:
            throw std::domain_error{"Average is not defined in Path::integration_limit"};
    }
}



// For integration region III, i.e. egg-like shape
// endpoints changes the sign of kappa
// average decides which angular average to use
Complex Path::s_plus_minus(double x, TypeOfAverage average,
                           PinocchioEndpoints endpoints)const{
    switch (average) {
        case TypeOfAverage::zero:
            switch (endpoints) {
                case PinocchioEndpoints::upper: return 0.5*(3.*s_0-x +1.i*kappa_pi(x)); break;
                case PinocchioEndpoints::lower: return 0.5*(3.*s_0-x -1.i*kappa_pi(x)); break;
                default: throw std::domain_error{"Endpoints is not defined"};
            }
            break;

        case TypeOfAverage::minus:
            switch (endpoints) {
                case PinocchioEndpoints::upper: return 0.5*(3.*t_0-x -delta/x +1.i*kappa_eta(x)); break;
                case PinocchioEndpoints::lower: return 0.5*(3.*t_0-x -delta/x -1.i*kappa_eta(x)); break;
                default: throw std::domain_error{"Endpoints is not defined"};
            }
            break;

        case TypeOfAverage::plus:
            switch (endpoints) {
                case PinocchioEndpoints::upper: return 0.5*(3.*t_0-x +delta/x +1.i*kappa_eta(x)); break;
                case PinocchioEndpoints::lower: return 0.5*(3.*t_0-x +delta/x -1.i*kappa_eta(x)); break;
                default: throw std::domain_error{"Endpoints is not defined"};
            }
            break;

        default:
            throw std::domain_error{"Average is not defined in Path::s_plus_minus"};
    }
}


double Path::t_to_yb(double yb)const{
    return -1.*std::pow(yb,2.)*(-2*mass_3-2*mass_4+std::pow(yb,2.));
}

double Path::t_to_ya(double ya)const{
    return
    std::pow(
             std::pow(mass_4-mass_3,2.) +std::pow(ya,4.),
             2.)
    /(4.*std::pow(ya,4.));
}


Complex Path::contour(double y, TypeOfAverage average,
                         PinocchioEndpoints endpoints, IntegrationRegion region)const{
    double t;
    switch (region) {
        case IntegrationRegion::a:
            t=t_to_ya(y);
            return s_plus_minus(t, average, endpoints);
            break;
        case IntegrationRegion::b:
            t=t_to_yb(y);
            return s_plus_minus(t, average, endpoints);
            break;
        default:
            throw std::domain_error{"Region is not defined in Path::contour"};
    }
}

// Write the contour in a more convenient form
Complex Path::operator()(double y, Setting evaluation)const{
    switch (evaluation) {
        case Setting::plus_upper_a:  return contour(y, TypeOfAverage::plus, PinocchioEndpoints::upper, IntegrationRegion::a);  break;
        case Setting::plus_upper_b:  return contour(y, TypeOfAverage::plus, PinocchioEndpoints::upper, IntegrationRegion::b);  break;
        case Setting::plus_lower_a:  return contour(y, TypeOfAverage::plus, PinocchioEndpoints::lower, IntegrationRegion::a);  break;
        case Setting::plus_lower_b:  return contour(y, TypeOfAverage::plus, PinocchioEndpoints::lower, IntegrationRegion::b);  break;
        case Setting::minus_upper_a: return contour(y, TypeOfAverage::minus, PinocchioEndpoints::upper, IntegrationRegion::a); break;
        case Setting::minus_upper_b: return contour(y, TypeOfAverage::minus, PinocchioEndpoints::upper, IntegrationRegion::b); break;
        case Setting::minus_lower_a: return contour(y, TypeOfAverage::minus, PinocchioEndpoints::lower, IntegrationRegion::a); break;
        case Setting::minus_lower_b: return contour(y, TypeOfAverage::minus, PinocchioEndpoints::lower, IntegrationRegion::b); break;
        case Setting::zero_upper_a:  return contour(y, TypeOfAverage::zero, PinocchioEndpoints::upper, IntegrationRegion::a);  break;
        case Setting::zero_upper_b:  return contour(y, TypeOfAverage::zero, PinocchioEndpoints::upper, IntegrationRegion::b);  break;
        case Setting::zero_lower_a:  return contour(y, TypeOfAverage::zero, PinocchioEndpoints::lower, IntegrationRegion::a);  break;
        case Setting::zero_lower_b:  return contour(y, TypeOfAverage::zero, PinocchioEndpoints::lower, IntegrationRegion::b);  break;
        default:
            throw std::domain_error{"The parameterization of the complex path in eta'-> eta pi pi does not need the demanded evaluation."};
            break;
    }
}

Complex Path::derivative(double y, Setting evaluation)const{
    // write path in terms of only one variable
    Setting eval=evaluation;
    const auto lambda{[eval, this](double x)
        { return this->operator()(x, eval); }};
    
    // this can now be differentiated easily
    return cauchy::derivative(lambda, y, 1e-8);
}


double Path::yb(double t)const{
    return
    std::sqrt(
              ((mass_4+mass_3)+std::sqrt(t))/2.
              )
    -
    std::sqrt(
              ((mass_4+mass_3)-std::sqrt(t))/2.
              );
}


double Path::ya(double t)const{
    return
    std::sqrt(
              (std::sqrt(t)+(mass_4-mass_3))/2.
              )
    -
    std::sqrt(
              (std::sqrt(t)-(mass_4-mass_3))/2.
              );
}


}// namespace path_eta_pi_pi
