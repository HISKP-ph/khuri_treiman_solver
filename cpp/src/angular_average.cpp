#include "angular_average.h"
#include "matching.h"
#include "array.h"
#include "path.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <complex>

using namespace angular_average;


//-- Constructors ---------------------------------------------------------------------------------

AngularAverage::AngularAverage(double mass_1, double mass_2, double mass_3, double mass_4,
                                        FunctionSet amplitude, int isospin,
                                        double cutoff, int iteration, double epsilon, double validity,
                                        SubtractionConstant sub)
:
mass_1{mass_1},mass_2{mass_2},mass_3{mass_3},mass_4{mass_4},
s_I{std::pow(mass_1+mass_2,2)},s_II{mass_1*(std::pow(mass_4,2.)-std::pow(mass_3,2.))/(mass_3+mass_2)},
s_III{std::pow(mass_4-mass_3,2)}, s_IV{std::pow(mass_4+mass_3,2)},
s_0{(std::pow(mass_1,2)+std::pow(mass_2,2)+std::pow(mass_3,2)+std::pow(mass_4,2))/3.},
t_I{std::pow(mass_3+mass_1,2)},
t_II{(std::pow(mass_4,2.)+std::pow(mass_3,2.)-2.*std::pow(mass_1,2.))/2.},
t_III{std::pow(mass_4-mass_1,2)}, t_IV{std::pow(mass_4+mass_1,2)},
t_0{(std::pow(mass_1,2)+std::pow(mass_2,2)+std::pow(mass_3,2)+std::pow(mass_4,2))/3.},
u_II{(  mass_3*(std::pow(mass_4,2.)-std::pow(mass_1,2.))
      -mass_1*(std::pow(mass_3,2.)-std::pow(mass_2,2.)))/(mass_3+mass_2)},
s_x{std::pow(mass_1-mass_2,2.)},
t_x{std::pow(mass_1-mass_3,2.)},
validity{validity}, epsilon{epsilon},
isospin{isospin}, iteration{iteration},
sub{sub}, amplitude{amplitude},
matched(mass_1,mass_1,mass_1,mass_4, epsilon, validity)
{}

AngularAverageEta3Pi::AngularAverageEta3Pi(double mass_1, double mass_4,
                                        FunctionSet amplitude, int isospin,
                                        double cutoff, int iteration, double epsilon, double validity,
                                        SubtractionConstant sub)
:
AngularAverage{mass_1,mass_1,mass_1,mass_4,amplitude,isospin,cutoff,iteration,epsilon,validity,sub},
mandelstam_list(mass_1, mass_4, cutoff),
polar_egg(mass_1, mass_4, cutoff),
egg_deriv(mass_1, mass_4, cutoff)
{build_tilde_functions(); build_matching_list(); make_splines(); build_matched_tilde(); make_splines_matched();}

AngularAverageEta3PiC::AngularAverageEta3PiC(double mass_1, double mass_4,
                                        FunctionSet amplitude, int isospin,
                                        double cutoff, int iteration, double epsilon, double validity,
                                        SubtractionConstant sub)
:
AngularAverage{mass_1,mass_1,mass_1,mass_4,amplitude,isospin,cutoff,iteration,epsilon,validity,sub},
mandelstam_list(mass_1, mass_4, cutoff),
polar_egg(mass_1, mass_4, cutoff),
egg_deriv(mass_1, mass_4, cutoff)
{build_tilde_functions(); build_matching_list(); make_splines(); build_matched_tilde(); make_splines_matched();}

AngularAverageEtap3Pi::AngularAverageEtap3Pi(double mass_1, double mass_4,
                                        FunctionSet amplitude, int isospin,
                                        double cutoff, int iteration, double epsilon, double validity,
                                        SubtractionConstant sub)
:
AngularAverage{mass_1,mass_1,mass_1,mass_4,amplitude,isospin,cutoff,iteration,epsilon,validity,sub},
mandelstam_list(mass_1, mass_4, cutoff),
polar_egg(mass_1, mass_4, cutoff),
egg_deriv(mass_1, mass_4, cutoff)
{build_tilde_functions(); build_matching_list(); make_splines(); build_matched_tilde(); make_splines_matched();}

AngularAverageEtap3PiC::AngularAverageEtap3PiC(double mass_1, double mass_4,
                                        FunctionSet amplitude, int isospin,
                                        double cutoff, int iteration, double epsilon, double validity,
                                        SubtractionConstant sub)
:
AngularAverage{mass_1,mass_1,mass_1,mass_4,amplitude,isospin,cutoff,iteration,epsilon,validity,sub},
mandelstam_list(mass_1, mass_4, cutoff),
polar_egg(mass_1, mass_4, cutoff),
egg_deriv(mass_1, mass_4, cutoff)
{build_tilde_functions(); build_matching_list(); make_splines(); build_matched_tilde(); make_splines_matched();}

AngularAverageEtapEtaPiPi::AngularAverageEtapEtaPiPi(double mass_1, double mass_3, double mass_4,
                                        FunctionSet amplitude, int isospin,
                                        double cutoff, int iteration, double epsilon, double validity,
                                        SubtractionConstant sub)
:
AngularAverage{mass_1,mass_1,mass_3,mass_4,amplitude,isospin,cutoff,iteration,epsilon,validity,sub},
mandelstam_list(mass_1, mass_3, mass_4, cutoff),
contour(mass_1, mass_1, mass_3, mass_4),
matched_t_channel(mass_1,mass_3,mass_1,mass_4, epsilon, validity)
{build_tilde_functions(); build_matching_list(); make_splines(); build_matched_tilde(); make_splines_matched();}

AngularAverageEtapEtaPiPiC::AngularAverageEtapEtaPiPiC(double mass_1, double mass_3, double mass_4,
                                        FunctionSet amplitude, int isospin,
                                        double cutoff, int iteration, double epsilon, double validity,
                                        SubtractionConstant sub)
:
AngularAverage{mass_1,mass_1,mass_3,mass_4,amplitude,isospin,cutoff,iteration,epsilon,validity,sub},
mandelstam_list(mass_1, mass_3, mass_4, cutoff),
contour(mass_1, mass_1, mass_3, mass_4),
matched_t_channel(mass_1,mass_3,mass_1,mass_4, epsilon, validity)
{build_tilde_functions(); build_matching_list(); make_splines(); build_matched_tilde(); make_splines_matched();}

AngularAverageV3Pi::AngularAverageV3Pi(double mass_1, double mass_4,
                                        FunctionSet amplitude, int isospin,
                                        double cutoff, int iteration, double epsilon, double validity,
                                        SubtractionConstant sub)
:
AngularAverage{mass_1,mass_1,mass_1,mass_4,amplitude,isospin,cutoff,iteration,epsilon,validity,sub},
mandelstam_list(mass_1, mass_4, cutoff),
polar_egg(mass_1, mass_4, cutoff),
egg_deriv(mass_1, mass_4, cutoff)
{build_tilde_functions(); build_matching_list(); make_splines(); build_matched_tilde(); make_splines_matched();}

AngularAverageX3Pi::AngularAverageX3Pi(double mass_1, double mass_4,
                                        FunctionSet amplitude, int isospin,
                                        double cutoff, int iteration, double epsilon, double validity,
                                        SubtractionConstant sub)
:
AngularAverage{mass_1,mass_1,mass_1,mass_4,amplitude,isospin,cutoff,iteration,epsilon,validity,sub},
mandelstam_list(mass_1, mass_4, cutoff),
polar_egg(mass_1, mass_4, cutoff),
egg_deriv(mass_1, mass_4, cutoff)
{build_tilde_functions(); build_matching_list(); make_splines(); build_matched_tilde(); make_splines_matched();}

AngularAverageT3Pi::AngularAverageT3Pi(double mass_1, double mass_4,
                                        FunctionSet amplitude, int isospin,
                                        double cutoff, int iteration, double epsilon, double validity,
                                        SubtractionConstant sub)
:
AngularAverage{mass_1,mass_1,mass_1,mass_4,amplitude,isospin,cutoff,iteration,epsilon,validity,sub},
mandelstam_list(mass_1, mass_4, cutoff),
polar_egg(mass_1, mass_4, cutoff),
egg_deriv(mass_1, mass_4, cutoff)
{build_tilde_functions(); build_matching_list(); make_splines(); build_matched_tilde(); make_splines_matched();}

//------------------------------------------------------------------------
//-- Integrands ----------------------------------------------------------
//------------------------------------------------------------------------

Complex AngularAverageEta3Pi::integrand(double s, double x, Setting evaluation, TypeOfAverage aver)const{
    double phi;
    if (isospin==0) {
        switch (evaluation) {
            case Setting::above: case Setting::below:
                return 1/9.* (  6.*amplitude(x, 0, sub, evaluation)
                             +18.*(s - s_0) *amplitude(x, 1, sub, evaluation)
                             +12.*(x-(3*s_0-s)/2.)*amplitude(x, 1, sub, evaluation)
                              +20.*amplitude(x, 2, sub, evaluation) );
                
            // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
            case Setting::egg:
                phi=x;
                return 1/9.* (  6.*amplitude(phi, 0, sub, evaluation)
                             +18.*(s - s_0) *amplitude(phi, 1, sub, evaluation)
                             +12.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 1, sub, evaluation)
                              +20.*amplitude(phi, 2, sub, evaluation) )
                * egg_deriv.spline_egg_derivative(phi);
                
            default:
                throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
        }
    }
    if (isospin==1) {
        switch (evaluation) {
            case Setting::above: case Setting::below:
                return 6.*(x-(3*s_0-s)/2.)*amplitude(x, 0, sub, evaluation)
                        +9.*(s - s_0)*(x-(3*s_0-s)/2.)*amplitude(x, 1, sub, evaluation)
                        +6.*std::pow(x-(3*s_0-s)/2.,2.)*amplitude(x, 1, sub, evaluation)
                        -10.*(x-(3*s_0-s)/2.)*amplitude(x, 2, sub, evaluation);
                
            // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
            case Setting::egg:
                phi=x;
                return (  6.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 0, sub, evaluation)
                        +9.*(s - s_0)*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 1, sub, evaluation)
                        +6.*std::pow(polar_egg.spline_egg(phi)-(3*s_0-s)/2.,2.)*amplitude(phi, 1, sub, evaluation)
                        -10.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 2, sub, evaluation))
                * egg_deriv.spline_egg_derivative(phi);
                
            default:
                throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
        }
    }
    if (isospin==2) {
        switch (evaluation) {
            case Setting::above: case Setting::below:
                return 1/6.* ( 6.*amplitude(x, 0, sub, evaluation)
                             -9.*(s - s_0) *amplitude(x, 1, sub, evaluation)
                             -6.*(x-(3*s_0-s)/2.)*amplitude(x, 1, sub, evaluation)
                             +2.*amplitude(x, 2, sub, evaluation) );
                
            // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
            case Setting::egg:
                phi=x;
                return 1/6.* ( 6.*amplitude(phi, 0, sub, evaluation)
                             -9.*(s - s_0) *amplitude(phi, 1, sub, evaluation)
                             -6.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 1, sub, evaluation)
                             +2.*amplitude(phi, 2, sub, evaluation) )
                * egg_deriv.spline_egg_derivative(phi);
                
            default:
                throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
        }
    }
    throw std::domain_error{"For total isospin 1 the two-pion state must have isospin 0, 1 or 2"};
}


Complex AngularAverageEta3PiC::integrand(double s, double x, Setting evaluation, TypeOfAverage aver)const{
    double phi;
    switch (sub) {
            // Beyond Standard Model (C-violating, total isospin=0)
            case SubtractionConstant::g1:
                if (isospin==1) {
                    switch (evaluation) {
                        case Setting::above: case Setting::below:
                            return -6.* ( 3.*(s-s_0)*(x-(3*s_0-s)/2.)
                                        + 2.*std::pow(x-(3*s_0-s)/2.,2.))*amplitude(x, 1, sub, evaluation);
                            
                        // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
                        case Setting::egg:
                            phi=x;
                            return -6.* ( 3.*(s-s_0)*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)
                                        + 2.*std::pow(polar_egg.spline_egg(phi)-(3*s_0-s)/2.,2.))*amplitude(phi, 1, sub, evaluation)
                                        * egg_deriv.spline_egg_derivative(phi);
                            
                        default:
                            throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
                    }
                }
                else{throw std::domain_error{"For total isospin 0 the two-pion state must have isospin 1"};};
                
            // Beyond Standard Model (C-violating, total isospin=2)
            case SubtractionConstant::h1:
                if (isospin==1) {
                    switch (evaluation) {
                        case Setting::above: case Setting::below:
                            return 3.* ( 3.*(s-s_0)*(x-(3*s_0-s)/2.)*amplitude(x, 1, sub, evaluation)
                                        +2.*std::pow(x-(3*s_0-s)/2.,2.)*amplitude(x, 1, sub, evaluation)
                                        +2.*(x-(3*s_0-s)/2.)*amplitude(x, 2, sub, evaluation));
                            
                        // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
                        case Setting::egg:
                            phi=x;
                            return 3.* ( 3.*(s-s_0)*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 1, sub, evaluation)
                                        +2.*std::pow(polar_egg.spline_egg(phi)-(3*s_0-s)/2.,2.)*amplitude(phi, 1, sub, evaluation)
                                        +2.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 2, sub, evaluation))
                            * egg_deriv.spline_egg_derivative(phi);
                            
                        default:
                            throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
                    }
                }
                if (isospin==2) {
                    switch (evaluation) {
                        case Setting::above: case Setting::below:
                            return 1/2.* ( 9.*(s-s_0)*amplitude(x, 1, sub, evaluation)
                                          +6.*(x-(3*s_0-s)/2.)*amplitude(x, 1, sub, evaluation)
                                          -2.*amplitude(x, 2, sub, evaluation));
                            
                        // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
                        case Setting::egg:
                            phi=x;
                            return 1/2.* ( 9.*(s-s_0)*amplitude(phi, 1, sub, evaluation)
                                          +6.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 1, sub, evaluation)
                                          -2.*amplitude(phi, 2, sub, evaluation))
                            * egg_deriv.spline_egg_derivative(phi);
                            
                        default:
                            throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
                    }
                }
                throw std::domain_error{"For total isospin 2 the two-pion state must have isospin 1 or 2"};
            
        default:
            throw std::domain_error{"Subtraction constant does not belong to the chosen integrand."};
    }
}

Complex AngularAverageEtap3Pi::integrand(double s, double x, Setting evaluation, TypeOfAverage aver)const{
    double phi;
    if (isospin==0) {
        switch (evaluation) {
            case Setting::above: case Setting::below:
                return 1/9.* (  6.*amplitude(x, 0, sub, evaluation)
                             +18.*(s - s_0) *amplitude(x, 1, sub, evaluation)
                             +12.*(x-(3*s_0-s)/2.)*amplitude(x, 1, sub, evaluation)
                              +20.*amplitude(x, 2, sub, evaluation) );
                
            // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
            case Setting::egg:
                phi=x;
                return 1/9.* (  6.*amplitude(phi, 0, sub, evaluation)
                             +18.*(s - s_0) *amplitude(phi, 1, sub, evaluation)
                             +12.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 1, sub, evaluation)
                              +20.*amplitude(phi, 2, sub, evaluation) )
                * egg_deriv.spline_egg_derivative(phi);
                
            default:
                throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
        }
    }
    if (isospin==1) {
        switch (evaluation) {
            case Setting::above: case Setting::below:
                return 6.*(x-(3*s_0-s)/2.)*amplitude(x, 0, sub, evaluation)
                        +9.*(s - s_0)*(x-(3*s_0-s)/2.)*amplitude(x, 1, sub, evaluation)
                        +6.*std::pow(x-(3*s_0-s)/2.,2.)*amplitude(x, 1, sub, evaluation)
                        -10.*(x-(3*s_0-s)/2.)*amplitude(x, 2, sub, evaluation);
                
            // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
            case Setting::egg:
                phi=x;
                return (  6.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 0, sub, evaluation)
                        +9.*(s - s_0)*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 1, sub, evaluation)
                        +6.*std::pow(polar_egg.spline_egg(phi)-(3*s_0-s)/2.,2.)*amplitude(phi, 1, sub, evaluation)
                        -10.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 2, sub, evaluation))
                * egg_deriv.spline_egg_derivative(phi);
                
            default:
                throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
        }
    }
    if (isospin==2) {
        switch (evaluation) {
            case Setting::above: case Setting::below:
                return 1/6.* ( 6.*amplitude(x, 0, sub, evaluation)
                             -9.*(s - s_0) *amplitude(x, 1, sub, evaluation)
                             -6.*(x-(3*s_0-s)/2.)*amplitude(x, 1, sub, evaluation)
                             +2.*amplitude(x, 2, sub, evaluation) );
                
            // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
            case Setting::egg:
                phi=x;
                return 1/6.* ( 6.*amplitude(phi, 0, sub, evaluation)
                             -9.*(s - s_0) *amplitude(phi, 1, sub, evaluation)
                             -6.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 1, sub, evaluation)
                             +2.*amplitude(phi, 2, sub, evaluation) )
                * egg_deriv.spline_egg_derivative(phi);
                
            default:
                throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
        }
    }
    throw std::domain_error{"For total isospin 1 the two-pion state must have isospin 0, 1 or 2"};
}


Complex AngularAverageEtap3PiC::integrand(double s, double x, Setting evaluation, TypeOfAverage aver)const{
    double phi;
    switch (sub) {
            // Beyond Standard Model (C-violating, total isospin=0)
            case SubtractionConstant::g1:
                if (isospin==1) {
                    switch (evaluation) {
                        case Setting::above: case Setting::below:
                            return  -6.* ( 3.*(s-s_0)*(x-(3*s_0-s)/2.)
                                        + 2.*std::pow(x-(3*s_0-s)/2.,2.))*amplitude(x, 1, sub, evaluation);
                            
                        // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
                        case Setting::egg:
                            phi=x;
                            return  -6.* ( 3.*(s-s_0)*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)
                                        + 2.*std::pow(polar_egg.spline_egg(phi)-(3*s_0-s)/2.,2.))*amplitude(phi, 1, sub, evaluation)
                                        * egg_deriv.spline_egg_derivative(phi);
                            
                        default:
                            throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
                    }
                }
                else{throw std::domain_error{"For total isospin 0 the two-pion state must have isospin 1"};};
                
            // Beyond Standard Model (C-violating, total isospin=2)
            case SubtractionConstant::h1:
                if (isospin==1) {
                    switch (evaluation) {
                        case Setting::above: case Setting::below:
                            return 3.* ( 3.*(s-s_0)*(x-(3*s_0-s)/2.)*amplitude(x, 1, sub, evaluation)
                                        +2.*std::pow(x-(3*s_0-s)/2.,2.)*amplitude(x, 1, sub, evaluation)
                                        +2.*(x-(3*s_0-s)/2.)*amplitude(x, 2, sub, evaluation));
                            
                        // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
                        case Setting::egg:
                            phi=x;
                            return 3.* ( 3.*(s-s_0)*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 1, sub, evaluation)
                                        +2.*std::pow(polar_egg.spline_egg(phi)-(3*s_0-s)/2.,2.)*amplitude(phi, 1, sub, evaluation)
                                        +2.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 2, sub, evaluation))
                            * egg_deriv.spline_egg_derivative(phi);
                            
                        default:
                            throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
                    }
                }
                if (isospin==2) {
                    switch (evaluation) {
                        case Setting::above: case Setting::below:
                            return 1/2.* ( 9.*(s-s_0)*amplitude(x, 1, sub, evaluation)
                                          +6.*(x-(3*s_0-s)/2.)*amplitude(x, 1, sub, evaluation)
                                          -2.*amplitude(x, 2, sub, evaluation));
                            
                        // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
                        case Setting::egg:
                            phi=x;
                            return 1/2.* ( 9.*(s-s_0)*amplitude(phi, 1, sub, evaluation)
                                          +6.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.)*amplitude(phi, 1, sub, evaluation)
                                          -2.*amplitude(phi, 2, sub, evaluation))
                            * egg_deriv.spline_egg_derivative(phi);
                            
                        default:
                            throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
                    }
                }
                throw std::domain_error{"For total isospin 2 the two-pion state must have isospin 1 or 2"};
            
        default:
            throw std::domain_error{"Subtraction constant does not belong to the chosen integrand."};
    }
}

Complex AngularAverageEtapEtaPiPi::integrand(double s, double x, Setting evaluation, TypeOfAverage aver)const{
    switch (aver) {
        case TypeOfAverage::zero:
            // since the amplitudes for this decay can not be uniquely identified by their isospin,
            // use instead the combination: total isospin, partial-wave, isospin of the two-particle state
            if (isospin==000) {
                switch (evaluation) {
                    case Setting::above: case Setting::below:
                        return 2.* amplitude(x, 010, sub, evaluation);
                        
                    case Setting::egg:
                        throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
                        
                    // along the egg-like part of the complex contour with curve-parameter y
                    default:
                        double y=x;
                        return 2.* amplitude(y, 010, sub, evaluation) * contour.derivative(y, evaluation);
                }
            }
            throw std::domain_error{"Isospin of two-particle state is not allowed for the chosen configuration."};
            
            
        case TypeOfAverage::minus:
            if (isospin==010) {
                switch (evaluation) {
                    case Setting::above: case Setting::below:
                        return amplitude(x, 000, sub, evaluation);
                        
                    case Setting::egg:
                        throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
                        
                    // along the egg-like part of the complex contour with curve-parameter y
                    default:
                        double y=x;
                        return amplitude(y, 000, sub, evaluation)
                        * contour.derivative(y, evaluation);
                }
            }
            throw std::domain_error{"Isospin of two-particle state is not allowed for the chosen configuration."};
            
            
        case TypeOfAverage::plus:
            if (isospin==010) {
                switch (evaluation) {
                    case Setting::above: case Setting::below:
                        return amplitude(x, 010, sub, evaluation);
                        
                    case Setting::egg:
                        throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
                        
                    // along the egg-like part of the complex contour with curve-parameter y
                    default:
                        double y=x;
                        return amplitude(y, 010, sub, evaluation)
                        * contour.derivative(y, evaluation);
                }
            }
            throw std::domain_error{"Isospin of two-particle state is not allowed for the chosen configuration."};

        default:
            throw std::domain_error{"Average not defined properly"};
    }
}


Complex AngularAverageEtapEtaPiPiC::integrand(double s, double x, Setting evaluation, TypeOfAverage aver)const{
    switch (aver) {
        case TypeOfAverage::zero:
            // since the amplitudes for this decay can not be uniquely identified by their isospin,
            // use instead the combination: total isospin, partial-wave, isospin of the two-particle state
            if (isospin==111) {
                switch (evaluation) {
                    case Setting::above: case Setting::below:
                        return 6. * (2.*x -3.*s_0 + s) * amplitude(x, 110, sub, evaluation);
                        
                    case Setting::egg:
                        throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
                        
                    // along the egg-like part of the complex contour with curve-parameter y
                    default:
                        double y=x;
                        return
                        6. * (2.*contour(y,evaluation) -3.*s_0 + s)
                        * amplitude(y, 110, sub, evaluation)
                        * contour.derivative(y, evaluation);
                }
            }
            throw std::domain_error{"Isospin of two-particle state is not allowed for the chosen configuration."};
            
            
        case TypeOfAverage::minus:
            if (isospin==110) {
                switch (evaluation) {
                    case Setting::above: case Setting::below:
                        return
                        -3./2.*(t_0 -s + contour.delta/(3.*s))
                        *amplitude(x, 111, sub, evaluation)
                        +1./2.*(2.*x -3.*t_0 + s + contour.delta/s)
                        *amplitude(x, 111, sub, evaluation);
                        
                    case Setting::egg:
                        throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
                        
                    // along the egg-like part of the complex contour with curve-parameter y
                    default:
                        double y=x;
                        return
                        (-3./2.*(t_0 -s + contour.delta/(3.*s))
                        *amplitude(y, 111, sub, evaluation)
                        +1./2.*(2.*contour(y,evaluation) -3.*t_0 + s +contour.delta/s)
                        *amplitude(y, 111, sub, evaluation))
                        * contour.derivative(y, evaluation);
                }
            }
            
            throw std::domain_error{"Isospin of two-particle state is not allowed for the chosen configuration."};
            
            
        case TypeOfAverage::plus:
            if (isospin==110) {
                switch (evaluation) {
                    case Setting::above: case Setting::below:
                        return -1.*amplitude(x, 110, sub, evaluation);
                        
                    case Setting::egg:
                        throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
                        
                    // along the egg-like part of the complex contour with curve-parameter y
                    default:
                        double y=x;
                        return
                        -1.*amplitude(y, 110, sub, evaluation)
                        * contour.derivative(y, evaluation);
                }
            }
            throw std::domain_error{"Isospin of two-particle state is not allowed for the chosen configuration."};

        default:
            throw std::domain_error{"Average not defined properly"};
    }
}

Complex AngularAverageV3Pi::integrand(double s, double x, Setting evaluation, TypeOfAverage aver)const{
    double phi;
    if (isospin==1) {
        switch (evaluation) {
            case Setting::above: case Setting::below:
                return (3.*(1-4*std::pow(mass_1,2)/s)*(kallen(s,mass_4,mass_1)) //multiply with kappa(s)^2 to render both terms equal in kappa
                        -12.*std::pow(x-(3*s_0-s)/2.,2.))*amplitude(x, 1, sub, evaluation);
                
            // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
            case Setting::egg:
                phi=x;
                return  ((3.*(1-4*std::pow(mass_1,2)/s)*(kallen(s,mass_4,mass_1)) //multiply with kappa(s)^2 to render both terms equal in kappa
                        -12.*std::pow(polar_egg.spline_egg(phi)-(3*s_0-s)/2.,2.))*amplitude(phi, 1, sub, evaluation) )
                        * egg_deriv.spline_egg_derivative(phi);
                
            default:
                throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
        }
    }
    throw std::domain_error{"For total isospin 1 the two-pion state must have isospin 1."};
}

Complex AngularAverageX3Pi::integrand(double s, double x, Setting evaluation, TypeOfAverage aver)const{
    double phi;
    if (isospin==1) {
        switch (evaluation) {
            case Setting::above: case Setting::below:
                return -1./2*(3.*(1-4*std::pow(mass_1,2)/s)*(kallen(s,mass_4,mass_1)) //multiply with kappa(s)^2 to render both terms equal in kappa
                        -12.*std::pow(x-(3*s_0-s)/2.,2.))*amplitude(x, 1, sub, evaluation);
                
            // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
            case Setting::egg:
                phi=x;
                return  -1./2*((3.*(1-4*std::pow(mass_1,2)/s)*(kallen(s,mass_4,mass_1)) //multiply with kappa(s)^2 to render both terms equal in kappa
                        -12.*std::pow(polar_egg.spline_egg(phi)-(3*s_0-s)/2.,2.))*amplitude(phi, 1, sub, evaluation) )
                        * egg_deriv.spline_egg_derivative(phi);
                
            default:
                throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
        }
    }
    throw std::domain_error{"For total isospin 1 the two-pion state must have isospin 1."};
}

Complex AngularAverageT3Pi::integrand(double s, double x, Setting evaluation, TypeOfAverage aver)const{
    double phi;
    Complex z;
    Complex kappa2;
    Complex prefactor;
    if (isospin==1) {
        kappa2 = (1-4*std::pow(mass_1,2)/s)*kallen(s,mass_4,mass_1);
        switch (evaluation) {
            case Setting::above: case Setting::below:
                z = 2.*(x-(3*s_0-s)/2.);
                prefactor = (1-4*std::pow(mass_1,2)/s)*(s+std::pow(mass_4,2)-std::pow(mass_1,2));
                return 3./4*(kappa2-std::pow(z,2.))*amplitude(x, 1, sub, evaluation)*kappa2 //multiply with kappa(s)^2 to render both terms equal in kappa
                - 3./4* prefactor * (kappa2-std::pow(z,2.)) * z * amplitude(x, 1, sub, evaluation);
                
            // along the egg-like part of the complex contour rewrite in terms of real-valued angle 'phi'
            case Setting::egg:
                phi=x;
                z = 2.*(polar_egg.spline_egg(phi)-(3*s_0-s)/2.);
                prefactor = (1-4*std::pow(mass_1,2)/s)*(s+std::pow(mass_4,2)-std::pow(mass_1,2));
                return  (3./4*(kappa2-std::pow(z,2.))*amplitude(phi, 1, sub, evaluation)*kappa2 //multiply with kappa(s)^2 to render both terms equal in kappa
                        - 3./4* prefactor * (kappa2-std::pow(z,2.)) * z * amplitude(phi, 1, sub, evaluation))
                        * egg_deriv.spline_egg_derivative(phi);
                
            default:
                throw std::domain_error{"Demanded evaluation is not needed for the basis amplitude you are calling."};
        }
    }
    throw std::domain_error{"For total isospin 1 the two-pion state must have isospin 1."};
}


//------------------------------------------------------------------------
//-- Integral ------------------------------------------------------------
//------------------------------------------------------------------------

// need to introduce lambda-expressions for the basis amplitudes, otherwise they can not be used
// as arguments of the member functions 'integrand_...'.

// Here one can switch between an adaptive integration routine and GaussLegrendre quadrature with
// an arbitrary number of integrations nodes
Complex AngularAverage::integral_3em(double s, const splined_path::SplinedPath &polar_egg)const{
    
    double mandelstam_s=s;
    
    if (s_I<=s && s<s_II) {
        // write in terms of only one variable
        const auto integrand_above{[mandelstam_s,this](double s_prime)
            { return this->integrand(mandelstam_s, s_prime, Setting::above); }};
            
        return complex_integration(integrand_above, polar_egg.s_lower(s), polar_egg.s_upper(s), false, 100);
    }

    if (s_II<=s && s<=s_III){
        // write in terms of only one variable
        const auto integrand_above{[mandelstam_s,this](double s_prime)
            { return this->integrand(mandelstam_s, s_prime, Setting::above); }};
            
        const auto integrand_below{[mandelstam_s,this](double s_prime)
            { return this->integrand(mandelstam_s, s_prime, Setting::below); }};

        return
        complex_integration(integrand_below, polar_egg.s_lower(s), s_I, false, 300)
        +complex_integration(integrand_above, s_I, polar_egg.s_upper(s), false, 300);
    }

    if (s_III<s && s<s_IV) {
        // write in terms of only one variable
        const auto integrand_egg{[mandelstam_s,this](double s_prime)
            { return this->integrand(mandelstam_s, s_prime, Setting::egg); }};

        Complex integral=complex_integration(integrand_egg, polar_egg.spline_phi_lower(mandelstam_s), polar_egg.spline_phi_upper(mandelstam_s), false, 300);

        // As the Omnès function is perfectly Schwartz the real part of the
        // angular average vanishes in this region for the first iteration
        if (iteration==0) {
            return 1.0i*imag(integral);
        }
        return integral;
    }
    if (s_IV<=s) {
        const auto integrand_above{[mandelstam_s,this](double s_prime)
            { return this->integrand(mandelstam_s, s_prime, Setting::above); }};
        
        Complex integral=complex_integration(integrand_above, polar_egg.s_lower(s), polar_egg.s_upper(s), false, 500);

        // As the Omnès function is perfectly Schwartz the real part of the
        // angular average vanishes in this region for the first iteration
        if (iteration==0) {
            return real(integral)+0.0i;
        }
        return integral;
    }
    throw std::domain_error{"Value for Mandelstam 's' in tilde_function is not allowed 1!"};
}


Complex AngularAverage::integral_2em(double s, TypeOfAverage aver, const path_eta_pi_pi::Path &contour)const{
    TypeOfAverage type=aver; // introduce type to capture s in lambda expression
    double x=s;// introduce x to capture s in lambda expression

    auto [x_I, x_II, x_III, x_IV] = [aver,this]() {
        switch (aver) {
            case TypeOfAverage::zero:
                return std::make_tuple(s_I, s_II, s_III, s_IV);
            case TypeOfAverage::minus:
                return std::make_tuple(t_I,t_II,t_III,t_IV);
            case TypeOfAverage::plus:
                return std::make_tuple(t_I,u_II,t_III,t_IV);
            default:
                throw std::domain_error{"Type of angular average has to be either zero, minus or plus."};
        }
    }();
    
    if (x_I<=x && x<x_II) {
        // write in terms of only one variable
        const auto integrand_above{[x,type,this](double s_prime)
            { return this->integrand(x, s_prime, Setting::above, type); }};

        return complex_integration(integrand_above,
                                   contour.integration_limit(x, type, PinocchioEndpoints::lower),
                                   contour.integration_limit(x, type, PinocchioEndpoints::upper), true, 100);
    }


    if (x_II<=x && x<=x_III){
        // write in terms of only one variable
        const auto integrand_above{[x,type,this](double s_prime)
            { return this->integrand(x, s_prime, Setting::above, type); }};

        const auto integrand_below{[x,type,this](double s_prime)
            { return this->integrand(x, s_prime, Setting::below, type); }};
        
        // Watch out for the different sign in +-i*epsilon for the different averages! (only need that distinction in this integration region)
        switch (type) {
            case TypeOfAverage::zero: case TypeOfAverage::plus:
                return
                complex_integration(integrand_below, contour.integration_limit(x, type, PinocchioEndpoints::lower), x_I, false, 100)
                +complex_integration(integrand_above, x_I, contour.integration_limit(x, type, PinocchioEndpoints::upper), false, 100);
                
                // HERE IS ACTUALLY A MISTAKE:
                // ONE SHOULD EVALUATE THE INTEGRAL ABOVE AND BELOW THE REAL AXIS.
                // BUT THIS RESULTS IN AN UNSTEADY IMAGINARY PART OF THE AVERAGE. EVALUATING THE INTEGRAL ONLY ABOVE THE REAL AXIS COMPENSATES THE SIGN ERROR, WHICH APPEARS SOMEWHERE IN THE REMAINING CODE.
                // NOTE THAT THIS ERROR ONLY AFFECTS THE eta'-> eta pi pi DECAY AND ONLY THIS INTEGRATION REGION.
            case TypeOfAverage::minus:
                return
                complex_integration(integrand_above, contour.integration_limit(x, type, PinocchioEndpoints::lower), x_I, true, 100)
                +complex_integration(integrand_above, x_I, contour.integration_limit(x, type, PinocchioEndpoints::upper), true, 100);
                
            default:
                throw std::domain_error{"Type of angular average has to be either zero, minus or plus."};
        }
        
    }

    
    if (x_III<x && x<x_IV) {

        if (type==TypeOfAverage::zero) {
                const auto integrand_zero_upper_a{[x,type,this](double s_prime)
                    { return this->integrand(x, s_prime, Setting::zero_upper_a, type); }};
                const auto integrand_zero_upper_b{[x,type,this](double s_prime)
                    { return this->integrand(x, s_prime, Setting::zero_upper_b, type); }};
                const auto integrand_zero_lower_a{[x,type,this](double s_prime)
                    { return this->integrand(x, s_prime, Setting::zero_lower_a, type); }};
                const auto integrand_zero_lower_b{[x,type,this](double s_prime)
                    { return this->integrand(x, s_prime, Setting::zero_lower_b, type); }};

                return helper(x, type, integrand_zero_lower_a, integrand_zero_upper_a, integrand_zero_lower_b, integrand_zero_upper_b, contour);
        }
        if (type==TypeOfAverage::plus) {
                const auto integrand_plus_upper_a{[x,type,this](double s_prime)
                    { return this->integrand(x, s_prime, Setting::plus_upper_a, type); }};
                const auto integrand_plus_upper_b{[x,type,this](double s_prime)
                    { return this->integrand(x, s_prime, Setting::plus_upper_b, type); }};
                const auto integrand_plus_lower_a{[x,type,this](double s_prime)
                    { return this->integrand(x, s_prime, Setting::plus_lower_a, type); }};
                const auto integrand_plus_lower_b{[x,type,this](double s_prime)
                    { return this->integrand(x, s_prime, Setting::plus_lower_b, type); }};

                return helper(x, type, integrand_plus_lower_a, integrand_plus_upper_a, integrand_plus_lower_b, integrand_plus_upper_b, contour);
        }
        if (type==TypeOfAverage::minus) {
                const auto integrand_minus_upper_a{[x,type,this](double s_prime)
                    { return this->integrand(x, s_prime, Setting::minus_upper_a, type); }};
                const auto integrand_minus_upper_b{[x,type,this](double s_prime)
                    { return this->integrand(x, s_prime, Setting::minus_upper_b, type); }};
                const auto integrand_minus_lower_a{[x,type,this](double s_prime)
                    { return this->integrand(x, s_prime, Setting::minus_lower_a, type); }};
                const auto integrand_minus_lower_b{[x,type,this](double s_prime)
                    { return this->integrand(x, s_prime, Setting::minus_lower_b, type); }};

                return helper(x, type, integrand_minus_lower_a, integrand_minus_upper_a, integrand_minus_lower_b, integrand_minus_upper_b, contour);
        }
        throw std::domain_error{"Type of angular average has to be either zero, minus or plus."};
    }


    if (x_IV<=x) {
        const auto integrand_above{[x,type,this](double s_prime)
            { return this->integrand(x, s_prime, Setting::above, type); }};

        Complex integral=complex_integration(integrand_above,
                                             contour.integration_limit(x, type, PinocchioEndpoints::lower),
                                             contour.integration_limit(x, type, PinocchioEndpoints::upper), false, 500);

        // As the Omnès function is perfectly Schwartz the real part of the
        // angular average vanishes in this region for the first iteration
//        if (iteration==0) {
//            return real(integral)+0.0i;
//        }
        return integral;
    }
    throw std::domain_error{"Value for Mandelstam 's' in tilde_function is not allowed!"};
}


Complex AngularAverage::helper(double x, TypeOfAverage aver, Curve integrand_lower_a, Curve integrand_upper_a, Curve integrand_lower_b, Curve integrand_upper_b, const path_eta_pi_pi::Path &contour)const{

    auto [x_I, x_II, x_III, x_IV] = [aver,this]() {
        switch (aver) {
            case TypeOfAverage::zero:
                return std::make_tuple(s_I, s_II, s_III, s_IV);
            case TypeOfAverage::minus:
                return std::make_tuple(t_I,t_II,t_III,t_IV);
            case TypeOfAverage::plus:
                return std::make_tuple(t_I,u_II,t_III,t_IV);
            default:
                throw std::domain_error{"Type of angular average has to be either zero, minus or plus."};
        }
    }();


    double p=(x_III+x_IV)/2.;
    if (p<=x && x< x_IV ) {
        
        Complex integral_lower=
        complex_integration(integrand_lower_b,
                            contour.yb(x),
                            contour.yb(x_IV), false, 100);
        Complex integral_upper=
        complex_integration(integrand_upper_b,
                            contour.yb(x_IV),
                            contour.yb(x), false, 100);

        // As the Omnès function is perfectly Schwartz the real part of the
        // angular average vanishes in this region for the first iteration
//            if (iteration==0) {
//                return 1.0i*imag(integral_lower+integral_upper);
//            }
        return integral_lower+integral_upper;
        
    }// if: p<=x && x< x_IV
    
    if (x_III<x && x< p) {
        
        Complex integral_lower_a=
        complex_integration(integrand_lower_a,
                            contour.ya(x),
                            contour.ya(p), false, 100);
        Complex integral_upper_a=
        complex_integration(integrand_upper_a,
                            contour.ya(p),
                            contour.ya(x), false, 100);
        Complex integral_lower_b=
        complex_integration(integrand_lower_b,
                            contour.yb(p),
                            contour.yb(x_IV), false, 100);
        Complex integral_upper_b=
        complex_integration(integrand_upper_b,
                            contour.yb(x_IV),
                            contour.yb(p), false, 100);

        // As the Omnès function is perfectly Schwartz the real part of the
        // angular average vanishes in this region for the first iteration
//            if (iteration==0) {
//                return 1.0i*imag(integral);
//            }
        return
         integral_lower_a + integral_upper_a
        +integral_lower_b + integral_upper_b;
        
    }// if: s_III<x && x< p
    throw std::domain_error{"Helper for angular averages is not defined for the inserted value of s"};
}
//------------------------------------------------------------------------
//------------------------------------------------------------------------
Complex AngularAverageEta3Pi::eval_int(double s, TypeOfAverage aver)const{
    return integral_3em(s, polar_egg);
}

Complex AngularAverageEta3PiC::eval_int(double s, TypeOfAverage aver)const{
    return integral_3em(s, polar_egg);
}

Complex AngularAverageEtap3Pi::eval_int(double s, TypeOfAverage aver)const{
    return integral_3em(s, polar_egg);
}

Complex AngularAverageEtap3PiC::eval_int(double s, TypeOfAverage aver)const{
    return integral_3em(s, polar_egg);
}

Complex AngularAverageEtapEtaPiPi::eval_int(double s, TypeOfAverage aver)const{
    return integral_2em(s, aver, contour);
}

Complex AngularAverageEtapEtaPiPiC::eval_int(double s, TypeOfAverage aver)const{
    return integral_2em(s, aver, contour);
}

Complex AngularAverageV3Pi::eval_int(double s, TypeOfAverage aver)const{
    return integral_3em(s, polar_egg);
}

Complex AngularAverageX3Pi::eval_int(double s, TypeOfAverage aver)const{
    return integral_3em(s, polar_egg);
}

Complex AngularAverageT3Pi::eval_int(double s, TypeOfAverage aver)const{
    return integral_3em(s, polar_egg);
}


//------------------------------------------------------------------------
//-- Build ---------------------------------------------------------------
//------------------------------------------------------------------------

// One could use 'interval_s_0_th' instead of 'interval_s_th' to use more points around the cusp
// in the isospin 0 phase shift (be careful to replace ALL 'interval_s_th' in this file)

// if the isospin and decay are not the desired ones given in the constructor fill hte corresponding lists with zeros
void AngularAverageEta3Pi::build_tilde_functions(){
    // evaluate integrals along the array in s
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); j++) {
        switch (isospin) {
            case 0: case 1: case 2:
                // for eta-> 3 pi
                angular_average_list.push_back(eval_int(mandelstam_list.interval_s_th[j]));
                break;
                        
            default:
                throw std::domain_error{"Wrong isospin for the chosen Process!"};
        }
    }
}

void AngularAverageEta3PiC::build_tilde_functions(){
    // evaluate integrals along the array in s
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); j++) {
        switch (isospin) {
            case 0: case 1: case 2:
                // for eta-> 3 pi
                angular_average_list.push_back(eval_int(mandelstam_list.interval_s_th[j]));
                break;
                        
            default:
                throw std::domain_error{"Wrong isospin for the chosen Process!"};
        }
    }
}

void AngularAverageEtap3Pi::build_tilde_functions(){
    // evaluate integrals along the array in s
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); j++) {
        switch (isospin) {
            case 0: case 1: case 2:
                // for eta-> 3 pi
                angular_average_list.push_back(eval_int(mandelstam_list.interval_s_th[j]));
                break;
                        
            default:
                throw std::domain_error{"Wrong isospin for the chosen Process!"};
        }
    }
}

void AngularAverageEtap3PiC::build_tilde_functions(){
    // evaluate integrals along the array in s
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); j++) {
        switch (isospin) {
            case 0: case 1: case 2:
                // for eta-> 3 pi
                angular_average_list.push_back(eval_int(mandelstam_list.interval_s_th[j]));
                break;
                        
            default:
                throw std::domain_error{"Wrong isospin for the chosen Process!"};
        }
    }
}

void AngularAverageEtapEtaPiPi::build_tilde_functions(){
    // evaluate integrals along the array in s
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); j++) {
        switch (isospin) {
            case 000: case 111:
                // for eta'-> eta pi pi with pi-pi intermediate state
                angular_average_list.push_back(eval_int(mandelstam_list.interval_s_th[j],TypeOfAverage::zero));
                break;
                        
                case 010: case 110: // these isospins are evaluated along another grid, therefore leave them out at this point
                    angular_average_list.push_back(100.);
                    break;
                        
            default:
                throw std::domain_error{"Wrong isospin for the chosen Process!"};
        }
    }
    
    for (std::size_t j=0; j<mandelstam_list.interval_t_th.size(); j++) {
        
        switch (isospin) {
            case 010: case 110:
                // for eta'-> eta pi pi with pi-eta intermediate state
                angular_average_intermediate_list.push_back(
                                                      eval_int(mandelstam_list.interval_t_th[j],TypeOfAverage::minus)
                                                      +eval_int(mandelstam_list.interval_t_th[j],TypeOfAverage::plus)
                                                      );
                break;
                
            default:
                angular_average_intermediate_list.push_back(100.);
                break;
        }
        
    }
    return;
}

void AngularAverageEtapEtaPiPiC::build_tilde_functions(){
    // evaluate integrals along the array in s
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); j++) {
        switch (isospin) {
            case 000: case 111:
                // for eta'-> eta pi pi with pi-pi intermediate state
                angular_average_list.push_back(eval_int(mandelstam_list.interval_s_th[j],TypeOfAverage::zero));
                break;
                        
                case 010: case 110: // these isospins are evaluated along another grid, therefore leave them out at this point
                    angular_average_list.push_back(100.);
                    break;
                        
            default:
                throw std::domain_error{"Wrong isospin for the chosen Process!"};
        }
    }
    
    for (std::size_t j=0; j<mandelstam_list.interval_t_th.size(); j++) {
        
        switch (isospin) {
            case 010: case 110:
                // for eta'-> eta pi pi with pi-eta intermediate state
                angular_average_intermediate_list.push_back(
                                                      eval_int(mandelstam_list.interval_t_th[j],TypeOfAverage::minus)
                                                      +eval_int(mandelstam_list.interval_t_th[j],TypeOfAverage::plus)
                                                      );
                break;
                
            default:
                angular_average_intermediate_list.push_back(100.);
                break;
        }
        
    }
    return;
}

void AngularAverageV3Pi::build_tilde_functions(){
    // evaluate integrals along the array in s
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); j++) {
        switch (isospin) {
            case 1:
                angular_average_list.push_back(eval_int(mandelstam_list.interval_s_th[j]));
                break;
                        
            default:
                throw std::domain_error{"Wrong isospin for the chosen Process!"};
        }
    }
}

void AngularAverageX3Pi::build_tilde_functions(){
    // evaluate integrals along the array in s
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); j++) {
        switch (isospin) {
            case 1:
                // for X-> 3 pi
                angular_average_list.push_back(eval_int(mandelstam_list.interval_s_th[j]));
                break;
                        
            default:
                throw std::domain_error{"Wrong isospin for the chosen Process!"};
        }
    }
}

void AngularAverageT3Pi::build_tilde_functions(){
    // evaluate integrals along the array in s
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); j++) {
        switch (isospin) {
            case 1:
                // for T-> 3 pi
                angular_average_list.push_back(eval_int(mandelstam_list.interval_s_th[j]));
                break;
                        
            default:
                throw std::domain_error{"Wrong isospin for the chosen Process!"};
        }
    }
}


void AngularAverageEta3Pi::build_matching_list(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        real_angular_average_list.push_back(std::real(angular_average_list[j]));
        imag_angular_average_list.push_back(std::imag(angular_average_list[j]));
    }
    return;
}

void AngularAverageEta3PiC::build_matching_list(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        real_angular_average_list.push_back(std::real(angular_average_list[j]));
        imag_angular_average_list.push_back(std::imag(angular_average_list[j]));
    }
    return;
}

void AngularAverageEtap3Pi::build_matching_list(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        real_angular_average_list.push_back(std::real(angular_average_list[j]));
        imag_angular_average_list.push_back(std::imag(angular_average_list[j]));
    }
    return;
}

void AngularAverageEtap3PiC::build_matching_list(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        real_angular_average_list.push_back(std::real(angular_average_list[j]));
        imag_angular_average_list.push_back(std::imag(angular_average_list[j]));
    }
    return;
}

void AngularAverageEtapEtaPiPi::build_matching_list(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        real_angular_average_list.push_back(std::real(angular_average_list[j]));
        imag_angular_average_list.push_back(std::imag(angular_average_list[j]));
    }
    for (std::size_t j=0; j<angular_average_intermediate_list.size(); ++j) {
        real_angular_average_intermediate_list.push_back(std::real(angular_average_intermediate_list[j]));
        imag_angular_average_intermediate_list.push_back(std::imag(angular_average_intermediate_list[j]));
    }
    return;
}

void AngularAverageEtapEtaPiPiC::build_matching_list(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        real_angular_average_list.push_back(std::real(angular_average_list[j]));
        imag_angular_average_list.push_back(std::imag(angular_average_list[j]));
    }
    for (std::size_t j=0; j<angular_average_intermediate_list.size(); ++j) {
        real_angular_average_intermediate_list.push_back(std::real(angular_average_intermediate_list[j]));
        imag_angular_average_intermediate_list.push_back(std::imag(angular_average_intermediate_list[j]));
    }
    return;
}

void AngularAverageV3Pi::build_matching_list(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        real_angular_average_list.push_back(std::real(angular_average_list[j]));
        imag_angular_average_list.push_back(std::imag(angular_average_list[j]));
    }
    return;
}

void AngularAverageX3Pi::build_matching_list(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        real_angular_average_list.push_back(std::real(angular_average_list[j]));
        imag_angular_average_list.push_back(std::imag(angular_average_list[j]));
    }
    return;
}

void AngularAverageT3Pi::build_matching_list(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        real_angular_average_list.push_back(std::real(angular_average_list[j]));
        imag_angular_average_list.push_back(std::imag(angular_average_list[j]));
    }
    return;
}
//------------------------------------------------------------------------
//------------------------------------------------------------------------




//------------------------------------------------------------------------
//-- Interpolation -------------------------------------------------------
//------------------------------------------------------------------------
Complex AngularAverageEta3Pi::angular_average(double s)const{
    return spline_angular_average(s);
}

Complex AngularAverageEta3PiC::angular_average(double s)const{
    return spline_angular_average(s);
}

Complex AngularAverageEtap3Pi::angular_average(double s)const{
    return spline_angular_average(s);
}

Complex AngularAverageEtap3PiC::angular_average(double s)const{
    return spline_angular_average(s);
}

Complex AngularAverageEtapEtaPiPi::angular_average(double s)const{
    switch (isospin) {
        // for pi-pi intermediate states
        case 000: case 111:
            return spline_angular_average(s);
        // for eta-pi intermediate states
        case 010: case 110:
            return spline_angular_average_intermediate(s);
        default:
            throw std::domain_error{"Unvalid isospin for the chosen configuration!"};
    }
}

Complex AngularAverageEtapEtaPiPiC::angular_average(double s)const{
    switch (isospin) {
        // for pi-pi intermediate states
        case 000: case 111:
            return spline_angular_average(s);
        // for eta-pi intermediate states
        case 010: case 110:
            return spline_angular_average_intermediate(s);
        default:
            throw std::domain_error{"Unvalid isospin for the chosen configuration!"};
    }
}

Complex AngularAverageV3Pi::angular_average(double s)const{
    return spline_angular_average(s);
}

Complex AngularAverageX3Pi::angular_average(double s)const{
    return spline_angular_average(s);
}

Complex AngularAverageT3Pi::angular_average(double s)const{
    return spline_angular_average(s);
}
//------------------------------------------------------------------------
//------------------------------------------------------------------------


Complex AngularAverage::kallen(double s, double mass1, double mass2)const{
    return std::pow(s,2.)+ std::pow(mass1,4.) +std::pow(mass2,4.)
    -2.*(s*std::pow(mass1,2.) +s*std::pow(mass2,2.) +std::pow(mass1,2.)*std::pow(mass2,2.));
}


Complex AngularAverage::nu(double x)const{
    double x_I;
    double x_IV;
    double x_0;
    switch (isospin) {
        case 010: case 110:
            x_I=t_I;
            x_IV=t_IV;
            x_0=t_x;
            break;
            
        default:
            x_I=s_I;
            x_IV=s_IV;
            x_0=s_x;
            break;
    }
    if (x<=x_IV) {
        // within the matching range use the following two if conditions
        if (x>=x_I-validity && x<=x_I+validity) {
            return std::sqrt(Complex(x_IV-x))/std::sqrt(Complex(x)) * std::sqrt(Complex(1.-x_0/x));
        }
        if (x>=x_IV-validity) {
            return std::sqrt(Complex(1.-x_I/x)) * std::sqrt(Complex(1.-x_0/x));
        }
        return std::sqrt(Complex(1.-x_I/x))*std::sqrt(Complex(x_IV-x)) * std::sqrt(Complex(1.-x_0/x));
    }
    if (x>x_IV) {
        // within the matching range use the following if condition
        if (x<=x_IV+validity) {
            return 1.i*std::sqrt(Complex(1.-x_I/x))
            * std::sqrt(Complex(1.-x_0/x))
            ;
        }
        return 1.i*std::sqrt(Complex(1.-x_I/x))*std::sqrt(Complex(x-x_IV))
        * std::sqrt(Complex(1.-x_0/x))
        ;
    }
    throw std::domain_error{"Value for Mandelstam 's' not allowed 3!"};
}



// singularities in the angular average at the thresholds are already analyticaly canceled with the
// corresponding factors that would appear in 'nu(s)'
void AngularAverageEta3Pi::build_matched_tilde(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        if (isospin==0 || isospin==2) {
            matched_angular_average.push_back(
                                (matched.s_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                + 1.i*matched.s_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                )/nu(mandelstam_list.interval_s_th[j]));
        }
        if (isospin==1) {
            matched_angular_average.push_back(
                                (matched.p_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                + 1.i*matched.p_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                )/std::pow(nu(mandelstam_list.interval_s_th[j]),3.0));
        }
    }
    return;
}

void AngularAverageEta3PiC::build_matched_tilde(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        if (isospin==0 || isospin==2) {
            matched_angular_average.push_back(
                                (matched.s_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                + 1.i*matched.s_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                )/nu(mandelstam_list.interval_s_th[j]));
        }
        if (isospin==1) {
            matched_angular_average.push_back(
                                (matched.p_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                + 1.i*matched.p_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                )/std::pow(nu(mandelstam_list.interval_s_th[j]),3.0));
        }
    }
    return;
}

void AngularAverageEtap3Pi::build_matched_tilde(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        if (isospin==0 || isospin==2) {
            matched_angular_average.push_back(
                                (matched.s_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                + 1.i*matched.s_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                )/nu(mandelstam_list.interval_s_th[j]));
        }
        if (isospin==1) {
            matched_angular_average.push_back(
                                (matched.p_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                + 1.i*matched.p_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                )/std::pow(nu(mandelstam_list.interval_s_th[j]),3.0));
        }
    }
    return;
}

void AngularAverageEtap3PiC::build_matched_tilde(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        if (isospin==0 || isospin==2) {
            matched_angular_average.push_back(
                                (matched.s_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                + 1.i*matched.s_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                )/nu(mandelstam_list.interval_s_th[j]));
        }
        if (isospin==1) {
            matched_angular_average.push_back(
                                (matched.p_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                + 1.i*matched.p_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                )/std::pow(nu(mandelstam_list.interval_s_th[j]),3.0));
        }
    }
    return;
}

//Caution!:
//for this Process (eta->3pi) Bose symmetry demands that isospin 0 and 2 states are in an S-wave,
//whereas isospin 1 states are in a P-wave
void AngularAverageEtapEtaPiPi::build_matched_tilde(){
    switch (isospin) {
        case 000:
            for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
                matched_angular_average.push_back(
                                                  (matched.s_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                                  + 1.i*matched.s_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                           )/nu(mandelstam_list.interval_s_th[j])
                                          );
            }
            // this list is not needed for this isospin, therefore fill it with arbitrary numbers
            // (otherwise the interpolation in the constructor fails)
            for (std::size_t j=0; j<mandelstam_list.interval_t_th.size(); ++j) {
                matched_angular_average_intermediate.push_back(100.);
            }
            break;
            
        case 111:
            for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
                matched_angular_average.push_back(
                                                  (matched.p_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                                  + 1.i*matched.p_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                           )/std::pow(nu(mandelstam_list.interval_s_th[j]),3.)
                                          );
            }
            // this list is not needed for this isospin, therefore fill it with arbitrary numbers
            // (otherwise the interpolation in the constructor fails)
            for (std::size_t j=0; j<mandelstam_list.interval_t_th.size(); ++j) {
                matched_angular_average_intermediate.push_back(100.);
            }
            break;
            
        case 010: case 110:
            for (std::size_t j=0; j<mandelstam_list.interval_t_th.size(); ++j) {
                matched_angular_average_intermediate.push_back(
                                           (matched_t_channel.s_match(mandelstam_list.interval_t_th[j], real_angular_average_intermediate)
                                          + 1.i*matched_t_channel.s_match(mandelstam_list.interval_t_th[j], imag_angular_average_intermediate)
                                           )
                                            /nu(mandelstam_list.interval_t_th[j])
                                          );
                
            }
            // this list is not needed for this isospin, therefore fill it with arbitrary numbers
            // (otherwise the interpolation in the constructor fails)
            for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
                matched_angular_average.push_back(100.);
            }
        break;
    }
    return;
}

void AngularAverageEtapEtaPiPiC::build_matched_tilde(){
    switch (isospin) {
        case 000:
            for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
                matched_angular_average.push_back(
                                                  (matched.s_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                                  + 1.i*matched.s_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                           )/nu(mandelstam_list.interval_s_th[j])
                                          );
            }
            // this list is not needed for this isospin, therefore fill it with arbitrary numbers
            // (otherwise the interpolation in the constructor fails)
            for (std::size_t j=0; j<mandelstam_list.interval_t_th.size(); ++j) {
                matched_angular_average_intermediate.push_back(100.);
            }
            break;
            
        case 111:
            for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
                matched_angular_average.push_back(
                                                  (matched.p_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                                  + 1.i*matched.p_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                           )/std::pow(nu(mandelstam_list.interval_s_th[j]),3.)
                                          );
            }
            // this list is not needed for this isospin, therefore fill it with arbitrary numbers
            // (otherwise the interpolation in the constructor fails)
            for (std::size_t j=0; j<mandelstam_list.interval_t_th.size(); ++j) {
                matched_angular_average_intermediate.push_back(100.);
            }
            break;
            
        case 010: case 110:
            for (std::size_t j=0; j<mandelstam_list.interval_t_th.size(); ++j) {
                matched_angular_average_intermediate.push_back(
                                           (matched_t_channel.s_match(mandelstam_list.interval_t_th[j], real_angular_average_intermediate)
                                          + 1.i*matched_t_channel.s_match(mandelstam_list.interval_t_th[j], imag_angular_average_intermediate)
                                           )
                                            /nu(mandelstam_list.interval_t_th[j])
                                          );
                
            }
            // this list is not needed for this isospin, therefore fill it with arbitrary numbers
            // (otherwise the interpolation in the constructor fails)
            for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
                matched_angular_average.push_back(100.);
            }
        break;
    }
    return;
}

void AngularAverageV3Pi::build_matched_tilde(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        if (isospin==1) {
            matched_angular_average.push_back(
                                (matched.p_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                + 1.i*matched.p_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                )/std::pow(nu(mandelstam_list.interval_s_th[j]),3.0));
        }
    }
    return;
}

void AngularAverageX3Pi::build_matched_tilde(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        if (isospin==1) {
            matched_angular_average.push_back(
                                (matched.p_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                + 1.i*matched.p_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                )/std::pow(nu(mandelstam_list.interval_s_th[j]),3.0));
        }
    }
    return;
}

void AngularAverageT3Pi::build_matched_tilde(){
    for (std::size_t j=0; j<mandelstam_list.interval_s_th.size(); ++j) {
        if (isospin==1) {
            matched_angular_average.push_back(
                                (matched.d_match(mandelstam_list.interval_s_th[j], real_angular_average)
                                + 1.i*matched.d_match(mandelstam_list.interval_s_th[j], imag_angular_average)
                                )/std::pow(nu(mandelstam_list.interval_s_th[j]),5.0));
        }
    }
    return;
}

void AngularAverageEta3Pi::make_splines(){
    real_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, real_angular_average_list,gsl::InterpolationMethod::cubic);
    imag_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, imag_angular_average_list,gsl::InterpolationMethod::cubic);
    spline_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th,angular_average_list,gsl::InterpolationMethod::cubic);
}

void AngularAverageEta3PiC::make_splines(){
    real_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, real_angular_average_list,gsl::InterpolationMethod::cubic);
    imag_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, imag_angular_average_list,gsl::InterpolationMethod::cubic);
    spline_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th,angular_average_list,gsl::InterpolationMethod::cubic);
}

void AngularAverageEtap3Pi::make_splines(){
    real_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, real_angular_average_list,gsl::InterpolationMethod::cubic);
    imag_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, imag_angular_average_list,gsl::InterpolationMethod::cubic);
    spline_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th,angular_average_list,gsl::InterpolationMethod::cubic);
}

void AngularAverageEtap3PiC::make_splines(){
    real_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, real_angular_average_list,gsl::InterpolationMethod::cubic);
    imag_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, imag_angular_average_list,gsl::InterpolationMethod::cubic);
    spline_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th,angular_average_list,gsl::InterpolationMethod::cubic);
}

void AngularAverageEtapEtaPiPi::make_splines(){
    real_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, real_angular_average_list,gsl::InterpolationMethod::cubic);
    imag_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, imag_angular_average_list,gsl::InterpolationMethod::cubic);
    spline_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th,angular_average_list,gsl::InterpolationMethod::cubic);
    real_angular_average_intermediate = gsl::Interpolate(mandelstam_list.interval_t_th, real_angular_average_intermediate_list,gsl::InterpolationMethod::cubic);
    imag_angular_average_intermediate = gsl::Interpolate(mandelstam_list.interval_t_th, imag_angular_average_intermediate_list,gsl::InterpolationMethod::cubic);
    spline_angular_average_intermediate = cauchy::Interpolate(mandelstam_list.interval_t_th,angular_average_intermediate_list,gsl::InterpolationMethod::cubic);
}

void AngularAverageEtapEtaPiPiC::make_splines(){
    real_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, real_angular_average_list,gsl::InterpolationMethod::cubic);
    imag_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, imag_angular_average_list,gsl::InterpolationMethod::cubic);
    spline_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th,angular_average_list,gsl::InterpolationMethod::cubic);
    real_angular_average_intermediate = gsl::Interpolate(mandelstam_list.interval_t_th, real_angular_average_intermediate_list,gsl::InterpolationMethod::cubic);
    imag_angular_average_intermediate = gsl::Interpolate(mandelstam_list.interval_t_th, imag_angular_average_intermediate_list,gsl::InterpolationMethod::cubic);
    spline_angular_average_intermediate = cauchy::Interpolate(mandelstam_list.interval_t_th,angular_average_intermediate_list,gsl::InterpolationMethod::cubic);
}

void AngularAverageV3Pi::make_splines(){
    real_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, real_angular_average_list,gsl::InterpolationMethod::cubic);
    imag_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, imag_angular_average_list,gsl::InterpolationMethod::cubic);
    spline_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th,angular_average_list,gsl::InterpolationMethod::cubic);
}

void AngularAverageX3Pi::make_splines(){
    real_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, real_angular_average_list,gsl::InterpolationMethod::cubic);
    imag_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, imag_angular_average_list,gsl::InterpolationMethod::cubic);
    spline_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th,angular_average_list,gsl::InterpolationMethod::cubic);
}

void AngularAverageT3Pi::make_splines(){
    real_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, real_angular_average_list,gsl::InterpolationMethod::cubic);
    imag_angular_average = gsl::Interpolate(mandelstam_list.interval_s_th, imag_angular_average_list,gsl::InterpolationMethod::cubic);
    spline_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th,angular_average_list,gsl::InterpolationMethod::cubic);
}

void AngularAverageEta3Pi::make_splines_matched(){
    spline_matched_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th, matched_angular_average, gsl::InterpolationMethod::cubic);
}

void AngularAverageEta3PiC::make_splines_matched(){
    spline_matched_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th, matched_angular_average, gsl::InterpolationMethod::cubic);
}

void AngularAverageEtap3Pi::make_splines_matched(){
    spline_matched_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th, matched_angular_average, gsl::InterpolationMethod::cubic);
}

void AngularAverageEtap3PiC::make_splines_matched(){
    spline_matched_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th, matched_angular_average, gsl::InterpolationMethod::cubic);
}

void AngularAverageEtapEtaPiPi::make_splines_matched(){
    spline_matched_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th, matched_angular_average, gsl::InterpolationMethod::cubic);
    spline_matched_angular_average_intermediate = cauchy::Interpolate(mandelstam_list.interval_t_th, matched_angular_average_intermediate, gsl::InterpolationMethod::cubic);
}

void AngularAverageEtapEtaPiPiC::make_splines_matched(){
    spline_matched_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th, matched_angular_average, gsl::InterpolationMethod::cubic);
    spline_matched_angular_average_intermediate = cauchy::Interpolate(mandelstam_list.interval_t_th, matched_angular_average_intermediate, gsl::InterpolationMethod::cubic);
}

void AngularAverageV3Pi::make_splines_matched(){
    spline_matched_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th, matched_angular_average, gsl::InterpolationMethod::cubic);
}

void AngularAverageX3Pi::make_splines_matched(){
    spline_matched_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th, matched_angular_average, gsl::InterpolationMethod::cubic);
}

void AngularAverageT3Pi::make_splines_matched(){
    spline_matched_angular_average = cauchy::Interpolate(mandelstam_list.interval_s_th, matched_angular_average, gsl::InterpolationMethod::cubic);
}



Complex AngularAverageEta3Pi::operator()(double s)const{
    return spline_matched_angular_average(s);
}

Complex AngularAverageEta3PiC::operator()(double s)const{
    return spline_matched_angular_average(s);
}

Complex AngularAverageEtap3Pi::operator()(double s)const{
    return spline_matched_angular_average(s);
}

Complex AngularAverageEtap3PiC::operator()(double s)const{
    return spline_matched_angular_average(s);
}

Complex AngularAverageEtapEtaPiPi::operator()(double s)const{
    switch (isospin) {
        case 000: case 111:
            return spline_matched_angular_average(s);
            
        case 010: case 110:
            return spline_matched_angular_average_intermediate(s);
            
        default:
            throw std::domain_error{"Unvalid isospin for the chosen configuration!"};
    }  
}

Complex AngularAverageEtapEtaPiPiC::operator()(double s)const{
    switch (isospin) {
        case 000: case 111:
            return spline_matched_angular_average(s);
            
        case 010: case 110:
            return spline_matched_angular_average_intermediate(s);
            
        default:
            throw std::domain_error{"Unvalid isospin for the chosen configuration!"};
    }  
}

Complex AngularAverageV3Pi::operator()(double s)const{
    return spline_matched_angular_average(s);
}

Complex AngularAverageX3Pi::operator()(double s)const{
    return spline_matched_angular_average(s);
}

Complex AngularAverageT3Pi::operator()(double s)const{
    return spline_matched_angular_average(s);
}