#include "splined_omnes.h"
#include <fstream>
using omnes::Omnes;
using namespace splined_omnes;
using namespace std::complex_literals;
using enums::Setting;




//-- Constructors ---------------------------------------------------------------------

//interpolate the Omnès functions
OmnesSpline::OmnesSpline(double mass_1, double mass_2, double mass_3, double mass_4, double cutoff)
:
s_I{pow(mass_1+mass_2,2.)},
t_I{pow(mass_1+mass_3,2.)},
polar_egg(mass_4,mass_1),
gamma(mass_1,mass_2,mass_3,mass_4)
{}


OmnesSplineEta3Pi::OmnesSplineEta3Pi(double mass_1, double mass_4, double cutoff, phase::PhaseEta3Pi phases)
:
OmnesSpline{mass_1,mass_1,mass_1,mass_4,cutoff},
grid(mass_1,mass_4,cutoff),
ph{phases}
{build_omnes(); spline_omnes();}

OmnesSplineEtap3Pi::OmnesSplineEtap3Pi(double mass_1, double mass_4, double cutoff, phase::PhaseEtap3Pi phases)
:
OmnesSpline{mass_1,mass_1,mass_1,mass_4,cutoff},
grid(mass_1,mass_4,cutoff),
ph{phases}
{build_omnes(); spline_omnes();}

OmnesSplineEtapEtaPiPi::OmnesSplineEtapEtaPiPi(double mass_1, double mass_3, double mass_4, double cutoff, phase::PhaseEtapEtaPiPi phases)
:
OmnesSpline{mass_1,mass_1,mass_1,mass_4,cutoff},
grid(mass_1,mass_3,mass_4,cutoff),
ph{phases}
{build_omnes(); spline_omnes();}

OmnesSplineV3Pi::OmnesSplineV3Pi(double mass_1, double mass_4, double cutoff, phase::PhaseV3Pi phases)
:
OmnesSpline{mass_1,mass_1,mass_1,mass_4,cutoff},
grid(mass_1,mass_4,cutoff),
ph{phases}
{build_omnes(); spline_omnes();}

OmnesSplineX3Pi::OmnesSplineX3Pi(double mass_1, double mass_4, double cutoff, phase::PhaseX3Pi phases)
:
OmnesSpline{mass_1,mass_1,mass_1,mass_4,cutoff},
grid(mass_1,mass_4,cutoff),
ph{phases}
{build_omnes(); spline_omnes();}

OmnesSplineT3Pi::OmnesSplineT3Pi(double mass_1, double mass_4, double cutoff, phase::PhaseT3Pi phases)
:
OmnesSpline{mass_1,mass_1,mass_1,mass_4,cutoff},
grid(mass_1,mass_4,cutoff),
ph{phases}
{build_omnes(); spline_omnes();}

//-- Lists for Omnes Function with complex arguments ----------------------
// need to introduce lambda-expressions for the phases, otherwise they can not be used
// to initialize the constructors of the corresponding Omnès functions


void OmnesSplineEta3Pi::build_omnes(){
    
    std::cout<<"initializing omnes functions and phases...\n";
    const auto phase0{[this](double s){ return ph.phase(s,0); }};
    const auto phase1{[this](double s){ return ph.phase(s,1); }};
    const auto phase2{[this](double s){ return ph.phase(s,2); }};
    
    Omnes omn0(phase0,s_I,infinitesimal);
    Omnes omn1(phase1,s_I,infinitesimal);
    Omnes omn2(phase2,s_I,infinitesimal);
    
    
    std::cout<<"for pi pi intermediate states in eta -> 3pi\n";
    // pi-pi intermediate states: infinitesimal above and below cut
    for (std::size_t i=0; i<grid.interval_s.size(); ++i) {
        omnes_0_above.push_back(omn0(grid.interval_s[i]+1.i*infinitesimal));
        omnes_1_above.push_back(omn1(grid.interval_s[i]+1.i*infinitesimal));
        omnes_2_above.push_back(omn2(grid.interval_s[i]+1.i*infinitesimal));
        
        omnes_0_below.push_back(omn0(grid.interval_s[i]-1.i*infinitesimal));
        omnes_1_below.push_back(omn1(grid.interval_s[i]-1.i*infinitesimal));
        omnes_2_below.push_back(omn2(grid.interval_s[i]-1.i*infinitesimal));
    }

    // pi-pi intermediate states: egg-like contour for three-equal masses (e.g. like in eta->3pi)
    // expressed in terms of curve parameter phi
    for (std::size_t i=0; i<grid.interval_phi.size(); ++i) {
        omnes_0_egg.push_back(omn0(polar_egg(grid.interval_phi[i])));
        omnes_1_egg.push_back(omn1(polar_egg(grid.interval_phi[i])));
        omnes_2_egg.push_back(omn2(polar_egg(grid.interval_phi[i])));
    }
    std::cout<<"Done!\n";
    return;
}

void OmnesSplineEtap3Pi::build_omnes(){
    
    std::cout<<"initializing omnes functions and phases...\n";
    const auto phase0{[this](double s){ return ph.phase(s,0); }};
    const auto phase1{[this](double s){ return ph.phase(s,1); }};
    const auto phase2{[this](double s){ return ph.phase(s,2); }};
    
    Omnes omn0(phase0,s_I,infinitesimal);
    Omnes omn1(phase1,s_I,infinitesimal);
    Omnes omn2(phase2,s_I,infinitesimal);
    
    
    std::cout<<"for pi pi intermediate states in eta' -> 3pi\n";
    // pi-pi intermediate states: infinitesimal above and below cut
    for (std::size_t i=0; i<grid.interval_s.size(); ++i) {
        omnes_0_above.push_back(omn0(grid.interval_s[i]+1.i*infinitesimal));
        omnes_1_above.push_back(omn1(grid.interval_s[i]+1.i*infinitesimal));
        omnes_2_above.push_back(omn2(grid.interval_s[i]+1.i*infinitesimal));
        
        omnes_0_below.push_back(omn0(grid.interval_s[i]-1.i*infinitesimal));
        omnes_1_below.push_back(omn1(grid.interval_s[i]-1.i*infinitesimal));
        omnes_2_below.push_back(omn2(grid.interval_s[i]-1.i*infinitesimal));
    }

    // pi-pi intermediate states: egg-like contour for three-equal masses (e.g. like in eta->3pi)
    // expressed in terms of curve parameter phi
    for (std::size_t i=0; i<grid.interval_phi.size(); ++i) {
        omnes_0_egg.push_back(omn0(polar_egg(grid.interval_phi[i])));
        omnes_1_egg.push_back(omn1(polar_egg(grid.interval_phi[i])));
        omnes_2_egg.push_back(omn2(polar_egg(grid.interval_phi[i])));
    }
    std::cout<<"Done!\n";
    return;
}

void OmnesSplineEtapEtaPiPi::build_omnes(){
    
    std::cout<<"initializing omnes functions and phases...\n";
    const auto phase0{[this](double s){ return ph.phase(s,0); }};
    const auto phase1{[this](double s){ return ph.phase(s,1); }};
    const auto phase2{[this](double s){ return ph.phase(s,2); }};
    const auto ph_eta_pi{[this](double s){ return ph.phase_eta_pi(s); }};
    
    Omnes omn0(phase0,s_I,infinitesimal);
    Omnes omn1(phase1,s_I,infinitesimal);
    Omnes omn2(phase2,s_I,infinitesimal);
    Omnes omn_eta_pi(ph_eta_pi,t_I,infinitesimal);
    
    
    std::cout<<"for pi pi intermediate states in eta' -> eta pi pi\n";
    // pi-pi intermediate states: infinitesimal above and below cut
    for (std::size_t i=0; i<grid.interval_s.size(); ++i) {
        omnes_0_above.push_back(omn0(grid.interval_s[i]+1.i*infinitesimal));
        omnes_1_above.push_back(omn1(grid.interval_s[i]+1.i*infinitesimal));
        omnes_2_above.push_back(omn2(grid.interval_s[i]+1.i*infinitesimal));
        
        omnes_0_below.push_back(omn0(grid.interval_s[i]-1.i*infinitesimal));
        omnes_1_below.push_back(omn1(grid.interval_s[i]-1.i*infinitesimal));
        omnes_2_below.push_back(omn2(grid.interval_s[i]-1.i*infinitesimal));
    }
    
    std::cout<<"for eta pi intermediate states in eta' -> eta pi pi\n";
    // eta-pi intermediate state (S-wave): infinitesimal above and below cut
    for (std::size_t i=0; i<grid.interval_t.size(); ++i) {
        omnes_eta_pi_above.push_back(omn_eta_pi(grid.interval_t[i]+1.i*infinitesimal));
        omnes_eta_pi_below.push_back(omn_eta_pi(grid.interval_t[i]-1.i*infinitesimal));
    }
    
    std::cout<<"along egg-like contour in eta' -> eta pi pi\n";
    // Evaluate Omnès along the curve parameters in eta'-> eta pi pi
    // take into account the cases for the three different angular averages (plus, minus, zero),
    // the ones for the upper and lower endpoints of Mandelstam t within the Pinocchio path  (upper, lower)
    // as well as the two different integration regions within the Pinocchio path to avoid end-point-singularities (a, b)
    for (std::size_t i=0; i<grid.interval_y.size(); ++i) {
        /// express in terms of curve parameter y
        double y= grid.interval_y[i];
        // eta-pi intermediate states
        omnes_eta_pi_plus_upper_a.push_back(omn_eta_pi(gamma(y,Setting::plus_upper_a)));
        omnes_eta_pi_plus_upper_b.push_back(omn_eta_pi(gamma(y,Setting::plus_upper_b)));
        omnes_eta_pi_plus_lower_a.push_back(omn_eta_pi(gamma(y,Setting::plus_lower_a)));
        omnes_eta_pi_plus_lower_b.push_back(omn_eta_pi(gamma(y,Setting::plus_lower_b)));

        omnes_eta_pi_zero_upper_a.push_back(omn_eta_pi(gamma(y,Setting::zero_upper_a)));
        omnes_eta_pi_zero_upper_b.push_back(omn_eta_pi(gamma(y,Setting::zero_upper_b)));
        omnes_eta_pi_zero_lower_a.push_back(omn_eta_pi(gamma(y,Setting::zero_lower_a)));
        omnes_eta_pi_zero_lower_b.push_back(omn_eta_pi(gamma(y,Setting::zero_lower_b)));


        // pi-pi intermediate states
        omnes_0_minus_upper_a.push_back(omn0(gamma(y,Setting::minus_upper_a)));
        omnes_0_minus_upper_b.push_back(omn0(gamma(y,Setting::minus_upper_b)));
        omnes_0_minus_lower_a.push_back(omn0(gamma(y,Setting::minus_lower_a)));
        omnes_0_minus_lower_b.push_back(omn0(gamma(y,Setting::minus_lower_b)));

        ///For C-violation
        omnes_1_minus_upper_a.push_back(omn1(gamma(y,Setting::minus_upper_a)));
        omnes_1_minus_upper_b.push_back(omn1(gamma(y,Setting::minus_upper_b)));
        omnes_1_minus_lower_a.push_back(omn1(gamma(y,Setting::minus_lower_a)));
        omnes_1_minus_lower_b.push_back(omn1(gamma(y,Setting::minus_lower_b)));
    }
    std::cout<<"Done!\n";
    return;
}

void OmnesSplineV3Pi::build_omnes(){
    
    std::cout<<"initializing omnes functions and phases...\n";
    const auto phase{[this](double s){ return ph.phase(s,1); }};
    
    Omnes omn(phase,s_I,infinitesimal);
    
    std::cout<<"for pi pi intermediate states in V -> 3pi\n";
    // pi-pi intermediate states: infinitesimal above and below cut
    for (std::size_t i=0; i<grid.interval_s.size(); ++i) {
        omnes_above.push_back(omn(grid.interval_s[i]+1.i*infinitesimal));
        
        omnes_below.push_back(omn(grid.interval_s[i]-1.i*infinitesimal));
    }

    // pi-pi intermediate states: egg-like contour for three-equal masses (e.g. like in eta->3pi)
    // expressed in terms of curve parameter phi
    for (std::size_t i=0; i<grid.interval_phi.size(); ++i) {
        omnes_egg.push_back(omn(polar_egg(grid.interval_phi[i])));
    }
    std::cout<<"Done!\n";
    return;
}

void OmnesSplineX3Pi::build_omnes(){
    
    std::cout<<"initializing omnes functions and phases...\n";
    const auto phase{[this](double s){ return ph.phase(s,1); }};
    
    Omnes omn(phase,s_I,infinitesimal);
    
    std::cout<<"for pi pi intermediate states in X -> 3pi\n";
    // pi-pi intermediate states: infinitesimal above and below cut
    for (std::size_t i=0; i<grid.interval_s.size(); ++i) {
        omnes_above.push_back(omn(grid.interval_s[i]+1.i*infinitesimal));
        
        omnes_below.push_back(omn(grid.interval_s[i]-1.i*infinitesimal));
    }

    // pi-pi intermediate states: egg-like contour for three-equal masses (e.g. like in eta->3pi)
    // expressed in terms of curve parameter phi
    for (std::size_t i=0; i<grid.interval_phi.size(); ++i) {
        omnes_egg.push_back(omn(polar_egg(grid.interval_phi[i])));
    }
    std::cout<<"Done!\n";
    return;
}

void OmnesSplineT3Pi::build_omnes(){
    
    std::cout<<"initializing omnes functions and phases...\n";
    const auto phase{[this](double s){ return ph.phase(s,1); }};
    
    Omnes omn(phase,s_I,infinitesimal);
    
    std::cout<<"for pi pi intermediate states in T -> 3pi\n";
    // pi-pi intermediate states: infinitesimal above and below cut
    for (std::size_t i=0; i<grid.interval_s.size(); ++i) {
        omnes_above.push_back(omn(grid.interval_s[i]+1.i*infinitesimal));
        
        omnes_below.push_back(omn(grid.interval_s[i]-1.i*infinitesimal));
    }

    // pi-pi intermediate states: egg-like contour for three-equal masses (e.g. like in eta->3pi)
    // expressed in terms of curve parameter phi
    for (std::size_t i=0; i<grid.interval_phi.size(); ++i) {
        omnes_egg.push_back(omn(polar_egg(grid.interval_phi[i])));
    }
    std::cout<<"Done!\n";
    return;
}

void OmnesSplineEta3Pi::spline_omnes(){
    spline_omnes_0_above = cauchy::Interpolate(grid.interval_s,omnes_0_above,gsl::InterpolationMethod::cubic);
    spline_omnes_1_above = cauchy::Interpolate(grid.interval_s,omnes_1_above,gsl::InterpolationMethod::cubic);
    spline_omnes_2_above = cauchy::Interpolate(grid.interval_s,omnes_2_above,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_0_below = cauchy::Interpolate(grid.interval_s,omnes_0_below,gsl::InterpolationMethod::cubic);
    spline_omnes_1_below = cauchy::Interpolate(grid.interval_s,omnes_1_below,gsl::InterpolationMethod::cubic);
    spline_omnes_2_below = cauchy::Interpolate(grid.interval_s,omnes_2_below,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_0_egg = cauchy::Interpolate(grid.interval_phi,omnes_0_egg,gsl::InterpolationMethod::cubic);
    spline_omnes_1_egg = cauchy::Interpolate(grid.interval_phi,omnes_1_egg,gsl::InterpolationMethod::cubic);
    spline_omnes_2_egg = cauchy::Interpolate(grid.interval_phi,omnes_2_egg,gsl::InterpolationMethod::cubic);
    return;
}

void OmnesSplineEtap3Pi::spline_omnes(){
    spline_omnes_0_above = cauchy::Interpolate(grid.interval_s,omnes_0_above,gsl::InterpolationMethod::cubic);
    spline_omnes_1_above = cauchy::Interpolate(grid.interval_s,omnes_1_above,gsl::InterpolationMethod::cubic);
    spline_omnes_2_above = cauchy::Interpolate(grid.interval_s,omnes_2_above,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_0_below = cauchy::Interpolate(grid.interval_s,omnes_0_below,gsl::InterpolationMethod::cubic);
    spline_omnes_1_below = cauchy::Interpolate(grid.interval_s,omnes_1_below,gsl::InterpolationMethod::cubic);
    spline_omnes_2_below = cauchy::Interpolate(grid.interval_s,omnes_2_below,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_0_egg = cauchy::Interpolate(grid.interval_phi,omnes_0_egg,gsl::InterpolationMethod::cubic);
    spline_omnes_1_egg = cauchy::Interpolate(grid.interval_phi,omnes_1_egg,gsl::InterpolationMethod::cubic);
    spline_omnes_2_egg = cauchy::Interpolate(grid.interval_phi,omnes_2_egg,gsl::InterpolationMethod::cubic);
    return;
}

void OmnesSplineEtapEtaPiPi::spline_omnes(){
    spline_omnes_0_above = cauchy::Interpolate(grid.interval_s,omnes_0_above,gsl::InterpolationMethod::cubic);
    spline_omnes_1_above = cauchy::Interpolate(grid.interval_s,omnes_1_above,gsl::InterpolationMethod::cubic);
    spline_omnes_2_above = cauchy::Interpolate(grid.interval_s,omnes_2_above,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_0_below = cauchy::Interpolate(grid.interval_s,omnes_0_below,gsl::InterpolationMethod::cubic);
    spline_omnes_1_below = cauchy::Interpolate(grid.interval_s,omnes_1_below,gsl::InterpolationMethod::cubic);
    spline_omnes_2_below = cauchy::Interpolate(grid.interval_s,omnes_2_below,gsl::InterpolationMethod::cubic);
    //
    //-- eta'-> eta pi pi ---------
    spline_omnes_eta_pi_above = cauchy::Interpolate(grid.interval_t,omnes_eta_pi_above,gsl::InterpolationMethod::cubic);
    spline_omnes_eta_pi_below = cauchy::Interpolate(grid.interval_t,omnes_eta_pi_below,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_eta_pi_plus_upper_a = cauchy::Interpolate(grid.interval_y,omnes_eta_pi_plus_upper_a,gsl::InterpolationMethod::cubic);
    spline_omnes_eta_pi_plus_upper_b = cauchy::Interpolate(grid.interval_y,omnes_eta_pi_plus_upper_b,gsl::InterpolationMethod::cubic);
    spline_omnes_eta_pi_plus_lower_a = cauchy::Interpolate(grid.interval_y,omnes_eta_pi_plus_lower_a,gsl::InterpolationMethod::cubic);
    spline_omnes_eta_pi_plus_lower_b = cauchy::Interpolate(grid.interval_y,omnes_eta_pi_plus_lower_b,gsl::InterpolationMethod::cubic);
    spline_omnes_eta_pi_zero_upper_a = cauchy::Interpolate(grid.interval_y,omnes_eta_pi_zero_upper_a,gsl::InterpolationMethod::cubic);
    spline_omnes_eta_pi_zero_upper_b = cauchy::Interpolate(grid.interval_y,omnes_eta_pi_zero_upper_b,gsl::InterpolationMethod::cubic);
    spline_omnes_eta_pi_zero_lower_a = cauchy::Interpolate(grid.interval_y,omnes_eta_pi_zero_lower_a,gsl::InterpolationMethod::cubic);
    spline_omnes_eta_pi_zero_lower_b = cauchy::Interpolate(grid.interval_y,omnes_eta_pi_zero_lower_b,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_0_minus_upper_a = cauchy::Interpolate(grid.interval_y,omnes_0_minus_upper_a,gsl::InterpolationMethod::cubic);
    spline_omnes_0_minus_upper_b = cauchy::Interpolate(grid.interval_y,omnes_0_minus_upper_b,gsl::InterpolationMethod::cubic);
    spline_omnes_0_minus_lower_a = cauchy::Interpolate(grid.interval_y,omnes_0_minus_lower_a,gsl::InterpolationMethod::cubic);
    spline_omnes_0_minus_lower_b = cauchy::Interpolate(grid.interval_y,omnes_0_minus_lower_b,gsl::InterpolationMethod::cubic);
    // C-violating
    spline_omnes_1_minus_upper_a = cauchy::Interpolate(grid.interval_y,omnes_1_minus_upper_a,gsl::InterpolationMethod::cubic);
    spline_omnes_1_minus_upper_b = cauchy::Interpolate(grid.interval_y,omnes_1_minus_upper_b,gsl::InterpolationMethod::cubic);
    spline_omnes_1_minus_lower_a = cauchy::Interpolate(grid.interval_y,omnes_1_minus_lower_a,gsl::InterpolationMethod::cubic);
    spline_omnes_1_minus_lower_b = cauchy::Interpolate(grid.interval_y,omnes_1_minus_lower_b,gsl::InterpolationMethod::cubic);
    return;
}

void OmnesSplineV3Pi::spline_omnes(){
    spline_omnes_above = cauchy::Interpolate(grid.interval_s,omnes_above,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_below = cauchy::Interpolate(grid.interval_s,omnes_below,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_egg = cauchy::Interpolate(grid.interval_phi,omnes_egg,gsl::InterpolationMethod::cubic);
    return;
}

void OmnesSplineX3Pi::spline_omnes(){
    spline_omnes_above = cauchy::Interpolate(grid.interval_s,omnes_above,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_below = cauchy::Interpolate(grid.interval_s,omnes_below,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_egg = cauchy::Interpolate(grid.interval_phi,omnes_egg,gsl::InterpolationMethod::cubic);
    return;
}

void OmnesSplineT3Pi::spline_omnes(){
    spline_omnes_above = cauchy::Interpolate(grid.interval_s,omnes_above,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_below = cauchy::Interpolate(grid.interval_s,omnes_below,gsl::InterpolationMethod::cubic);
    //
    spline_omnes_egg = cauchy::Interpolate(grid.interval_phi,omnes_egg,gsl::InterpolationMethod::cubic);
    return;
}

Complex OmnesSplineEta3Pi::operator()(double curve_parameter, int isospin, Setting evaluation)const{
    double x=curve_parameter;
    switch (evaluation) {
        case Setting::above:
            if (isospin==0) {return spline_omnes_0_above(x);}
            if (isospin==1) {return spline_omnes_1_above(x);}
            if (isospin==2) {return spline_omnes_2_above(x);}
            throw std::domain_error{"Isospin has to be either 0, 1 or 2.\n"};
        case Setting::below:
            if (isospin==0) {return spline_omnes_0_below(x);}
            if (isospin==1) {return spline_omnes_1_below(x);}
            if (isospin==2) {return spline_omnes_2_below(x);}
            throw std::domain_error{"Isospin has to be either 0, 1 or 2.\n"};
        case Setting::egg:
            if (isospin==0) {return spline_omnes_0_egg(x);}
            if (isospin==1) {return spline_omnes_1_egg(x);}
            if (isospin==2) {return spline_omnes_2_egg(x);}
            throw std::domain_error{"Isospin has to be either 0, 1 or 2.\n"};
        default:
            throw std::domain_error{"The demanded evaluation of the Omnès funtion"
                " is not needed for the chosen Process. Please check your input.\n"};
    }
}

Complex OmnesSplineEtap3Pi::operator()(double curve_parameter, int isospin, Setting evaluation)const{
    double x=curve_parameter;
    switch (evaluation) {
        case Setting::above:
            if (isospin==0) {return spline_omnes_0_above(x);}
            if (isospin==1) {return spline_omnes_1_above(x);}
            if (isospin==2) {return spline_omnes_2_above(x);}
            throw std::domain_error{"Isospin has to be either 0, 1 or 2.\n"};
        case Setting::below:
            if (isospin==0) {return spline_omnes_0_below(x);}
            if (isospin==1) {return spline_omnes_1_below(x);}
            if (isospin==2) {return spline_omnes_2_below(x);}
            throw std::domain_error{"Isospin has to be either 0, 1 or 2.\n"};
        case Setting::egg:
            if (isospin==0) {return spline_omnes_0_egg(x);}
            if (isospin==1) {return spline_omnes_1_egg(x);}
            if (isospin==2) {return spline_omnes_2_egg(x);}
            throw std::domain_error{"Isospin has to be either 0, 1 or 2.\n"};
        default:
            throw std::domain_error{"The demanded evaluation of the Omnès funtion"
                " is not needed for the chosen Process. Please check your input.\n"};
    }
}

Complex OmnesSplineEtapEtaPiPi::operator()(double curve_parameter, int isospin, Setting evaluation)const{
    double x=curve_parameter;
    switch (evaluation) {
        case Setting::above:
            if (isospin==000) {return spline_omnes_0_above(x);}
            if (isospin==111) {return spline_omnes_1_above(x);}
            if (isospin==010 || isospin==110) {return spline_omnes_eta_pi_above(x);}
            throw std::domain_error{"In eta'-> eta pi pi: wrong isospin chosen for Omnès function.\n"};
        case Setting::below:
            if (isospin==000) {return spline_omnes_0_below(x);}
            if (isospin==111) {return spline_omnes_1_below(x);}
            if (isospin==010 || isospin==110) {return spline_omnes_eta_pi_below(x);}
            throw std::domain_error{"In eta'-> eta pi pi: wrong isospin chosen for Omnès function.\n"};
        // only the pi-pi intermediate state contributes to angular average <>^-
        case Setting::minus_upper_a:
            if (isospin==000) {return spline_omnes_0_minus_upper_a(x);}
            if (isospin==111) {return spline_omnes_1_minus_upper_a(x);}
            throw std::domain_error{"In eta'-> eta pi pi the pi-pi intermediate state has to have isospin 000 or 111.\n"};
        case Setting::minus_upper_b:
            if (isospin==000) {return spline_omnes_0_minus_upper_b(x);}
            if (isospin==111) {return spline_omnes_1_minus_upper_b(x);}
            throw std::domain_error{"In eta'-> eta pi pi the pi-pi intermediate state has to have isospin 000 or 111.\n"};
        case Setting::minus_lower_a:
            if (isospin==000) {return spline_omnes_0_minus_lower_a(x);}
            if (isospin==111) {return spline_omnes_1_minus_lower_a(x);}
            throw std::domain_error{"In eta'-> eta pi pi the pi-pi intermediate state has to have isospin 000 or 111.\n"};
        case Setting::minus_lower_b:
            if (isospin==000) {return spline_omnes_0_minus_lower_b(x);}
            if (isospin==111) {return spline_omnes_1_minus_lower_b(x);}
            throw std::domain_error{"In eta'-> eta pi pi the pi-pi intermediate state has to have isospin 000 or 111.\n"};
        // only the eta-pi intermediate state contributes to angular average <>^0
        case Setting::zero_upper_a:
            if (isospin==110 || isospin==010) {return spline_omnes_eta_pi_zero_upper_a(x);}
            throw std::domain_error{"In eta'-> eta pi pi the eta-pi intermediate state has to have isospin 110.\n"};
        case Setting::zero_upper_b:
            if (isospin==110 || isospin==010) {return spline_omnes_eta_pi_zero_upper_b(x);}
            throw std::domain_error{"In eta'-> eta pi pi the eta-pi intermediate state has to have isospin 110.\n"};
        case Setting::zero_lower_a:
            if (isospin==110 || isospin==010) {return spline_omnes_eta_pi_zero_lower_a(x);}
            throw std::domain_error{"In eta'-> eta pi pi the eta-pi intermediate state has to have isospin 110.\n"};
        case Setting::zero_lower_b:
            if (isospin==110 || isospin==010) {return spline_omnes_eta_pi_zero_lower_b(x);}
            throw std::domain_error{"In eta'-> eta pi pi the eta-pi intermediate state has to have isospin 110.\n"};
        // only the eta-pi intermediate state contributes to angular average <>^+
        case Setting::plus_upper_a:
            if (isospin==110 || isospin==010) {return spline_omnes_eta_pi_plus_upper_a(x);}
            throw std::domain_error{"In eta'-> eta pi pi the eta-pi intermediate state has to have isospin 110.\n"};
        case Setting::plus_upper_b:
            if (isospin==110 || isospin==010) {return spline_omnes_eta_pi_plus_upper_b(x);}
            throw std::domain_error{"In eta'-> eta pi pi the eta-pi intermediate state has to have isospin 110.\n"};
        case Setting::plus_lower_a:
            if (isospin==110 || isospin==010) {return spline_omnes_eta_pi_plus_lower_a(x);}
            throw std::domain_error{"In eta'-> eta pi pi the eta-pi intermediate state has to have isospin 110.\n"};
        case Setting::plus_lower_b:
            if (isospin==110 || isospin==010) {return spline_omnes_eta_pi_plus_lower_b(x);}
            throw std::domain_error{"In eta'-> eta pi pi the eta-pi intermediate state has to have isospin 110.\n"};
        default:
            throw std::domain_error{"The demanded evaluation of the Omnès funtion"
                " is not needed for the chosen Process. Please check your input.\n"};
    }
}

Complex OmnesSplineV3Pi::operator()(double curve_parameter, int isospin, Setting evaluation)const{
    double x=curve_parameter;
    switch (evaluation) {
        case Setting::above:
            if (isospin==1) {return spline_omnes_above(x);}
            throw std::domain_error{"Isospin has to be 1.\n"};
        case Setting::below:
            if (isospin==1) {return spline_omnes_below(x);}
            throw std::domain_error{"Isospin has to be 1.\n"};
        case Setting::egg:
            if (isospin==1) {return spline_omnes_egg(x);}
            throw std::domain_error{"Isospin has to be 1.\n"};
        default:
            throw std::domain_error{"The demanded evaluation of the Omnès funtion"
                " is not needed for the chosen Process. Please check your input.\n"};
    }
}

Complex OmnesSplineX3Pi::operator()(double curve_parameter, int isospin, Setting evaluation)const{
    double x=curve_parameter;
    switch (evaluation) {
        case Setting::above:
            if (isospin==1) {return spline_omnes_above(x);}
            throw std::domain_error{"Isospin has to be 1.\n"};
        case Setting::below:
            if (isospin==1) {return spline_omnes_below(x);}
            throw std::domain_error{"Isospin has to be 1.\n"};
        case Setting::egg:
            if (isospin==1) {return spline_omnes_egg(x);}
            throw std::domain_error{"Isospin has to be 1.\n"};
        default:
            throw std::domain_error{"The demanded evaluation of the Omnès funtion"
                " is not needed for the chosen Process. Please check your input.\n"};
    }
}

Complex OmnesSplineT3Pi::operator()(double curve_parameter, int isospin, Setting evaluation)const{
    double x=curve_parameter;
    switch (evaluation) {
        case Setting::above:
            if (isospin==1) {return spline_omnes_above(x);}
            throw std::domain_error{"Isospin has to be 1.\n"};
        case Setting::below:
            if (isospin==1) {return spline_omnes_below(x);}
            throw std::domain_error{"Isospin has to be 1.\n"};
        case Setting::egg:
            if (isospin==1) {return spline_omnes_egg(x);}
            throw std::domain_error{"Isospin has to be 1.\n"};
        default:
            throw std::domain_error{"The demanded evaluation of the Omnès funtion"
                " is not needed for the chosen Process. Please check your input.\n"};
    }
}