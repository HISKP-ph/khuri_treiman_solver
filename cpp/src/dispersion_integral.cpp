#include "dispersion_integral.h"
#include "constants.h"
#include "array.h"
#include "matching.h"
#include <fstream>



namespace disp {


//-------------------------------------------------------------------------
//-- Constructors ---------------------------------------------------------
//-------------------------------------------------------------------------

Numerator::Numerator(double mass_1, double mass_2, double mass_3, double mass_4,
                                       FunctionSetOmnes omnes,
                                       FunctionSet matched_tilde,
                                       int isospin, double epsilon, double validity, double cutoff,
                                       SubtractionConstant sub, int max_subs)
:
s_I{pow(mass_1+mass_2,2)},s_III{pow(mass_3-mass_4,2)},
s_0{(pow(mass_1,2)+pow(mass_2,2)+pow(mass_3,2)+pow(mass_4,2))/3.},
t_I{std::pow(mass_3+mass_1,2)}, t_III{std::pow(mass_4-mass_2,2)},
t_0{(pow(mass_1,2)+pow(mass_2,2)+pow(mass_3,2)+pow(mass_4,2))/3.},
isospin{isospin}, max_subs{max_subs},
matched_tilde{matched_tilde},
omnes{omnes},
sub{sub},
matched(mass_1, mass_2, mass_3, mass_4, epsilon, validity),
matched_t_channel(mass_1,mass_3,mass_2,mass_4, epsilon, validity)
{}


NumeratorEta3Pi::NumeratorEta3Pi(double mass_1, double mass_4,
                                       FunctionSetOmnes omnes,
                                       FunctionSet matched_tilde, phase::PhaseEta3Pi phases,
                                       int isospin, double epsilon, double validity, double cutoff,
                                       SubtractionConstant sub, int max_subs)
:
Numerator(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,isospin,epsilon,validity,cutoff,sub,max_subs),
phases{phases},
mandelstam_list(mass_1, mass_4, cutoff)
{build_numerator_real_imag_parts();build_splines();build_integrand_expansions();build_numerator_deriv();}


NumeratorEtap3Pi::NumeratorEtap3Pi(double mass_1, double mass_4,
                                       FunctionSetOmnes omnes,
                                       FunctionSet matched_tilde, phase::PhaseEtap3Pi phases,
                                       int isospin, double epsilon, double validity, double cutoff,
                                       SubtractionConstant sub, int max_subs)
:
Numerator(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,isospin,epsilon,validity,cutoff,sub,max_subs),
phases{phases},
mandelstam_list(mass_1, mass_4, cutoff)
{build_numerator_real_imag_parts();build_splines();build_integrand_expansions();build_numerator_deriv();}


NumeratorEtapEtaPiPi::NumeratorEtapEtaPiPi(double mass_1, double mass_3, double mass_4,
                                       FunctionSetOmnes omnes,
                                       FunctionSet matched_tilde, phase::PhaseEtapEtaPiPi phases,
                                       int isospin, double epsilon, double validity, double cutoff,
                                       SubtractionConstant sub, int max_subs)
:
Numerator(mass_1,mass_1,mass_3,mass_4,omnes,matched_tilde,isospin,epsilon,validity,cutoff,sub,max_subs),
phases{phases},
mandelstam_list(mass_1, mass_3, mass_4, cutoff)
{build_numerator_real_imag_parts();build_splines();build_integrand_expansions();build_numerator_deriv();}

NumeratorV3Pi::NumeratorV3Pi(double mass_1, double mass_4,
                                       FunctionSetOmnes omnes,
                                       FunctionSet matched_tilde, phase::PhaseV3Pi phases,
                                       int isospin, double epsilon, double validity, double cutoff,
                                       SubtractionConstant sub, int max_subs)
:
Numerator(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,isospin,epsilon,validity,cutoff,sub,max_subs),
phases{phases},
mandelstam_list(mass_1, mass_4, cutoff)
{build_numerator_real_imag_parts();build_splines();build_integrand_expansions();build_numerator_deriv();}

NumeratorX3Pi::NumeratorX3Pi(double mass_1, double mass_4,
                                       FunctionSetOmnes omnes,
                                       FunctionSet matched_tilde, phase::PhaseX3Pi phases,
                                       int isospin, double epsilon, double validity, double cutoff,
                                       SubtractionConstant sub, int max_subs)
:
Numerator(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,isospin,epsilon,validity,cutoff,sub,max_subs),
phases{phases},
mandelstam_list(mass_1, mass_4, cutoff)
{build_numerator_real_imag_parts();build_splines();build_integrand_expansions();build_numerator_deriv();}

NumeratorT3Pi::NumeratorT3Pi(double mass_1, double mass_4,
                                       FunctionSetOmnes omnes,
                                       FunctionSet matched_tilde, phase::PhaseT3Pi phases,
                                       int isospin, double epsilon, double validity, double cutoff,
                                       SubtractionConstant sub, int max_subs)
:
Numerator(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,isospin,epsilon,validity,cutoff,sub,max_subs),
phases{phases},
mandelstam_list(mass_1, mass_4, cutoff)
{build_numerator_real_imag_parts();build_splines();build_integrand_expansions();build_numerator_deriv();}


Dispersion::Dispersion(double mass_1, double mass_2, double mass_3, double mass_4,
                                 FunctionSetOmnes omnes,
                                 FunctionSet matched_tilde,
                                 int isospin, double epsilon, double validity, double cutoff,
                                 SubtractionConstant sub, int max_subs)
:
mass_1{mass_1},
isospin{isospin},
s_I{std::pow(mass_1+mass_2,2)},s_III{std::pow(mass_4-mass_3,2)},
s_0{(std::pow(mass_1,2)+std::pow(mass_2,2)+std::pow(mass_3,2)+std::pow(mass_4,2))/3.},
t_I{std::pow(mass_3+mass_1,2)}, t_III{std::pow(mass_4-mass_2,2)},
t_0{(pow(mass_1,2)+pow(mass_2,2)+pow(mass_3,2)+pow(mass_4,2))/3.},
validity{validity}, cutoff{cutoff}, max_subs{max_subs}
{}

DispersionEta3Pi::DispersionEta3Pi(double mass_1, double mass_4,
                                 FunctionSetOmnes omnes,
                                 FunctionSet matched_tilde,
                                 phase::PhaseEta3Pi phases,
                                 int isospin, double epsilon, double validity, double cutoff,
                                 SubtractionConstant sub, int max_subs)
:
Dispersion(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,isospin,epsilon,validity,cutoff,sub,max_subs),
phases{phases},
polar_egg(mass_4,mass_1),
mandelstam_list(mass_1, mass_4, cutoff),
numeratorc(mass_1,mass_4,omnes,matched_tilde,phases,
          isospin, epsilon, validity, cutoff,sub,max_subs),
numerator_1_derivative_s_above{numeratorc.numerator_1_derivative_s_above},
numerator_1_derivative_s_below{numeratorc.numerator_1_derivative_s_below},
integrand_expansion_s{mandelstam_list.interval_disp_expansion,
    numeratorc.integrand_expansion_list_s,gsl::InterpolationMethod::cubic}
{build_dispersion_integral_lists();build_splines();}

DispersionEtap3Pi::DispersionEtap3Pi(double mass_1, double mass_4,
                                 FunctionSetOmnes omnes,
                                 FunctionSet matched_tilde,
                                 phase::PhaseEtap3Pi phases,
                                 int isospin, double epsilon, double validity, double cutoff,
                                 SubtractionConstant sub, int max_subs)
:
Dispersion(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,isospin,epsilon,validity,cutoff,sub,max_subs),
phases{phases},
polar_egg(mass_4,mass_1),
mandelstam_list(mass_1, mass_4, cutoff),
numeratorc(mass_1,mass_4,omnes,matched_tilde,phases,
          isospin, epsilon, validity, cutoff,sub,max_subs),
numerator_1_derivative_s_above{numeratorc.numerator_1_derivative_s_above},
numerator_1_derivative_s_below{numeratorc.numerator_1_derivative_s_below},
integrand_expansion_s{mandelstam_list.interval_disp_expansion,
    numeratorc.integrand_expansion_list_s,gsl::InterpolationMethod::cubic}
{build_dispersion_integral_lists();build_splines();}

DispersionEtapEtaPiPi::DispersionEtapEtaPiPi(double mass_1, double mass_3, double mass_4,
                                 FunctionSetOmnes omnes,
                                 FunctionSet matched_tilde,
                                 phase::PhaseEtapEtaPiPi phases,
                                 int isospin, double epsilon, double validity, double cutoff,
                                 SubtractionConstant sub, int max_subs)
:
Dispersion(mass_1,mass_1,mass_3,mass_4,omnes,matched_tilde,isospin,epsilon,validity,cutoff,sub,max_subs),
phases{phases},
gamma(mass_1,mass_1,mass_3,mass_4),
mandelstam_list(mass_1, mass_3, mass_4, cutoff),
numeratorc(mass_1,mass_3,mass_4,omnes,matched_tilde,phases,
          isospin, epsilon, validity, cutoff,sub,max_subs),
numerator_1_derivative_s_above{numeratorc.numerator_1_derivative_s_above},
numerator_1_derivative_s_below{numeratorc.numerator_1_derivative_s_below},
integrand_expansion_s{mandelstam_list.interval_disp_expansion,
    numeratorc.integrand_expansion_list_s,gsl::InterpolationMethod::cubic},
integrand_expansion_t{mandelstam_list.interval_t_disp_expansion,
    numeratorc.integrand_expansion_list_t,gsl::InterpolationMethod::cubic}
{build_dispersion_integral_lists();build_splines();}

DispersionV3Pi::DispersionV3Pi(double mass_1, double mass_4,
                                 FunctionSetOmnes omnes,
                                 FunctionSet matched_tilde,
                                 phase::PhaseV3Pi phases,
                                 int isospin, double epsilon, double validity, double cutoff,
                                 SubtractionConstant sub, int max_subs)
:
Dispersion(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,isospin,epsilon,validity,cutoff,sub,max_subs),
phases{phases},
polar_egg(mass_4,mass_1),
mandelstam_list(mass_1, mass_4, cutoff),
numeratorc(mass_1,mass_4,omnes,matched_tilde,phases,
          isospin, epsilon, validity, cutoff,sub,max_subs),
numerator_1_derivative_s_above{numeratorc.numerator_1_derivative_s_above},
numerator_1_derivative_s_below{numeratorc.numerator_1_derivative_s_below},
integrand_expansion_s{mandelstam_list.interval_disp_expansion,
    numeratorc.integrand_expansion_list_s,gsl::InterpolationMethod::cubic}
{build_dispersion_integral_lists();build_splines();}

DispersionX3Pi::DispersionX3Pi(double mass_1, double mass_4,
                                 FunctionSetOmnes omnes,
                                 FunctionSet matched_tilde,
                                 phase::PhaseX3Pi phases,
                                 int isospin, double epsilon, double validity, double cutoff,
                                 SubtractionConstant sub, int max_subs)
:
Dispersion(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,isospin,epsilon,validity,cutoff,sub,max_subs),
phases{phases},
polar_egg(mass_4,mass_1),
mandelstam_list(mass_1, mass_4, cutoff),
numeratorc(mass_1,mass_4,omnes,matched_tilde,phases,
          isospin, epsilon, validity, cutoff,sub,max_subs),
numerator_1_derivative_s_above{numeratorc.numerator_1_derivative_s_above},
numerator_1_derivative_s_below{numeratorc.numerator_1_derivative_s_below},
integrand_expansion_s{mandelstam_list.interval_disp_expansion,
    numeratorc.integrand_expansion_list_s,gsl::InterpolationMethod::cubic}
{build_dispersion_integral_lists();build_splines();}

DispersionT3Pi::DispersionT3Pi(double mass_1, double mass_4,
                                 FunctionSetOmnes omnes,
                                 FunctionSet matched_tilde,
                                 phase::PhaseT3Pi phases,
                                 int isospin, double epsilon, double validity, double cutoff,
                                 SubtractionConstant sub, int max_subs)
:
Dispersion(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,isospin,epsilon,validity,cutoff,sub,max_subs),
phases{phases},
polar_egg(mass_4,mass_1),
mandelstam_list(mass_1, mass_4, cutoff),
numeratorc(mass_1,mass_4,omnes,matched_tilde,phases,
          isospin, epsilon, validity, cutoff,sub,max_subs),
numerator_1_derivative_s_above{numeratorc.numerator_1_derivative_s_above},
numerator_1_derivative_s_below{numeratorc.numerator_1_derivative_s_below},
numerator_2_derivative_s_above{numeratorc.numerator_2_derivative_s_above},
numerator_2_derivative_s_below{numeratorc.numerator_2_derivative_s_below},
integrand_expansion_s{mandelstam_list.interval_disp_expansion,
    numeratorc.integrand_expansion_list_s,gsl::InterpolationMethod::cubic}
{build_dispersion_integral_lists();build_splines();}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------




//-------------------------------------------------------------------------
//-- Helpful functions ----------------------------------------------------
//-------------------------------------------------------------------------

// as we use the absolute value of the omnes function, it does not matter if we choose 'above' or 'below'
Complex NumeratorEta3Pi::numerator(double s)const{
    // ADJUST FOR DIFFERENT SUBTRAKTIONSCHEME
    // the power of 's' in the denominator determines the convergence
    if (isospin==0) {
        return std::sin(phases.phase(s,isospin))/(std::pow(s,2.)*std::abs(omnes(s,isospin, Setting::above)))
        *matched_tilde(s, isospin, sub);
    }
    // For the chosen subtraction-scheme the latter two cases are the same for the C-conserving and C-violating Process
    if (isospin==1) {
        return std::sin(phases.phase(s,isospin))/(std::pow(s,1.)*std::abs(omnes(s,isospin, Setting::above)))
        * matched_tilde(s, isospin, sub);
    }
    if (isospin==2) {
        return std::sin(phases.phase(s,isospin))/(std::pow(s,1.)*std::abs(omnes(s,isospin, Setting::above)))
        * matched_tilde(s, isospin, sub);
    }
    throw std::domain_error{"Isospin not allowed for the chosen decay, it should be either 0, 1 or 2."};
}

Complex NumeratorEtap3Pi::numerator(double s)const{
    // ADJUST FOR DIFFERENT SUBTRAKTIONSCHEME
    // the power of 's' in the denominator determines the convergence
    if (isospin==0) {
        return std::sin(phases.phase(s,isospin))/(std::pow(s,2.)*std::abs(omnes(s,isospin, Setting::above)))
        *matched_tilde(s, isospin, sub);
    }
    // For the chosen subtraction-scheme the latter two cases are the same for the C-conserving and C-violating Process
    if (isospin==1) {
        return std::sin(phases.phase(s,isospin))/(std::pow(s,1.)*std::abs(omnes(s,isospin, Setting::above)))
        * matched_tilde(s, isospin, sub);
    }
    if (isospin==2) {
        return std::sin(phases.phase(s,isospin))/(std::pow(s,1.)*std::abs(omnes(s,isospin, Setting::above)))
        * matched_tilde(s, isospin, sub);
    }
    throw std::domain_error{"Isospin not allowed for the chosen decay, it should be either 0, 1 or 2."};
}


Complex NumeratorEtapEtaPiPi::numerator(double s)const{
    // ADJUST FOR DIFFERENT SUBTRAKTIONSCHEME
    // the power of 's' in the denominator determines the convergence
    if (isospin==000) {
        return std::sin(phases.phase(s,0))/(std::pow(s,3.)*std::abs(omnes(s,isospin,Setting::above)))
        *matched_tilde(s, isospin, sub);
    }
    if (isospin==010) {
        return std::sin(phases.phase_eta_pi(s))/(std::pow(s,3.)*std::abs(omnes(s,isospin,Setting::above)))
        * matched_tilde(s, isospin, sub);
    }
    if (isospin==111) {
        return std::sin(phases.phase(s,1))/(std::pow(s,1.)*std::abs(omnes(s,isospin,Setting::above)))
        *matched_tilde(s, isospin, sub);
    }
    if (isospin==110) {
        return std::sin(phases.phase_eta_pi(s))/(std::pow(s,2.)*std::abs(omnes(s,isospin,Setting::above)))
        * matched_tilde(s, isospin, sub);
    }
    throw std::domain_error{"Isospin not allowed for the chosen decay, should be either 000, 111, 010 or 110."};
}

Complex NumeratorV3Pi::numerator(double s)const{
    if (isospin==1) {
        return std::sin(phases.phase(s,isospin))/(std::pow(s,max_subs)*std::abs(omnes(s,isospin, Setting::above)))
        * matched_tilde(s, isospin, sub);
    }
    throw std::domain_error{"Isospin not allowed for the chosen decay, it should be 1."};
}

Complex NumeratorX3Pi::numerator(double s)const{
    if (isospin==1) {
        return std::sin(phases.phase(s,isospin))/(std::pow(s,max_subs)*std::abs(omnes(s,isospin, Setting::above)))
        * matched_tilde(s, isospin, sub);
    }
    throw std::domain_error{"Isospin not allowed for the chosen decay, it should be 1."};
}

Complex NumeratorT3Pi::numerator(double s)const{
    if (isospin==1) {
        return std::sin(phases.phase(s,isospin))/(std::pow(s,max_subs)*std::abs(omnes(s,isospin, Setting::above)))
        * matched_tilde(s, isospin, sub);
    }
    throw std::domain_error{"Isospin not allowed for the chosen decay, it should be 1."};
}


void NumeratorEta3Pi::build_numerator_real_imag_parts(){
    for (std::size_t i=0; i<mandelstam_list.interval_disp_expansion.size(); i++) {
        real_numerator_s.push_back(real(numerator(mandelstam_list.interval_disp_expansion[i])));
        imag_numerator_s.push_back(imag(numerator(mandelstam_list.interval_disp_expansion[i])));
    }
    return;
}

void NumeratorEtap3Pi::build_numerator_real_imag_parts(){
    for (std::size_t i=0; i<mandelstam_list.interval_disp_expansion.size(); i++) {
        real_numerator_s.push_back(real(numerator(mandelstam_list.interval_disp_expansion[i])));
        imag_numerator_s.push_back(imag(numerator(mandelstam_list.interval_disp_expansion[i])));
    }
    return;
}

void NumeratorEtapEtaPiPi::build_numerator_real_imag_parts(){
    for (std::size_t i=0; i<mandelstam_list.interval_disp_expansion.size(); i++) {
        real_numerator_s.push_back(real(numerator(mandelstam_list.interval_disp_expansion[i])));
        imag_numerator_s.push_back(imag(numerator(mandelstam_list.interval_disp_expansion[i])));
    }
    for (std::size_t i=0; i<mandelstam_list.interval_t_disp_expansion.size(); i++) {
        real_numerator_t.push_back(real(numerator(mandelstam_list.interval_t_disp_expansion[i])));
        imag_numerator_t.push_back(imag(numerator(mandelstam_list.interval_t_disp_expansion[i])));
    }
    return;
}

void NumeratorV3Pi::build_numerator_real_imag_parts(){
    for (std::size_t i=0; i<mandelstam_list.interval_disp_expansion.size(); i++) {
        real_numerator_s.push_back(real(numerator(mandelstam_list.interval_disp_expansion[i])));
        imag_numerator_s.push_back(imag(numerator(mandelstam_list.interval_disp_expansion[i])));
    }
    return;
}

void NumeratorX3Pi::build_numerator_real_imag_parts(){
    for (std::size_t i=0; i<mandelstam_list.interval_disp_expansion.size(); i++) {
        real_numerator_s.push_back(real(numerator(mandelstam_list.interval_disp_expansion[i])));
        imag_numerator_s.push_back(imag(numerator(mandelstam_list.interval_disp_expansion[i])));
    }
    return;
}

void NumeratorT3Pi::build_numerator_real_imag_parts(){
    for (std::size_t i=0; i<mandelstam_list.interval_disp_expansion.size(); i++) {
        real_numerator_s.push_back(real(numerator(mandelstam_list.interval_disp_expansion[i])));
        imag_numerator_s.push_back(imag(numerator(mandelstam_list.interval_disp_expansion[i])));
    }
    return;
}


void NumeratorEta3Pi::build_splines(){
    real_num_s = gsl::Interpolate(mandelstam_list.interval_disp_expansion,real_numerator_s,gsl::InterpolationMethod::cubic);
    imag_num_s = gsl::Interpolate(mandelstam_list.interval_disp_expansion,imag_numerator_s,gsl::InterpolationMethod::cubic);
    return;
}

void NumeratorEtap3Pi::build_splines(){
    real_num_s = gsl::Interpolate(mandelstam_list.interval_disp_expansion,real_numerator_s,gsl::InterpolationMethod::cubic);
    imag_num_s = gsl::Interpolate(mandelstam_list.interval_disp_expansion,imag_numerator_s,gsl::InterpolationMethod::cubic);
    return;
}

void NumeratorEtapEtaPiPi::build_splines(){
    real_num_s = gsl::Interpolate(mandelstam_list.interval_disp_expansion,real_numerator_s,gsl::InterpolationMethod::cubic);
    imag_num_s = gsl::Interpolate(mandelstam_list.interval_disp_expansion,imag_numerator_s,gsl::InterpolationMethod::cubic);
    real_num_t = gsl::Interpolate(mandelstam_list.interval_t_disp_expansion,real_numerator_t,gsl::InterpolationMethod::cubic);
    imag_num_t = gsl::Interpolate(mandelstam_list.interval_t_disp_expansion,imag_numerator_t,gsl::InterpolationMethod::cubic);
    return;
}

void NumeratorV3Pi::build_splines(){
    real_num_s = gsl::Interpolate(mandelstam_list.interval_disp_expansion,real_numerator_s,gsl::InterpolationMethod::cubic);
    imag_num_s = gsl::Interpolate(mandelstam_list.interval_disp_expansion,imag_numerator_s,gsl::InterpolationMethod::cubic);
    return;
}

void NumeratorX3Pi::build_splines(){
    real_num_s = gsl::Interpolate(mandelstam_list.interval_disp_expansion,real_numerator_s,gsl::InterpolationMethod::cubic);
    imag_num_s = gsl::Interpolate(mandelstam_list.interval_disp_expansion,imag_numerator_s,gsl::InterpolationMethod::cubic);
    return;
}

void NumeratorT3Pi::build_splines(){
    real_num_s = gsl::Interpolate(mandelstam_list.interval_disp_expansion,real_numerator_s,gsl::InterpolationMethod::cubic);
    imag_num_s = gsl::Interpolate(mandelstam_list.interval_disp_expansion,imag_numerator_s,gsl::InterpolationMethod::cubic);
    return;
}


void Numerator::build_numerator_deriv(){
    // the value +0.1 is arbitrary. The only important fact is that one evaluates slightly above and below s_III
    numerator_1_derivative_s_above = matched.p_wave_parameter(s_III+0.1, real_num_s) + 1.i*matched.p_wave_parameter(s_III+0.1, imag_num_s);
    numerator_1_derivative_s_below = matched.p_wave_parameter(s_III-0.1, real_num_s) + 1.i*matched.p_wave_parameter(s_III-0.1, imag_num_s);
}

void NumeratorT3Pi::build_numerator_deriv(){
    // the value +0.1 is arbitrary. The only important fact is that one evaluates slightly above and below s_III
    auto [numerator_1_derivative_s_above_real,numerator_2_derivative_s_above_real] = matched.d_wave_parameters(s_III+0.1, real_num_s);
    auto [numerator_1_derivative_s_above_imag,numerator_2_derivative_s_above_imag] = matched.d_wave_parameters(s_III+0.1, imag_num_s);
    auto [numerator_1_derivative_s_below_real,numerator_2_derivative_s_below_real] = matched.d_wave_parameters(s_III-0.1, real_num_s);
    auto [numerator_1_derivative_s_below_imag,numerator_2_derivative_s_below_imag] = matched.d_wave_parameters(s_III-0.1, imag_num_s);
    numerator_1_derivative_s_above = numerator_1_derivative_s_above_real + 1.i* numerator_1_derivative_s_above_imag;
    numerator_2_derivative_s_above = numerator_2_derivative_s_above_real + 1.i* numerator_2_derivative_s_above_imag;
    numerator_1_derivative_s_below = numerator_1_derivative_s_below_real + 1.i* numerator_1_derivative_s_below_imag;
    numerator_2_derivative_s_below = numerator_2_derivative_s_below_real + 1.i* numerator_2_derivative_s_below_imag;
}



void NumeratorEta3Pi::build_integrand_expansions(){
    
    // Fill here the respective lists with pi-pi intermediate states
    for (std::size_t j=0; j<mandelstam_list.interval_disp_expansion.size(); ++j) {
        if (isospin==0 || isospin==2) {
            integrand_expansion_list_s.push_back(matched.s_wave_integrand_expansion
            (mandelstam_list.interval_disp_expansion[j], real_num_s)
            +1.i*matched.s_wave_integrand_expansion
            (mandelstam_list.interval_disp_expansion[j], imag_num_s));
        }
        // match a P-wave
        if (isospin==1) {
            integrand_expansion_list_s.push_back
            (matched.p_wave_integrand_expansion(mandelstam_list.interval_disp_expansion[j], real_num_s)
            +1.i*matched.p_wave_integrand_expansion
            (mandelstam_list.interval_disp_expansion[j], imag_num_s));
        }
    }
    return;
}

void NumeratorEtap3Pi::build_integrand_expansions(){
    
    // Fill here the respective lists with pi-pi intermediate states
    for (std::size_t j=0; j<mandelstam_list.interval_disp_expansion.size(); ++j) {
        if (isospin==0 || isospin==2) {
            integrand_expansion_list_s.push_back(matched.s_wave_integrand_expansion
            (mandelstam_list.interval_disp_expansion[j], real_num_s)
            +1.i*matched.s_wave_integrand_expansion
            (mandelstam_list.interval_disp_expansion[j], imag_num_s));
        }
        // match a P-wave
        if (isospin==1) {
            integrand_expansion_list_s.push_back
            (matched.p_wave_integrand_expansion(mandelstam_list.interval_disp_expansion[j], real_num_s)
            +1.i*matched.p_wave_integrand_expansion
            (mandelstam_list.interval_disp_expansion[j], imag_num_s));
        }
    }
    return;
}

void NumeratorEtapEtaPiPi::build_integrand_expansions(){
    
    // Fill here the respective lists with pi-pi intermediate states
    for (std::size_t j=0; j<mandelstam_list.interval_disp_expansion.size(); ++j) {
        if (isospin==000) {
            integrand_expansion_list_s.push_back
            (matched.s_wave_integrand_expansion(mandelstam_list.interval_disp_expansion[j], real_num_s)
            +1.i*matched.s_wave_integrand_expansion(mandelstam_list.interval_disp_expansion[j], imag_num_s));
        }
        if (isospin==111) {
            integrand_expansion_list_s.push_back
            (matched.p_wave_integrand_expansion(mandelstam_list.interval_disp_expansion[j], real_num_s)
            +1.i*matched.p_wave_integrand_expansion(mandelstam_list.interval_disp_expansion[j], imag_num_s));
        }
        if (isospin==010 || isospin==110) {
            integrand_expansion_list_s.push_back(42.);
                // 'integrand_expansion_list_s' is not needed for isospin 010 and 110 (i.e. eta-pi inetermediate states),
                // therefore fill the list with some arbitrary number,
                // otherwise the interpolation within the constructor will complain
        }
    }
    
    // Fill here the respective lists with eta-pi intermediate states
    for (std::size_t j=0; j<mandelstam_list.interval_t_disp_expansion.size(); ++j) {
        if (isospin==010 || isospin==110) {
            integrand_expansion_list_t.push_back
            (matched_t_channel.s_wave_integrand_expansion(mandelstam_list.interval_t_disp_expansion[j], real_num_t)
            +1.i*matched_t_channel.s_wave_integrand_expansion(mandelstam_list.interval_t_disp_expansion[j], imag_num_t));
        }
        if (isospin==000 || isospin==111) {
            integrand_expansion_list_t.push_back(42.);
                // 'integrand_expansion_list_t' is not needed for isospin 000 and 111 (i.e. pi-pi inetermediate states),
                // therefore fill the list with some arbitrary number,
                // otherwise the interpolation within the constructor will complain
        }
    }
    return;
}

void NumeratorV3Pi::build_integrand_expansions(){
    
    // Fill here the respective lists with pi-pi intermediate states
    for (std::size_t j=0; j<mandelstam_list.interval_disp_expansion.size(); ++j) {
        // match a P-wave
        if (isospin==1) {
            integrand_expansion_list_s.push_back
            (matched.p_wave_integrand_expansion(mandelstam_list.interval_disp_expansion[j], real_num_s)
            +1.i*matched.p_wave_integrand_expansion
            (mandelstam_list.interval_disp_expansion[j], imag_num_s));
        }
        else{
            std::domain_error{"Isospin not allowed for the chosen decay, it should be 1."};
        }
    }
    return;
}

void NumeratorX3Pi::build_integrand_expansions(){
    
    // Fill here the respective lists with pi-pi intermediate states
    for (std::size_t j=0; j<mandelstam_list.interval_disp_expansion.size(); ++j) {
        // match a P-wave
        if (isospin==1) {
            integrand_expansion_list_s.push_back
            (matched.p_wave_integrand_expansion(mandelstam_list.interval_disp_expansion[j], real_num_s)
            +1.i*matched.p_wave_integrand_expansion
            (mandelstam_list.interval_disp_expansion[j], imag_num_s));
        }
        else{
            std::domain_error{"Isospin not allowed for the chosen decay, it should be 1."};
        }
    }
    return;
}

void NumeratorT3Pi::build_integrand_expansions(){
    
    // Fill here the respective lists with pi-pi intermediate states
    for (std::size_t j=0; j<mandelstam_list.interval_disp_expansion.size(); ++j) {
        // match a D-wave
        if (isospin==1) {
            integrand_expansion_list_s.push_back
            (matched.d_wave_integrand_expansion(mandelstam_list.interval_disp_expansion[j], real_num_s)
            +1.i*matched.d_wave_integrand_expansion
            (mandelstam_list.interval_disp_expansion[j], imag_num_s));
        }
        else{
            std::domain_error{"Isospin not allowed for the chosen decay, it should be 1."};
        }
    }
    return;
}
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------










//------------------------------------------------------------------------------------
//-- Numeric parts of dispersion integrals -------------------------------------------
//------------------------------------------------------------------------------------

Complex DispersionEta3Pi::numerator(double s)const{
    return numeratorc.numerator(s);
}

Complex DispersionEtap3Pi::numerator(double s)const{
    return numeratorc.numerator(s);
}

Complex DispersionEtapEtaPiPi::numerator(double s)const{
    return numeratorc.numerator(s);
}

Complex DispersionV3Pi::numerator(double s)const{
    return numeratorc.numerator(s);
}

Complex DispersionX3Pi::numerator(double s)const{
    return numeratorc.numerator(s);
}

Complex DispersionT3Pi::numerator(double s)const{
    return numeratorc.numerator(s);
}
// Integrands with Cauchy singularity

Complex DispersionEta3Pi::integrand_cauchy(double s, double s_prime)const{
    // for S-waves
    if (isospin==0 || isospin==2) {
        return (numerator(s_prime)-numerator(s)
                )/(std::sqrt(Complex(s_III-s_prime))*(s_prime-s));
    }
    // for P-waves
    if (isospin==1) {
        return (numerator(s_prime)-numerator(s)
                )/(std::pow(std::sqrt(Complex(s_III-s_prime)),3.0)*(s_prime-s));
    }
    throw std::domain_error{"Wrong isospin inserted."};
}

Complex DispersionEtap3Pi::integrand_cauchy(double s, double s_prime)const{
    // for S-waves
    if (isospin==0 || isospin==2) {
        return (numerator(s_prime)-numerator(s)
                )/(std::sqrt(Complex(s_III-s_prime))*(s_prime-s));
    }
    // for P-waves
    if (isospin==1) {
        return (numerator(s_prime)-numerator(s)
                )/(std::pow(std::sqrt(Complex(s_III-s_prime)),3.0)*(s_prime-s));
    }
    throw std::domain_error{"Wrong isospin inserted."};
}


Complex DispersionEtapEtaPiPi::integrand_cauchy(double s, double s_prime)const{
    // for S-waves
    if (isospin==000) {
        return (numerator(s_prime)-numerator(s)
                )/(std::sqrt(Complex(s_III-s_prime))*(s_prime-s));
    }
    // This isospin configuration is for the t-channel in eta'-> eta pi pi with a different pseudo-threshold t_III
    // (still s-wave)
    if (isospin==010 || isospin==110) {
        double t=s;
        return (numerator(s_prime)-numerator(t)
                )/(std::sqrt(Complex(t_III-s_prime))*(s_prime-t));
    }
    // for P-waves
    if (isospin==111) {
        return (numerator(s_prime)-numerator(s)
                )/(std::pow(std::sqrt(Complex(s_III-s_prime)),3.0)*(s_prime-s));
    }
    throw std::domain_error{"Wrong isospin inserted."};
}

Complex DispersionV3Pi::integrand_cauchy(double s, double s_prime)const{
    // for P-waves
    if (isospin==1) {
        return (numerator(s_prime)-numerator(s)
                )/(std::pow(std::sqrt(Complex(s_III-s_prime)),3.0)*(s_prime-s));
    }
    throw std::domain_error{"Wrong isospin inserted."};
}

Complex DispersionX3Pi::integrand_cauchy(double s, double s_prime)const{
    // for P-waves
    if (isospin==1) {
        return (numerator(s_prime)-numerator(s)
                )/(std::pow(std::sqrt(Complex(s_III-s_prime)),3.0)*(s_prime-s));
    }
    throw std::domain_error{"Wrong isospin inserted."};
}

Complex DispersionT3Pi::integrand_cauchy(double s, double s_prime)const{
    // for P-waves
    if (isospin==1) {
        return (numerator(s_prime)-numerator(s)
                )/(std::pow(std::sqrt(Complex(s_III-s_prime)),5.0)*(s_prime-s));
    }
    throw std::domain_error{"Wrong isospin inserted."};
}


Complex DispersionEta3Pi::integrand_pseudo(Complex s, double s_prime)const{
    double x_III = s_III;
    // for S-waves
    if (isospin==0 || isospin==2) {
        // close to singularity s_III insert the well-defined expansion
        if (s_prime>=x_III-validity && s_prime<=x_III+validity) {
            return integrand_expansion_s(s_prime)/(s_prime-s); // use here the expansion around the pseudothreshold s_III
        }
        return (numerator(s_prime)-numerator(x_III)
                )/(std::sqrt(Complex(x_III-s_prime))*(s_prime-s));
    }
    // for P-waves
    if (isospin==1) {
        // close to singularity s_III insert the well-defined expansion
        if (s_prime>=x_III-validity && s_prime<=x_III+validity) {
            return integrand_expansion_s(s_prime)/(s_prime-s);
        }
        if (s_prime<x_III-validity) {
            return (
                    numerator(s_prime)
                    -numerator(x_III)
                    -(x_III-s_prime) * numerator_1_derivative_s_below
                    )/(std::pow(std::sqrt(Complex(x_III-s_prime)),3.0)*(s_prime-s));
        }
        if (s_prime>x_III+validity) {
            return (
                    numerator(s_prime)
                    -numerator(x_III)
                    -(x_III-s_prime) * numerator_1_derivative_s_above
                    )/(std::pow(std::sqrt(Complex(x_III-s_prime)),3.0)*(s_prime-s));
        }  
        else{throw std::domain_error{"Value of Mandelstam s is not allowed."};}
    }    
    else{throw std::domain_error{"Something goes wrong in integrand_pseudo for P-wave."};}
}

Complex DispersionEtap3Pi::integrand_pseudo(Complex s, double s_prime)const{
    double x_III = s_III;
    // for S-waves
    if (isospin==0 || isospin==2) {
        // close to singularity s_III insert the well-defined expansion
        if (s_prime>=x_III-validity && s_prime<=x_III+validity) {
            return integrand_expansion_s(s_prime)/(s_prime-s); // use here the expansion around the pseudothreshold s_III
        }
        return (numerator(s_prime)-numerator(x_III)
                )/(std::sqrt(Complex(x_III-s_prime))*(s_prime-s));
    }
    // for P-waves
    if (isospin==1) {
        // close to singularity s_III insert the well-defined expansion
        if (s_prime>=x_III-validity && s_prime<=x_III+validity) {
            return integrand_expansion_s(s_prime)/(s_prime-s);
        }
        if (s_prime<x_III-validity) {
            return (
                    numerator(s_prime)
                    -numerator(x_III)
                    -(x_III-s_prime) * numerator_1_derivative_s_below
                    )/(std::pow(std::sqrt(Complex(x_III-s_prime)),3.0)*(s_prime-s));
        }
        if (s_prime>x_III+validity) {
            return (
                    numerator(s_prime)
                    -numerator(x_III)
                    -(x_III-s_prime) * numerator_1_derivative_s_above
                    )/(std::pow(std::sqrt(Complex(x_III-s_prime)),3.0)*(s_prime-s));
        }  
        else{throw std::domain_error{"Value of Mandelstam s is not allowed."};}
    }    
    else{throw std::domain_error{"Something goes wrong in integrand_pseudo for P-wave."};}
}

Complex DispersionEtapEtaPiPi::integrand_pseudo(Complex s, double s_prime)const{
    double x_III;
    // Consider a different pseudo-threshold t_III that belongs to the t-channel in eta'-> eta pi pi
    // (still s-wave)
    switch (isospin) {
        case 010: case 110:   // isospin configurations for t-channel in eta'-> eta pi pi
            x_III=t_III;
            break;
            
        default:
            x_III=s_III;
            break;
    }
    // for S-waves
    if (isospin==000) {
        // close to singularity s_III insert the well-defined expansion
        if (s_prime>=x_III-validity && s_prime<=x_III+validity) {
            return integrand_expansion_s(s_prime)/(s_prime-s); // use here the expansion around the pseudothreshold s_III
        }
        return (numerator(s_prime)-numerator(x_III)
                )/(std::sqrt(Complex(x_III-s_prime))*(s_prime-s));
    }
    if (isospin==010 || isospin==110) {
        // close to singularity t_III insert the well-defined expansion
        if (s_prime>=x_III-validity && s_prime<=x_III+validity) {
            return integrand_expansion_t(s_prime)/(s_prime-s); // use here the expansion around the pseudothreshold t_III
        }
        return (numerator(s_prime)-numerator(x_III)
                )/(std::sqrt(Complex(x_III-s_prime))*(s_prime-s));
    }
    // for P-waves
    if (isospin==111) {
        // close to singularity s_III insert the well-defined expansion
        if (s_prime>=x_III-validity && s_prime<=x_III+validity) {
            return integrand_expansion_s(s_prime)/(s_prime-s);
        }
        if (s_prime<x_III-validity) {
            return (
                    numerator(s_prime)
                    -numerator(x_III)
                    -(x_III-s_prime) * numerator_1_derivative_s_below
                    )/(std::pow(std::sqrt(Complex(x_III-s_prime)),3.0)*(s_prime-s));
        }
        if (s_prime>x_III+validity) {
            return (
                    numerator(s_prime)
                    -numerator(x_III)
                    -(x_III-s_prime) * numerator_1_derivative_s_above
                    )/(std::pow(std::sqrt(Complex(x_III-s_prime)),3.0)*(s_prime-s));
        }
        
        else{throw std::domain_error{"Value of Mandelstam s is not allowed."};}
    }
    
    else{throw std::domain_error{"Something goes wrong in integrand_pseudo for P-wave."};}
}

Complex DispersionV3Pi::integrand_pseudo(Complex s, double s_prime)const{
    double x_III = s_III;
    // for P-waves
    if (isospin==1) {
        // close to singularity s_III insert the well-defined expansion
        if (s_prime>=x_III-validity && s_prime<=x_III+validity) {
            return integrand_expansion_s(s_prime)/(s_prime-s);
        }
        if (s_prime<x_III-validity) {
            return (
                    numerator(s_prime)
                    -numerator(x_III)
                    -(x_III-s_prime) * numerator_1_derivative_s_below
                    )/(std::pow(std::sqrt(Complex(x_III-s_prime)),3.0)*(s_prime-s));
        }
        if (s_prime>x_III+validity) {
            return (
                    numerator(s_prime)
                    -numerator(x_III)
                    -(x_III-s_prime) * numerator_1_derivative_s_above
                    )/(std::pow(std::sqrt(Complex(x_III-s_prime)),3.0)*(s_prime-s));
        }  
        else{throw std::domain_error{"Value of Mandelstam s is not allowed."};}
    }    
    else{throw std::domain_error{"Something goes wrong in integrand_pseudo for P-wave."};}
}

Complex DispersionX3Pi::integrand_pseudo(Complex s, double s_prime)const{
    double x_III = s_III;
    // for P-waves
    if (isospin==1) {
        // close to singularity s_III insert the well-defined expansion
        if (s_prime>=x_III-validity && s_prime<=x_III+validity) {
            return integrand_expansion_s(s_prime)/(s_prime-s);
        }
        if (s_prime<x_III-validity) {
            return (
                    numerator(s_prime)
                    -numerator(x_III)
                    -(x_III-s_prime) * numerator_1_derivative_s_below
                    )/(std::pow(std::sqrt(Complex(x_III-s_prime)),3.0)*(s_prime-s));
        }
        if (s_prime>x_III+validity) {
            return (
                    numerator(s_prime)
                    -numerator(x_III)
                    -(x_III-s_prime) * numerator_1_derivative_s_above
                    )/(std::pow(std::sqrt(Complex(x_III-s_prime)),3.0)*(s_prime-s));
        }  
        else{throw std::domain_error{"Value of Mandelstam s is not allowed."};}
    }    
    else{throw std::domain_error{"Something goes wrong in integrand_pseudo for P-wave."};}
}

Complex DispersionT3Pi::integrand_pseudo(Complex s, double s_prime)const{
    double x_III = s_III;
    // for D-waves
    if (isospin==1) {
        // close to singularity s_III insert the well-defined expansion
        if (s_prime>=x_III-validity && s_prime<=x_III+validity) {
            return integrand_expansion_s(s_prime)/(s_prime-s);
        }
        if (s_prime<x_III-validity) {
            return (
                    numerator(s_prime)
                    -numerator(x_III)
                    -(x_III-s_prime) * numerator_1_derivative_s_below
                    -std::pow(x_III-s_prime,2.) * numerator_2_derivative_s_below
                    )/(std::pow(std::sqrt(Complex(x_III-s_prime)),5.0)*(s_prime-s));
        }
        if (s_prime>x_III+validity) {
            return (
                    numerator(s_prime)
                    -numerator(x_III)
                    -(x_III-s_prime) * numerator_1_derivative_s_above
                    -std::pow(x_III-s_prime,2.) * numerator_2_derivative_s_below
                    )/(std::pow(std::sqrt(Complex(x_III-s_prime)),5.0)*(s_prime-s));
        }  
        else{throw std::domain_error{"Value of Mandelstam s is not allowed."};}
    }    
    else{throw std::domain_error{"Something goes wrong in integrand_pseudo for P-wave."};}
}





//------------------------------------------------------------------------------------
//-- Integrals -----------------------------------------------------------------------
//------------------------------------------------------------------------------------

// Adjust here the number of integrations nodes for Gaussian quadrature.
// It is recommended to use a high precision, i.e. the adaptive routine by Setting 'adaptive=true', close to critical points
// like pseudo-threshold and cusp due to I=0 phase shift

Complex DispersionEta3Pi::numeric_integral_pseudo(Complex s, double lower_limit, double upper_limit, bool adaptive)const{
    Complex mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_pseudo(mandelstam_s,s_prime); }};
    
    int nodes = 500;

    return complex_integration(integrand, lower_limit, upper_limit, adaptive, nodes);
}

Complex DispersionEtap3Pi::numeric_integral_pseudo(Complex s, double lower_limit, double upper_limit, bool adaptive)const{
    Complex mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_pseudo(mandelstam_s,s_prime); }};
    
    int nodes = 200;
    
    return complex_integration(integrand, lower_limit, upper_limit, adaptive, nodes);
}

Complex DispersionEtapEtaPiPi::numeric_integral_pseudo(Complex s, double lower_limit, double upper_limit, bool adaptive)const{
    Complex mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_pseudo(mandelstam_s,s_prime); }};
    
    int nodes = 100;
    
    return complex_integration(integrand, lower_limit, upper_limit, adaptive, nodes);
}

Complex DispersionV3Pi::numeric_integral_pseudo(Complex s, double lower_limit, double upper_limit, bool adaptive)const{
    Complex mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_pseudo(mandelstam_s,s_prime); }};
    
    int nodes = 150;

    return complex_integration(integrand, lower_limit, upper_limit, adaptive, nodes);
}

Complex DispersionX3Pi::numeric_integral_pseudo(Complex s, double lower_limit, double upper_limit, bool adaptive)const{
    Complex mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_pseudo(mandelstam_s,s_prime); }};
    
    int nodes = 150;

    return complex_integration(integrand, lower_limit, upper_limit, adaptive, nodes);
}

Complex DispersionT3Pi::numeric_integral_pseudo(Complex s, double lower_limit, double upper_limit, bool adaptive)const{
    Complex mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_pseudo(mandelstam_s,s_prime); }};
    
    int nodes = 150;

    return complex_integration(integrand, lower_limit, upper_limit, adaptive, nodes);
}
    


Complex DispersionEta3Pi::numeric_integral_cauchy(double s, double lower_limit, double upper_limit, bool adaptive)const{
    double mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_cauchy(mandelstam_s,s_prime); }};
    
    int nodes = 150;
    
    return complex_integration(integrand, lower_limit, upper_limit, adaptive, nodes);
}

Complex DispersionEtap3Pi::numeric_integral_cauchy(double s, double lower_limit, double upper_limit, bool adaptive)const{
    double mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_cauchy(mandelstam_s,s_prime); }};
    
    int nodes = 200;
    
    return complex_integration(integrand, lower_limit, upper_limit, adaptive, nodes);
}


Complex DispersionEtapEtaPiPi::numeric_integral_cauchy(double s, double lower_limit, double upper_limit, bool adaptive)const{
    double mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_cauchy(mandelstam_s,s_prime); }};
    
    int nodes = 100;
    
    return complex_integration(integrand, lower_limit, upper_limit, adaptive, nodes);
}


Complex DispersionV3Pi::numeric_integral_cauchy(double s, double lower_limit, double upper_limit, bool adaptive)const{
    double mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_cauchy(mandelstam_s,s_prime); }};
    
    int nodes = 150;
    
    return complex_integration(integrand, lower_limit, upper_limit, adaptive, nodes);
}

Complex DispersionX3Pi::numeric_integral_cauchy(double s, double lower_limit, double upper_limit, bool adaptive)const{
    double mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_cauchy(mandelstam_s,s_prime); }};
    
    int nodes = 150;
    
    return complex_integration(integrand, lower_limit, upper_limit, adaptive, nodes);
}

Complex DispersionT3Pi::numeric_integral_cauchy(double s, double lower_limit, double upper_limit, bool adaptive)const{
    double mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_cauchy(mandelstam_s,s_prime); }};
    
    int nodes = 150;
    
    return complex_integration(integrand, lower_limit, upper_limit, adaptive, nodes);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------





//------------------------------------------------------------------------------------
//-- Analytic parts of dispersion integrals ------------------------------------------
//------------------------------------------------------------------------------------

// errors in these functions have a huge effect for Mandelstam s close to
// pseudo-threshold, but may not be visible at other energies


Complex Dispersion::analytic_integral_R12(Complex s, double min, double max, Setting evaluation)const{
    double x_III;
    // t-channel in eta'-> eta pi pi has a different pseudo-threshold t_III
    switch (isospin) {
        case 010: case 110: x_III=t_III; break;
        default: x_III=s_III; break;
    }
    
    if (evaluation==Setting::above) {
        return ( std::log(
                     Complex(
                                  (
                                   std::sqrt(Complex(x_III-min))
                                   +std::sqrt(Complex(x_III-s))
                                  )/(
                                     std::sqrt(Complex(x_III-min))
                                     -std::sqrt(Complex(x_III-s))
                                     )
                                  )
                     )
                +
                std::log(
                    Complex(
                                 (
                                  std::sqrt(Complex(x_III-s))
                                  -std::sqrt(Complex(x_III-max))
                                  )/(
                                     std::sqrt(Complex(x_III-s))
                                     +std::sqrt(Complex(x_III-max))
                                     )
                                 )
                    )
                + 1.i*pi()
                )/std::sqrt(Complex(x_III-s));
    }
    if (evaluation==Setting::below) {
        return ( std::log(
                     Complex(
                                  (
                                   std::sqrt(Complex(x_III-min))
                                   +std::sqrt(Complex(x_III-s))
                                  )/(
                                     std::sqrt(Complex(x_III-min))
                                     -std::sqrt(Complex(x_III-s))
                                     )
                                  )
                     )
                +
                std::log(
                    Complex(
                                 (
                                  std::sqrt(Complex(x_III-s))
                                  -std::sqrt(Complex(x_III-max))
                                  )/(
                                     std::sqrt(Complex(x_III-s))
                                     +std::sqrt(Complex(x_III-max))
                                     )
                                 )
                    )
                - 1.i*pi()
                )/std::sqrt(Complex(x_III-s));
    }
    else throw std::domain_error{" 'evaluation' should be either 'above' or 'below'."};
}




Complex Dispersion::analytic_integral_Q12(Complex s, double min, double max)const{
    double x_III;
    // t-channel in eta'-> eta pi pi has a different pseudo-threshold t_III
    switch (isospin) {
        case 010: case 110: x_III=t_III; break;
        default: x_III=s_III; break;
    }
    
        return ( std::log(
                     Complex(
                                  (
                                   std::sqrt(Complex(x_III-s))
                                   +std::sqrt(Complex(x_III-min))
                                   )/(
                                      std::sqrt(Complex(x_III-s))
                                      -std::sqrt(Complex(x_III-min))
                                      )
                                  )
                     )
                -2.i*std::atan(
                          Complex(
                                       std::sqrt(Complex(max-x_III))
                                       /std::sqrt(Complex(x_III-s))
                                       )
                          )
                )/std::sqrt(Complex(x_III-s));
}

// the order of the R integral is determined by (2n+1)/2
Complex Dispersion::analytic_integral_R(Complex s, double min, double max, Setting evaluation, int n)const{
    double x_III;
    // t-channel in eta'-> eta pi pi has a different pseudo-threshold t_III
    switch (isospin) {
        case 010: case 110: x_III=t_III; break;
        default: x_III=s_III; break;
    }
    
    if (n==0){
        return analytic_integral_R12(s, min, max, evaluation);
    }
    else if (n<0){
        throw std::domain_error{"Order n of R-type integral needs to be positive."};
    }
    else{
        return (
                2.*(
                    1./((2.*n-1.)*pow(std::sqrt(Complex(x_III-max)),(2.*n-1.)))
                    -1./((2.*n-1.)*pow(std::sqrt(Complex(x_III-min)),(2.*n-1.)))
                    )
                +analytic_integral_R(s, min, max, evaluation, n-1)
                )/(Complex(x_III-s));
    }
}


// the order of the Q integral is determined by (2n+1)/2
Complex Dispersion::analytic_integral_Q(Complex s, double min, double max, int n)const{
    double x_III;
    // t-channel in eta'-> eta pi pi has a different pseudo-threshold t_III
    switch (isospin) {
        case 010: case 110: x_III=t_III; break;
        default: x_III=s_III; break;
    }

    if (n==0){
        return analytic_integral_Q12(s, min, max);
    }
    else if (n<0){
        throw std::domain_error{"Order n of Q-type integral needs to be positive."};
    }
    else{
        return (
            -2.*(
                1.i/((2.*n-1.)*pow(std::sqrt(Complex(max-x_III)),(2.*n-1.)))
                +1./((2.*n-1.)*pow(std::sqrt(Complex(x_III-min)),(2.*n-1.)))
                )
                +analytic_integral_Q(s, min, max, n-1)
                )/(Complex(x_III-s));
    }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------




// For real-valued s with infinitesimal imaginary part one can split the
// dispersion integral in the following way (not a unique representation).
// This ensures a numerically stable result.
Complex DispersionEta3Pi::disp_integral_below_pseudo(Complex s,  bool adaptive, Setting evaluation)const{
    double x=real(s);
    double x_I=s_I;
    double x_III=s_III;
    
    double p=(x_III+x)/2.0;
    // For S-waves
    if (isospin==0 || isospin==2) {
        return numeric_integral_cauchy(x, x_I, p, adaptive)
        +numeric_integral_pseudo(x, p, cutoff, adaptive)
        +numerator(x) * analytic_integral_R(x, x_I, p, evaluation, 0)
        +numerator(x_III) * analytic_integral_Q(x, p, cutoff, 0)
        ;
    }
    // For P-waves
    if (isospin==1) {
        return numeric_integral_cauchy(x, x_I, p, adaptive)
        +numeric_integral_pseudo(x, p, cutoff, adaptive)
        +numerator(x) * analytic_integral_R(x, x_I, p, evaluation, 1)
        +numerator(x_III) * analytic_integral_Q(x, p, cutoff, 1)
        +numerator_1_derivative_s_below * analytic_integral_Q(x, p, cutoff, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 0, 1 or 2."};
}

Complex DispersionEtap3Pi::disp_integral_below_pseudo(Complex s,  bool adaptive, Setting evaluation)const{
    double x=real(s);
    double x_I=s_I;
    double x_III=s_III;
    
    double p=(x_III+x)/2.0;
    // For S-waves
    if (isospin==0 || isospin==2) {
        return numeric_integral_cauchy(x, x_I, p, adaptive)
        +numeric_integral_pseudo(x, p, cutoff, adaptive)
        +numerator(x) * analytic_integral_R(x, x_I, p, evaluation, 0)
        +numerator(x_III) * analytic_integral_Q(x, p, cutoff, 0)
        ;
    }
    // For P-waves
    if (isospin==1) {
        return numeric_integral_cauchy(x, x_I, p, adaptive)
        +numeric_integral_pseudo(x, p, cutoff, adaptive)
        +numerator(x) * analytic_integral_R(x, x_I, p, evaluation, 1)
        +numerator(x_III) * analytic_integral_Q(x, p, cutoff, 1)
        +numerator_1_derivative_s_below * analytic_integral_Q(x, p, cutoff, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 0, 1 or 2."};
}

Complex DispersionEtapEtaPiPi::disp_integral_below_pseudo(Complex s,  bool adaptive, Setting evaluation)const{
    double x=real(s);
    double x_I;
    double x_III;
    // t-channel in eta'-> eta pi pi has a different thresholds t_I and t_III
    switch (isospin) {
        case 010: case 110: x_III=t_III; x_I=t_I; break;
        default: x_III=s_III; x_I=s_I; break;
    }
    
    double p=(x_III+x)/2.0;
    // For S-waves
    if (isospin==000 || isospin==010 || isospin==110) {
        return numeric_integral_cauchy(x, x_I, p, adaptive)
        +numeric_integral_pseudo(x, p, cutoff, adaptive)
        +numerator(x) * analytic_integral_R(x, x_I, p, evaluation, 0)
        +numerator(x_III) * analytic_integral_Q(x, p, cutoff, 0)
        ;
    }
    // For P-waves
    if (isospin==111) {
        return numeric_integral_cauchy(x, x_I, p, adaptive)
        +numeric_integral_pseudo(x, p, cutoff, adaptive)
        +numerator(x) * analytic_integral_R(x, x_I, p, evaluation, 1)
        +numerator(x_III) * analytic_integral_Q(x, p, cutoff, 1)
        +numerator_1_derivative_s_below * analytic_integral_Q(x, p, cutoff, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 0, 1 or 2."};
}

Complex DispersionV3Pi::disp_integral_below_pseudo(Complex s,  bool adaptive, Setting evaluation)const{
    double x=real(s);
    double x_I=s_I;
    double x_III=s_III;
    
    double p=(x_III+x)/2.0;
    // For P-waves
    if (isospin==1) {
        return numeric_integral_cauchy(x, x_I, p, adaptive)
        +numeric_integral_pseudo(x, p, cutoff, adaptive)
        +numerator(x) * analytic_integral_R(x, x_I, p, evaluation, 1)
        +numerator(x_III) * analytic_integral_Q(x, p, cutoff, 1)
        +numerator_1_derivative_s_below * analytic_integral_Q(x, p, cutoff, 0)
        ;
    }
    throw std::domain_error{"Isospin should be 1."};
}

Complex DispersionX3Pi::disp_integral_below_pseudo(Complex s,  bool adaptive, Setting evaluation)const{
    double x=real(s);
    double x_I=s_I;
    double x_III=s_III;
    
    double p=(x_III+x)/2.0;
    // For P-waves
    if (isospin==1) {
        return numeric_integral_cauchy(x, x_I, p, adaptive)
        +numeric_integral_pseudo(x, p, cutoff, adaptive)
        +numerator(x) * analytic_integral_R(x, x_I, p, evaluation, 1)
        +numerator(x_III) * analytic_integral_Q(x, p, cutoff, 1)
        +numerator_1_derivative_s_below * analytic_integral_Q(x, p, cutoff, 0)
        ;
    }
    throw std::domain_error{"Isospin should be 1."};
}

Complex DispersionT3Pi::disp_integral_below_pseudo(Complex s,  bool adaptive, Setting evaluation)const{
    double x=real(s);
    double x_I=s_I;
    double x_III=s_III;
    
    double p=(x_III+x)/2.0;
    // For P-waves
    if (isospin==1) {
        return numeric_integral_cauchy(x, x_I, p, adaptive)
        +numeric_integral_pseudo(x, p, cutoff, adaptive)
        +numerator(x) * analytic_integral_R(x, x_I, p, evaluation, 2)
        +numerator(x_III) * analytic_integral_Q(x, p, cutoff, 2)
        +numerator_1_derivative_s_below * analytic_integral_Q(x, p, cutoff, 1)
        +numerator_2_derivative_s_below * analytic_integral_Q(x, p, cutoff, 0)
        ;
    }
    throw std::domain_error{"Isospin should be 1."};
}

Complex DispersionEta3Pi::disp_integral_above_pseudo(Complex s,  bool adaptive, Setting evaluation)const{
    double x=real(s);
    double x_I=s_I;
    double x_III=s_III;
    
    double p=(x_III+x)/2.0;
    // For S-waves
    if (isospin==0 || isospin==2) {
        return numeric_integral_cauchy(x, p, cutoff, adaptive)
        +numeric_integral_pseudo(x, x_I, p, adaptive)
        +numerator(x_III) * analytic_integral_Q(x, x_I, p, 0)
        +numerator(x) * analytic_integral_R(x, p, cutoff, evaluation, 0)
        ;
    }
    // For P-waves
    if (isospin==1) {
        return numeric_integral_cauchy(x, p, cutoff, adaptive)
        +numeric_integral_pseudo(x, x_I, p, adaptive)
        +numerator(x) * analytic_integral_R(x, p, cutoff, evaluation, 1)
        +numerator(x_III) * analytic_integral_Q(x, x_I, p, 1)
        +numerator_1_derivative_s_above * analytic_integral_Q(x, x_I, p, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 0, 1 or 2."};
}

Complex DispersionEtap3Pi::disp_integral_above_pseudo(Complex s,  bool adaptive, Setting evaluation)const{
    double x=real(s);
    double x_I=s_I;
    double x_III=s_III;
    
    double p=(x_III+x)/2.0;
    // For S-waves
    if (isospin==0 || isospin==2) {
        return numeric_integral_cauchy(x, p, cutoff, adaptive)
        +numeric_integral_pseudo(x, x_I, p, adaptive)
        +numerator(x_III) * analytic_integral_Q(x, x_I, p, 0)
        +numerator(x) * analytic_integral_R(x, p, cutoff, evaluation, 0)
        ;
    }
    // For P-waves
    if (isospin==1) {
        return numeric_integral_cauchy(x, p, cutoff, adaptive)
        +numeric_integral_pseudo(x, x_I, p, adaptive)
        +numerator(x) * analytic_integral_R(x, p, cutoff, evaluation, 1)
        +numerator(x_III) * analytic_integral_Q(x, x_I, p, 1)
        +numerator_1_derivative_s_above * analytic_integral_Q(x, x_I, p, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 0, 1 or 2."};
}

Complex DispersionEtapEtaPiPi::disp_integral_above_pseudo(Complex s,  bool adaptive, Setting evaluation)const{
    double x=real(s);
    double x_I;
    double x_III;
    // t-channel in eta'-> eta pi pi has a different pseudo-threshold t_III
    switch (isospin) {
        case 010: case 110: x_III=t_III; x_I=t_I; break;
        default: x_III=s_III; x_I=s_I; break;
    }
    
    double p=(x_III+x)/2.0;
    // For S-waves
    if (isospin==000 || isospin==010 || isospin==110) {
        return numeric_integral_cauchy(x, p, cutoff, adaptive)
        +numeric_integral_pseudo(x, x_I, p, adaptive)
        +numerator(x_III) * analytic_integral_Q(x, x_I, p, 0)
        +numerator(x) * analytic_integral_R(x, p, cutoff, evaluation, 0)
        ;
    }
    // For P-waves
    if (isospin==111) {
        return numeric_integral_cauchy(x, p, cutoff, adaptive)
        +numeric_integral_pseudo(x, x_I, p, adaptive)
        +numerator(x) * analytic_integral_R(x, p, cutoff, evaluation, 1)
        +numerator(x_III) * analytic_integral_Q(x, x_I, p, 1)
        +numerator_1_derivative_s_above * analytic_integral_Q(x, x_I, p, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 0, 1 or 2."};
}

Complex DispersionV3Pi::disp_integral_above_pseudo(Complex s,  bool adaptive, Setting evaluation)const{
    double x=real(s);
    double x_I=s_I;
    double x_III=s_III;
    
    double p=(x_III+x)/2.0;
    // For P-waves
    if (isospin==1) {
        return numeric_integral_cauchy(x, p, cutoff, adaptive)
        +numeric_integral_pseudo(x, x_I, p, adaptive)
        +numerator(x) * analytic_integral_R(x, p, cutoff, evaluation, 1)
        +numerator(x_III) * analytic_integral_Q(x, x_I, p, 1)
        +numerator_1_derivative_s_above * analytic_integral_Q(x, x_I, p, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 1."};
}

Complex DispersionX3Pi::disp_integral_above_pseudo(Complex s,  bool adaptive, Setting evaluation)const{
    double x=real(s);
    double x_I=s_I;
    double x_III=s_III;
    
    double p=(x_III+x)/2.0;
    // For P-waves
    if (isospin==1) {
        return numeric_integral_cauchy(x, p, cutoff, adaptive)
        +numeric_integral_pseudo(x, x_I, p, adaptive)
        +numerator(x) * analytic_integral_R(x, p, cutoff, evaluation, 1)
        +numerator(x_III) * analytic_integral_Q(x, x_I, p, 1)
        +numerator_1_derivative_s_above * analytic_integral_Q(x, x_I, p, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 1."};
}

Complex DispersionT3Pi::disp_integral_above_pseudo(Complex s,  bool adaptive, Setting evaluation)const{
    double x=real(s);
    double x_I=s_I;
    double x_III=s_III;
    
    double p=(x_III+x)/2.0;
    // For P-waves
    if (isospin==1) {
        return numeric_integral_cauchy(x, p, cutoff, adaptive)
        +numeric_integral_pseudo(x, x_I, p, adaptive)
        +numerator(x) * analytic_integral_R(x, p, cutoff, evaluation, 2)
        +numerator(x_III) * analytic_integral_Q(x, x_I, p, 2)
        +numerator_1_derivative_s_above * analytic_integral_Q(x, x_I, p, 1)
        +numerator_2_derivative_s_above * analytic_integral_Q(x, x_I, p, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 1."};
}


Complex DispersionEta3Pi::disp_integral_trivial(Complex s,  bool adaptive)const{
    double x_I=s_I;
    double x_III=s_III;
    
    // For S-waves
    if (isospin==0 || isospin==2) {
        return numeric_integral_pseudo(s, x_I, cutoff, adaptive)
        +numerator(x_III) * analytic_integral_Q(s, x_I, cutoff, 0);
    }
    // For P-waves
    if (isospin==1) {
        return numeric_integral_pseudo(s, x_I, cutoff, adaptive)
        +numerator(x_III) * analytic_integral_Q(s, x_I, cutoff, 1)
        +numerator_1_derivative_s_below * analytic_integral_Q(s, x_I, cutoff, 0);
    }
    throw std::domain_error{"Isospin should be either 0, 1 or 2."};
}

Complex DispersionEtap3Pi::disp_integral_trivial(Complex s,  bool adaptive)const{
    double x_I=s_I;
    double x_III=s_III;
    
    // For S-waves
    if (isospin==0 || isospin==2) {
        return numeric_integral_pseudo(s, x_I, cutoff, adaptive)
        +numerator(x_III) * analytic_integral_Q(s, x_I, cutoff, 0);
    }
    // For P-waves
    if (isospin==1) {
        return numeric_integral_pseudo(s, x_I, cutoff, adaptive)
        +numerator(x_III) * analytic_integral_Q(s, x_I, cutoff, 1)
        +numerator_1_derivative_s_below * analytic_integral_Q(s, x_I, cutoff, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 0, 1 or 2."};
}

Complex DispersionEtapEtaPiPi::disp_integral_trivial(Complex s,  bool adaptive)const{
    double x_I; 
    double x_III;
    // t-channel in eta'-> eta pi pi has a different pseudo-threshold t_III
    switch (isospin) {
        case 010: case 110: x_III=t_III; x_I=t_I; break;
        default: x_III=s_III; x_I=s_I; break;
    }
    
    // For S-waves
    if (isospin==000 || isospin==010 || isospin==110) {
        return numeric_integral_pseudo(s, x_I, cutoff, adaptive)
        +numerator(x_III) * analytic_integral_Q(s, x_I, cutoff, 0);
    }
    // For P-waves
    if (isospin==111) {
        return numeric_integral_pseudo(s, x_I, cutoff, adaptive)
        +numerator(x_III) * analytic_integral_Q(s, x_I, cutoff, 1)
        +numerator_1_derivative_s_below * analytic_integral_Q(s, x_I, cutoff, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 0, 1 or 2."};
}

Complex DispersionV3Pi::disp_integral_trivial(Complex s,  bool adaptive)const{
    double x_I=s_I;
    double x_III=s_III;

    // For P-waves
    if (isospin==1) {
        return numeric_integral_pseudo(s, x_I, cutoff, adaptive)
        +numerator(x_III) * analytic_integral_Q(s, x_I, cutoff, 1)
        +numerator_1_derivative_s_below * analytic_integral_Q(s, x_I, cutoff, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 1."};
}

Complex DispersionX3Pi::disp_integral_trivial(Complex s,  bool adaptive)const{
    double x_I=s_I;
    double x_III=s_III;

    // For P-waves
    if (isospin==1) {
        return numeric_integral_pseudo(s, x_I, cutoff, adaptive)
        +numerator(x_III) * analytic_integral_Q(s, x_I, cutoff, 1)
        +numerator_1_derivative_s_below * analytic_integral_Q(s, x_I, cutoff, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 1."};
}

Complex DispersionT3Pi::disp_integral_trivial(Complex s,  bool adaptive)const{
    double x_I=s_I;
    double x_III=s_III;

    // For P-waves
    if (isospin==1) {
        return numeric_integral_pseudo(s, x_I, cutoff, adaptive)
        +numerator(x_III) * analytic_integral_Q(s, x_I, cutoff, 2)
        +numerator_1_derivative_s_below * analytic_integral_Q(s, x_I, cutoff, 1)
        +numerator_2_derivative_s_below * analytic_integral_Q(s, x_I, cutoff, 0)
        ;
    }
    throw std::domain_error{"Isospin should be either 1."};
}



//------------------------------------------------------------------------------------
//-- Dispersion integral -------------------------------------------------------------
//------------------------------------------------------------------------------------

Complex DispersionEta3Pi::dispersion_integral(Complex s, Setting evaluation)const{
    // Increase the precision of the integration close to critical points.
    // Feel free to adjust the range around each point within which the high-precision
    // shall be used.
    
    double x=real(s);
    
    if (isospin!=0 && isospin!=1 && isospin!=2) {
        throw std::domain_error{"Value for isospin is not allowed."};
    }
    
    double x_I=s_I;
    double x_III=s_III;

    // To obtain the best trade-off between speed and precision one can:
    // 1) play around with 'eps_pseudo' and 'eps_cusp' (i.e. define where to use the adaptive integration routines)
    // 2) switch between 'true' and 'false' in 'disp_integral_above_pseudo' and 'disp_integral_below_pseudo'
    double eps_pseudo=2;
    double eps_cusp=1;
    
    switch (evaluation) {
        case Setting::above: case Setting::below:
            if (x>=x_I-eps_cusp) {
                if (x<x_III-eps_pseudo) {
                    return
                      disp_integral_below_pseudo(s, true, evaluation);
                }
                // use a high precision close to pseudo-threshold x_III
                if (x>=x_III-eps_pseudo&& x<=x_III) {
                    return
                    disp_integral_below_pseudo(s, true, evaluation);
                }
                // use a high precision close to pseudo-threshold x_III
                if (x>x_III && x<=x_III+eps_pseudo){
                    return
                    disp_integral_above_pseudo(s, true, evaluation);
                }
                // use a high precision close to cusp in the phase shift
                if (isospin==0) {
                    if (x>x_III+eps_pseudo && x<cusp-eps_cusp){
                        return
                        disp_integral_above_pseudo(s, false, evaluation);
                    }
                    if (x>=cusp-eps_cusp && x<=cusp+3.*eps_cusp){
                        return
                        disp_integral_above_pseudo(s, true, evaluation);
                    }
                    if (x>cusp+3.*eps_cusp){
                        return
                        disp_integral_above_pseudo(s, false, evaluation);
                    }
                    throw std::domain_error{"Value for Mandelstam 's' not allowed."};
                }
                // no cusp in phase shifts for these isospins
                if (isospin==1 || isospin==2) {
                    if (x>x_III+eps_pseudo){
                        return
                        disp_integral_above_pseudo(s, false, evaluation);
                    }
                    throw std::domain_error{"Value for Mandelstam 's' not allowed."};
                }
                throw std::domain_error{"Value for Mandelstam 's' not allowed."};
            }
            // for s below threshold one can evaluate the integral straightforwardly
            else{
                return
                disp_integral_trivial(s, false);
            }
            break;
            
            // for complex s one can evaluate the integral straightforwardly
        default:
            return
            disp_integral_trivial(s, false);
            break;
    }
}

Complex DispersionEtap3Pi::dispersion_integral(Complex s, Setting evaluation)const{
    // Increase the precision of the integration close to critical points.
    // Feel free to adjust the range around each point within which the high-precision
    // shall be used.
    
    double x=real(s);
    
    if (isospin!=0 && isospin!=1 && isospin!=2) {
        throw std::domain_error{"Value for isospin is not allowed."};
    }
    
    double x_I=s_I;
    double x_III=s_III;

    // To obtain the best trade-off between speed and precision one can:
    // 1) play around with 'eps_pseudo' and 'eps_cusp' (i.e. define where to use the adaptive integration routines)
    // 2) switch between 'true' and 'false' in 'disp_integral_above_pseudo' and 'disp_integral_below_pseudo'
    double eps_pseudo=2;
    double eps_cusp=1;
    
    switch (evaluation) {
        case Setting::above: case Setting::below:
            if (x>=x_I-eps_cusp) {
                if (x<x_III-eps_pseudo) {
                    return
                      disp_integral_below_pseudo(s, true, evaluation);
                }
                // use a high precision close to pseudo-threshold x_III
                if (x>=x_III-eps_pseudo&& x<=x_III) {
                    return
                    disp_integral_below_pseudo(s, true, evaluation);
                }
                // use a high precision close to pseudo-threshold x_III
                if (x>x_III && x<=x_III+eps_pseudo){
                    return
                    disp_integral_above_pseudo(s, true, evaluation);
                }
                // use a high precision close to cusp in the phase shift
                if (isospin==0) {
                    if (x>x_III+eps_pseudo && x<cusp-eps_cusp){
                        return
                        disp_integral_above_pseudo(s, false, evaluation);
                    }
                    if (x>=cusp-eps_cusp && x<=cusp+3.*eps_cusp){
                        return
                        disp_integral_above_pseudo(s, true, evaluation);
                    }
                    if (x>cusp+3.*eps_cusp){
                        return
                        disp_integral_above_pseudo(s, false, evaluation);
                    }
                    throw std::domain_error{"Value for Mandelstam 's' not allowed."};
                }
                // no cusp in phase shifts for these isospins
                if (isospin==1 || isospin==2) {
                    if (x>x_III+eps_pseudo){
                        return
                        disp_integral_above_pseudo(s, false, evaluation);
                    }
                    throw std::domain_error{"Value for Mandelstam 's' not allowed."};
                }
                throw std::domain_error{"Value for Mandelstam 's' not allowed."};
            }
            // for s below threshold one can evaluate the integral straightforwardly
            else{
                return
                disp_integral_trivial(s, false);
            }
            break;
            
            // for complex s one can evaluate the integral straightforwardly
        default:
            return
            disp_integral_trivial(s, false);
            break;
    }
}

Complex DispersionEtapEtaPiPi::dispersion_integral(Complex s, Setting evaluation)const{
    // Increase the precision of the integration close to critical points.
    // Feel free to adjust the range around each point within which the high-precision
    // shall be used.
    
    double x=real(s);
    
    if (isospin!=000 && isospin!=010 && isospin!=110 && isospin!=111) {
        throw std::domain_error{"Value for isospin is not allowed."};
    }
    
    double x_I;
    double x_III;
    // t-channel in eta'-> eta pi pi has different thresholds t_I and t_III
    switch (isospin) {
        case 010: case 110: x_III=t_III; x_I=t_I; break;
        default: x_III=s_III; x_I=s_I; break;
    }
    
    // To obtain the best trade-off between speed and precision one can:
    // 1) play around with 'eps_pseudo' and 'eps_cusp' (i.e. define where to use the adaptive integration routines)
    // 2) switch between 'true' and 'false' in 'disp_integral_above_pseudo' and 'disp_integral_below_pseudo'
    double eps_pseudo = 5;
    double eps_cusp=1;
    
    switch (evaluation) {
        case Setting::above: case Setting::below:
            if (x>=x_I-eps_cusp) {
                if (x<x_III-eps_pseudo) {
                    return
                      disp_integral_below_pseudo(s, true, evaluation);
                }
                // use a high precision close to pseudo-threshold x_III
                if (x>=x_III-eps_pseudo&& x<=x_III) {
                    return
                    disp_integral_below_pseudo(s, true, evaluation);
                }
                // use a high precision close to pseudo-threshold x_III
                if (x>x_III && x<=x_III+eps_pseudo){
                    return
                    disp_integral_above_pseudo(s, true, evaluation);
                }
                // use a high precision close to cusp in the phase shift
                if (isospin==000 || isospin==010 || isospin==110) {
                    if (x>x_III+eps_pseudo && x<cusp-eps_cusp){
                        return
                        disp_integral_above_pseudo(s, false, evaluation);
                    }
                    if (x>=cusp-eps_cusp && x<=cusp+3.*eps_cusp){
                        return
                        disp_integral_above_pseudo(s, true, evaluation);
                    }
                    if (x>cusp+3.*eps_cusp){
                        return
                        disp_integral_above_pseudo(s, false, evaluation);
                    }
                    throw std::domain_error{"Value for Mandelstam 's' not allowed."};
                }
                // no cusp in phase shifts for these isospins
                if (isospin==111) {
                    if (x>x_III+eps_pseudo){
                        return
                        disp_integral_above_pseudo(s, false, evaluation);
                    }
                    throw std::domain_error{"Value for Mandelstam 's' not allowed."};
                }
                throw std::domain_error{"Value for Mandelstam 's' not allowed."};
            }
            // for s below threshold one can evaluate the integral straightforwardly
            else{
                return
                disp_integral_trivial(s, false);
            }
            break;
            
            // for complex s one can evaluate the integral straightforwardly
        default:
            return
            disp_integral_trivial(s, false);
            break;
    }
}

Complex DispersionV3Pi::dispersion_integral(Complex s, Setting evaluation)const{
    // Increase the precision of the integration close to critical points.
    // Feel free to adjust the range around each point within which the high-precision
    // shall be used.
    
    double x=real(s);
    
    if (isospin!=1) {
        throw std::domain_error{"Value for isospin is not allowed."};
    }
    
    double x_I=s_I;
    double x_III=s_III;

    // To obtain the best trade-off between speed and precision one can:
    // 1) play around with 'eps_pseudo' and 'eps_cusp' (i.e. define where to use the adaptive integration routines)
    // 2) switch between 'true' and 'false' in 'disp_integral_above_pseudo' and 'disp_integral_below_pseudo'
    double eps_pseudo=2;
    double eps_cusp=1;
    
    switch (evaluation) {
        case Setting::above: case Setting::below:
            if (x>=x_I-eps_cusp) {
                if (x<x_III-eps_pseudo) {
                    return
                      disp_integral_below_pseudo(s, true, evaluation);
                }
                // use a high precision close to pseudo-threshold x_III
                if (x>=x_III-eps_pseudo&& x<=x_III) {
                    return
                    disp_integral_below_pseudo(s, true, evaluation);
                }
                // use a high precision close to pseudo-threshold x_III
                if (x>x_III && x<=x_III+eps_pseudo){
                    return
                    disp_integral_above_pseudo(s, true, evaluation);
                }
                // no cusp in phase shifts for these isospins
                if (isospin==1) {
                    if (x>x_III+eps_pseudo){
                        return
                        disp_integral_above_pseudo(s, false, evaluation);
                    }
                    throw std::domain_error{"Value for Mandelstam 's' not allowed."};
                }
                throw std::domain_error{"Value for Mandelstam 's' not allowed."};
            }
            // for s below threshold one can evaluate the integral straightforwardly
            else{
                return
                disp_integral_trivial(s, false);
            }
            break;
            
            // for complex s one can evaluate the integral straightforwardly
        default:
            return
            disp_integral_trivial(s, false);
            break;
    }
}

Complex DispersionX3Pi::dispersion_integral(Complex s, Setting evaluation)const{
    // Increase the precision of the integration close to critical points.
    // Feel free to adjust the range around each point within which the high-precision
    // shall be used.
    
    double x=real(s);
    
    if (isospin!=1) {
        throw std::domain_error{"Value for isospin is not allowed."};
    }
    
    double x_I=s_I;
    double x_III=s_III;

    // To obtain the best trade-off between speed and precision one can:
    // 1) play around with 'eps_pseudo' and 'eps_cusp' (i.e. define where to use the adaptive integration routines)
    // 2) switch between 'true' and 'false' in 'disp_integral_above_pseudo' and 'disp_integral_below_pseudo'
    double eps_pseudo=2;
    double eps_cusp=1;
    
    switch (evaluation) {
        case Setting::above: case Setting::below:
            if (x>=x_I-eps_cusp) {
                if (x<x_III-eps_pseudo) {
                    return
                      disp_integral_below_pseudo(s, true, evaluation);
                }
                // use a high precision close to pseudo-threshold x_III
                if (x>=x_III-eps_pseudo&& x<=x_III) {
                    return
                    disp_integral_below_pseudo(s, true, evaluation);
                }
                // use a high precision close to pseudo-threshold x_III
                if (x>x_III && x<=x_III+eps_pseudo){
                    return
                    disp_integral_above_pseudo(s, true, evaluation);
                }
                // no cusp in phase shifts for these isospins
                if (isospin==1) {
                    if (x>x_III+eps_pseudo){
                        return
                        disp_integral_above_pseudo(s, false, evaluation);
                    }
                    throw std::domain_error{"Value for Mandelstam 's' not allowed."};
                }
                throw std::domain_error{"Value for Mandelstam 's' not allowed."};
            }
            // for s below threshold one can evaluate the integral straightforwardly
            else{
                return
                disp_integral_trivial(s, false);
            }
            break;
            
            // for complex s one can evaluate the integral straightforwardly
        default:
            return
            disp_integral_trivial(s, false);
            break;
    }
}

Complex DispersionT3Pi::dispersion_integral(Complex s, Setting evaluation)const{
    // Increase the precision of the integration close to critical points.
    // Feel free to adjust the range around each point within which the high-precision
    // shall be used.
    
    double x=real(s);
    
    if (isospin!=1) {
        throw std::domain_error{"Value for isospin is not allowed."};
    }
    
    double x_I=s_I;
    double x_III=s_III;

    // To obtain the best trade-off between speed and precision one can:
    // 1) play around with 'eps_pseudo' and 'eps_cusp' (i.e. define where to use the adaptive integration routines)
    // 2) switch between 'true' and 'false' in 'disp_integral_above_pseudo' and 'disp_integral_below_pseudo'
    double eps_pseudo=2;
    double eps_cusp=1;
    
    switch (evaluation) {
        case Setting::above: case Setting::below:
            if (x>=x_I-eps_cusp) {
                if (x<x_III-eps_pseudo) {
                    return
                      disp_integral_below_pseudo(s, true, evaluation);
                }
                // use a high precision close to pseudo-threshold x_III
                if (x>=x_III-eps_pseudo&& x<=x_III) {
                    return
                    disp_integral_below_pseudo(s, true, evaluation);
                }
                // use a high precision close to pseudo-threshold x_III
                if (x>x_III && x<=x_III+eps_pseudo){
                    return
                    disp_integral_above_pseudo(s, true, evaluation);
                }
                // no cusp in phase shifts for these isospins
                if (isospin==1) {
                    if (x>x_III+eps_pseudo){
                        return
                        disp_integral_above_pseudo(s, false, evaluation);
                    }
                    throw std::domain_error{"Value for Mandelstam 's' not allowed."};
                }
                throw std::domain_error{"Value for Mandelstam 's' not allowed."};
            }
            // for s below threshold one can evaluate the integral straightforwardly
            else{
                return
                disp_integral_trivial(s, false);
            }
            break;
            
            // for complex s one can evaluate the integral straightforwardly
        default:
            return
            disp_integral_trivial(s, false);
            break;
    }
}
//-----------------------------------------------------------------
//-----------------------------------------------------------------




//-----------------------------------------------------------------------------------
//-- Build lists --------------------------------------------------------------------
//-----------------------------------------------------------------------------------
void DispersionEta3Pi::build_dispersion_integral_lists(){
    for (std::size_t j=0; j<mandelstam_list.interval_disp.size(); j++) {
        dispersion_integral_above.push_back(dispersion_integral(mandelstam_list.interval_disp[j]+1.i*infinitesimal, Setting::above));
        dispersion_integral_below.push_back(dispersion_integral(mandelstam_list.interval_disp[j]-1.i*infinitesimal, Setting::below));
    }
    for (std::size_t j=0; j<mandelstam_list.interval_phi.size(); j++) {
        dispersion_integral_egg.push_back(dispersion_integral(polar_egg(mandelstam_list.interval_phi[j]), Setting::egg));
    }
    return;
}

void DispersionEtap3Pi::build_dispersion_integral_lists(){
    for (std::size_t j=0; j<mandelstam_list.interval_disp.size(); j++) {
        dispersion_integral_above.push_back(dispersion_integral(mandelstam_list.interval_disp[j]+1.i*infinitesimal, Setting::above));
        dispersion_integral_below.push_back(dispersion_integral(mandelstam_list.interval_disp[j]-1.i*infinitesimal, Setting::below));
    }
    for (std::size_t j=0; j<mandelstam_list.interval_phi.size(); j++) {
        dispersion_integral_egg.push_back(dispersion_integral(polar_egg(mandelstam_list.interval_phi[j]), Setting::egg));
    }
    return;
}

void DispersionEtapEtaPiPi::build_dispersion_integral_lists(){
    for (std::size_t j=0; j<mandelstam_list.interval_disp.size(); j++) {
        dispersion_integral_above.push_back(dispersion_integral(mandelstam_list.interval_disp[j]+1.i*infinitesimal, Setting::above));
        dispersion_integral_below.push_back(dispersion_integral(mandelstam_list.interval_disp[j]-1.i*infinitesimal, Setting::below));
    }    
    for (std::size_t j=0; j<mandelstam_list.interval_t_disp.size(); j++) {
        dispersion_integral_t_channel_above.push_back(dispersion_integral(mandelstam_list.interval_t_disp[j]+1.i*infinitesimal, Setting::above));
        dispersion_integral_t_channel_below.push_back(dispersion_integral(mandelstam_list.interval_t_disp[j]-1.i*infinitesimal, Setting::below));
    }
    for (std::size_t j=0; j<mandelstam_list.interval_y.size(); ++j) {
        double y=mandelstam_list.interval_y[j];
        

        dispersion_integral_plus_upper_a.push_back(dispersion_integral(gamma(y,Setting::plus_upper_a),Setting::plus_upper_a));
        dispersion_integral_plus_upper_b.push_back(dispersion_integral(gamma(y,Setting::plus_upper_b),Setting::plus_upper_b));
        dispersion_integral_plus_lower_a.push_back(dispersion_integral(gamma(y,Setting::plus_lower_a),Setting::plus_lower_a));
        dispersion_integral_plus_lower_b.push_back(dispersion_integral(gamma(y,Setting::plus_lower_b),Setting::plus_lower_b));

        dispersion_integral_zero_upper_a.push_back(dispersion_integral(gamma(y,Setting::zero_upper_a),Setting::zero_upper_a));
        dispersion_integral_zero_upper_b.push_back(dispersion_integral(gamma(y,Setting::zero_upper_b),Setting::zero_upper_b));
        dispersion_integral_zero_lower_a.push_back(dispersion_integral(gamma(y,Setting::zero_lower_a),Setting::zero_lower_a));
        dispersion_integral_zero_lower_b.push_back(dispersion_integral(gamma(y,Setting::zero_lower_b),Setting::zero_lower_b));

        dispersion_integral_minus_upper_a.push_back(dispersion_integral(gamma(y,Setting::minus_upper_a),Setting::minus_upper_a));
        dispersion_integral_minus_upper_b.push_back(dispersion_integral(gamma(y,Setting::minus_upper_b),Setting::minus_upper_b));
        dispersion_integral_minus_lower_a.push_back(dispersion_integral(gamma(y,Setting::minus_lower_a),Setting::minus_lower_a));
        dispersion_integral_minus_lower_b.push_back(dispersion_integral(gamma(y,Setting::minus_lower_b),Setting::minus_lower_b));

    }
    return;
}

void DispersionV3Pi::build_dispersion_integral_lists(){
    for (std::size_t j=0; j<mandelstam_list.interval_disp.size(); j++) {
        dispersion_integral_above.push_back(dispersion_integral(mandelstam_list.interval_disp[j]+1.i*infinitesimal, Setting::above));
        dispersion_integral_below.push_back(dispersion_integral(mandelstam_list.interval_disp[j]-1.i*infinitesimal, Setting::below));
    }
    for (std::size_t j=0; j<mandelstam_list.interval_phi.size(); j++) {
        dispersion_integral_egg.push_back(dispersion_integral(polar_egg(mandelstam_list.interval_phi[j]), Setting::egg));
    }
    return;
}

void DispersionX3Pi::build_dispersion_integral_lists(){
    for (std::size_t j=0; j<mandelstam_list.interval_disp.size(); j++) {
        dispersion_integral_above.push_back(dispersion_integral(mandelstam_list.interval_disp[j]+1.i*infinitesimal, Setting::above));
        dispersion_integral_below.push_back(dispersion_integral(mandelstam_list.interval_disp[j]-1.i*infinitesimal, Setting::below));
    }
    for (std::size_t j=0; j<mandelstam_list.interval_phi.size(); j++) {
        dispersion_integral_egg.push_back(dispersion_integral(polar_egg(mandelstam_list.interval_phi[j]), Setting::egg));
    }
    return;
}

void DispersionT3Pi::build_dispersion_integral_lists(){
    for (std::size_t j=0; j<mandelstam_list.interval_disp.size(); j++) {
        dispersion_integral_above.push_back(dispersion_integral(mandelstam_list.interval_disp[j]+1.i*infinitesimal, Setting::above));
        dispersion_integral_below.push_back(dispersion_integral(mandelstam_list.interval_disp[j]-1.i*infinitesimal, Setting::below));
    }
    for (std::size_t j=0; j<mandelstam_list.interval_phi.size(); j++) {
        dispersion_integral_egg.push_back(dispersion_integral(polar_egg(mandelstam_list.interval_phi[j]), Setting::egg));
    }
    return;
}
//-----------------------------------------------------------------
//-----------------------------------------------------------------

void DispersionEta3Pi::build_splines(){
    spline_disp_above= cauchy::Interpolate(mandelstam_list.interval_disp,dispersion_integral_above,gsl::InterpolationMethod::cubic);
    spline_disp_below= cauchy::Interpolate(mandelstam_list.interval_disp,dispersion_integral_below,gsl::InterpolationMethod::cubic);
    spline_disp_egg= cauchy::Interpolate(mandelstam_list.interval_phi,dispersion_integral_egg,gsl::InterpolationMethod::cubic);
}

void DispersionEtap3Pi::build_splines(){
    spline_disp_above= cauchy::Interpolate(mandelstam_list.interval_disp,dispersion_integral_above,gsl::InterpolationMethod::cubic);
    spline_disp_below= cauchy::Interpolate(mandelstam_list.interval_disp,dispersion_integral_below,gsl::InterpolationMethod::cubic);
    spline_disp_egg= cauchy::Interpolate(mandelstam_list.interval_phi,dispersion_integral_egg,gsl::InterpolationMethod::cubic);
}

void DispersionEtapEtaPiPi::build_splines(){
    spline_disp_above= cauchy::Interpolate(mandelstam_list.interval_disp,dispersion_integral_above,gsl::InterpolationMethod::cubic);
    spline_disp_below= cauchy::Interpolate(mandelstam_list.interval_disp,dispersion_integral_below,gsl::InterpolationMethod::cubic);

    spline_disp_t_channel_above= cauchy::Interpolate(mandelstam_list.interval_t_disp,dispersion_integral_t_channel_above,gsl::InterpolationMethod::cubic);
    spline_disp_t_channel_below= cauchy::Interpolate(mandelstam_list.interval_t_disp,dispersion_integral_t_channel_below,gsl::InterpolationMethod::cubic);
    spline_disp_minus_upper_a= cauchy::Interpolate(mandelstam_list.interval_y,dispersion_integral_minus_upper_a,gsl::InterpolationMethod::cubic);
    spline_disp_minus_upper_b= cauchy::Interpolate(mandelstam_list.interval_y,dispersion_integral_minus_upper_b,gsl::InterpolationMethod::cubic);
    spline_disp_minus_lower_a= cauchy::Interpolate(mandelstam_list.interval_y,dispersion_integral_minus_lower_a,gsl::InterpolationMethod::cubic);
    spline_disp_minus_lower_b= cauchy::Interpolate(mandelstam_list.interval_y,dispersion_integral_minus_lower_b,gsl::InterpolationMethod::cubic);
    spline_disp_plus_upper_a = cauchy::Interpolate(mandelstam_list.interval_y,dispersion_integral_plus_upper_a ,gsl::InterpolationMethod::cubic);
    spline_disp_plus_upper_b = cauchy::Interpolate(mandelstam_list.interval_y,dispersion_integral_plus_upper_b ,gsl::InterpolationMethod::cubic);
    spline_disp_plus_lower_a = cauchy::Interpolate(mandelstam_list.interval_y,dispersion_integral_plus_lower_a ,gsl::InterpolationMethod::cubic);
    spline_disp_plus_lower_b = cauchy::Interpolate(mandelstam_list.interval_y,dispersion_integral_plus_lower_b ,gsl::InterpolationMethod::cubic);
    spline_disp_zero_upper_a = cauchy::Interpolate(mandelstam_list.interval_y,dispersion_integral_zero_upper_a ,gsl::InterpolationMethod::cubic);
    spline_disp_zero_upper_b = cauchy::Interpolate(mandelstam_list.interval_y,dispersion_integral_zero_upper_b ,gsl::InterpolationMethod::cubic);
    spline_disp_zero_lower_a = cauchy::Interpolate(mandelstam_list.interval_y,dispersion_integral_zero_lower_a ,gsl::InterpolationMethod::cubic);
    spline_disp_zero_lower_b = cauchy::Interpolate(mandelstam_list.interval_y,dispersion_integral_zero_lower_b ,gsl::InterpolationMethod::cubic);
}

void DispersionV3Pi::build_splines(){
    spline_disp_above= cauchy::Interpolate(mandelstam_list.interval_disp,dispersion_integral_above,gsl::InterpolationMethod::cubic);
    spline_disp_below= cauchy::Interpolate(mandelstam_list.interval_disp,dispersion_integral_below,gsl::InterpolationMethod::cubic);
    spline_disp_egg= cauchy::Interpolate(mandelstam_list.interval_phi,dispersion_integral_egg,gsl::InterpolationMethod::cubic);
}

void DispersionX3Pi::build_splines(){
    spline_disp_above= cauchy::Interpolate(mandelstam_list.interval_disp,dispersion_integral_above,gsl::InterpolationMethod::cubic);
    spline_disp_below= cauchy::Interpolate(mandelstam_list.interval_disp,dispersion_integral_below,gsl::InterpolationMethod::cubic);
    spline_disp_egg= cauchy::Interpolate(mandelstam_list.interval_phi,dispersion_integral_egg,gsl::InterpolationMethod::cubic);
}

void DispersionT3Pi::build_splines(){
    spline_disp_above= cauchy::Interpolate(mandelstam_list.interval_disp,dispersion_integral_above,gsl::InterpolationMethod::cubic);
    spline_disp_below= cauchy::Interpolate(mandelstam_list.interval_disp,dispersion_integral_below,gsl::InterpolationMethod::cubic);
    spline_disp_egg= cauchy::Interpolate(mandelstam_list.interval_phi,dispersion_integral_egg,gsl::InterpolationMethod::cubic);
}

//------------------------------------------------------------------------
//-- Output --------------------------------------------------------------
//------------------------------------------------------------------------

// distinguish between all possible cases and return the interpolation of the desired function
Complex DispersionEta3Pi::operator()(double curve_parameter, Setting evaluation)const{
    double x=curve_parameter;
    switch (evaluation) {
        case Setting::above: return spline_disp_above(x);
        case Setting::below: return spline_disp_below(x);
        case Setting::egg:   return spline_disp_egg(x);  
        default:
            throw std::domain_error{"The demanded evaluation of the dispersion integral"
                        " is not needed for the chosen Process. Please check your input.\n"};
    }
}

Complex DispersionEtap3Pi::operator()(double curve_parameter, Setting evaluation)const{
    double x=curve_parameter;
    switch (evaluation) {
        case Setting::above: return spline_disp_above(x);
        case Setting::below: return spline_disp_below(x);
        case Setting::egg:   return spline_disp_egg(x);  
        default:
            throw std::domain_error{"The demanded evaluation of the dispersion integral"
                        " is not needed for the chosen Process. Please check your input.\n"};
    }
}

Complex DispersionEtapEtaPiPi::operator()(double curve_parameter, Setting evaluation)const{
    double x=curve_parameter;
    switch (evaluation) {
        case Setting::above:
            if(isospin==010 || 110){return spline_disp_t_channel_above(x);}
            else{return spline_disp_above(x);}
            break;
        case Setting::below:
            if(isospin==010 || 110){return spline_disp_t_channel_below(x);}
            else{return spline_disp_below(x);}
            break;
            // these evaluations are
        case Setting::minus_upper_a: return spline_disp_minus_upper_a(x);
        case Setting::minus_upper_b: return spline_disp_minus_upper_b(x);
        case Setting::minus_lower_a: return spline_disp_minus_lower_a(x);
        case Setting::minus_lower_b: return spline_disp_minus_lower_b(x);
        case Setting::plus_upper_a:  return spline_disp_plus_upper_a(x); 
        case Setting::plus_upper_b:  return spline_disp_plus_upper_b(x); 
        case Setting::plus_lower_a:  return spline_disp_plus_lower_a(x); 
        case Setting::plus_lower_b:  return spline_disp_plus_lower_b(x); 
        case Setting::zero_upper_a:  return spline_disp_zero_upper_a(x); 
        case Setting::zero_upper_b:  return spline_disp_zero_upper_b(x); 
        case Setting::zero_lower_a:  return spline_disp_zero_lower_a(x); 
        case Setting::zero_lower_b:  return spline_disp_zero_lower_b(x); 
        default:
            throw std::domain_error{"The demanded evaluation of the dispersion integral"
                " is not needed for the chosen Process. Please check your input.\n"};
    }
}

Complex DispersionV3Pi::operator()(double curve_parameter, Setting evaluation)const{
    double x=curve_parameter;
    switch (evaluation) {
        case Setting::above: return spline_disp_above(x);
        case Setting::below: return spline_disp_below(x);
        case Setting::egg:   return spline_disp_egg(x);  
        default:
            throw std::domain_error{"The demanded evaluation of the dispersion integral"
                        " is not needed for the chosen Process. Please check your input.\n"};
    }
}

Complex DispersionX3Pi::operator()(double curve_parameter, Setting evaluation)const{
    double x=curve_parameter;
    switch (evaluation) {
        case Setting::above: return spline_disp_above(x);
        case Setting::below: return spline_disp_below(x);
        case Setting::egg:   return spline_disp_egg(x);  
        default:
            throw std::domain_error{"The demanded evaluation of the dispersion integral"
                        " is not needed for the chosen Process. Please check your input.\n"};
    }
}

Complex DispersionT3Pi::operator()(double curve_parameter, Setting evaluation)const{
    double x=curve_parameter;
    switch (evaluation) {
        case Setting::above: return spline_disp_above(x);
        case Setting::below: return spline_disp_below(x);
        case Setting::egg:   return spline_disp_egg(x);  
        default:
            throw std::domain_error{"The demanded evaluation of the dispersion integral"
                        " is not needed for the chosen Process. Please check your input.\n"};
    }
}
//------------------------------------------------------------------------
//------------------------------------------------------------------------
}
