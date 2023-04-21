#ifndef DISPERSION_INTEGRAL_H
#define DISPERSION_INTEGRAL_H


#include "angular_average.h"
#include "type_aliases.h"
#include "facilities.h"
#include "gsl_interface.h"
#include "enums.h"
#include "phase.h"

#include <complex>
#include <cmath>
#include <vector>

/// Calculate the dispersion integrals for a given (fixed) subtraction scheme.
namespace disp{

using type_aliases::Complex;
using type_aliases::CFunction;
using namespace std::complex_literals;
using gsl::GaussLegendre;
using array::Array;
using cauchy::complex_integration;
using namespace constants;
using namespace enums;


using constants::pi;

using FunctionSetOmnes=const std::function<Complex(double, int, Setting)>&;
using FunctionSet=const std::function<Complex(double, int, SubtractionConstant)>&;

/// Parent class defining helpful functions for the numerator and creating lists of those.
class Numerator{
public:
    Numerator(double mass_1, double mass_2, double mass_3, double mass_4,
                       FunctionSetOmnes omnes,
                       FunctionSet matched_tilde,
                       int isospin, double epsilon, double validity, double cutoff,
                       SubtractionConstant sub, int max_subs);
    ///< @param mass_1, mass_2, mass_3 masses of the decay products
    ///< @param mass_4 decay mass
    ///<@param omnes Omnès functions, which can be switched to an evaluation above the cut, below the cut or along the egg-like contour.
    /// Note that the imported FunctionSet depends on the isospin I.
    ///<@param matched_tilde angular average multiplied with kinematic factor 'nu(s)' such that no singularities are contained.
    /// Note that the importet FunctionSet depends on the isospin I.
    ///<@param epsilon, validity matching point and range within which an expansion within the dispersion integral is used.
    ///<@param cutoff numeric cutoff of the dispersion integral
    ///<@param max_subs maximum number of subtraction constants (not valid for all of the processes)
    
    
    ///@brief Function which appears in the numerator of the dispersion integral.
    ///
    ///Note that here the power of s explicitly appears. If the subtraction scheme is changed this needs to be adjusted.
    ///For the eta processes this is fixed, while for the other processes it is controlled via max_subs.
    virtual Complex numerator(double s)const=0;
    
    //thresholds
    double s_I, s_III, s_0;
    double t_I, t_III, t_0;
    
    int isospin;
    int max_subs;
    
    ///@brief Creates lists for the real and imaginary part of the function numerator.
    ///
    ///These are only needed for the matching process, therefore it is totally sufficient to evaluate them on a smaller array of mandelstam s.
    virtual void build_numerator_real_imag_parts()=0;
    /// Now match the numerator and put real and imaginary parts together.
    virtual void build_integrand_expansions()=0;
    ///Splines the real and imaginary part of the lists from build_numerator_real_imag_parts.
    virtual void build_splines()=0;
    ///Saves the first (and second) derivative of the numerator, calculated in the match::Matching class.
    virtual void build_numerator_deriv();
    
    FunctionSet matched_tilde;
    FunctionSetOmnes omnes;
    SubtractionConstant sub;

    match::Matching matched;
    match::Matching matched_t_channel;

    /// Array to save the real part of the numerator function.
    std::vector<double> real_numerator_s;
    /// Array to save the imaginary part of the numerator function.
    std::vector<double> imag_numerator_s;

    gsl::Interpolate real_num_s;
    gsl::Interpolate imag_num_s;
    
    /// Matched expansions of the function numerator close to pseudo-thresholds.
    std::vector<Complex> integrand_expansion_list_s;
    ///@brief First derivative of the numerator above the pseudo-threshold.
    ///
    /// For the P-wave one needs the first derivative of numerator evaluated
    /// at pseudo-threshold. This derivative is constant for 's' below and above
    /// pseudo-threshold and is automatically determined by the matching Process.
    Complex numerator_1_derivative_s_above;
    ///@brief First derivative of the numerator below the pseudo-threshold.
    ///
    /// For the P-wave one needs the first derivative of numerator evaluated
    /// at pseudo-threshold. This derivative is constant for 's' below and above
    /// pseudo-threshold and is automatically determined by the matching Process.
    Complex numerator_1_derivative_s_below;

    ///@brief Second derivative of the numerator above the pseudo-threshold.
    ///
    /// For the D-wave one additionally needs the second derivative of numerator evaluated
    /// at pseudo-threshold. This derivative is constant for 's' below and above
    /// pseudo-threshold and is automatically determined by the matching Process.
    Complex numerator_2_derivative_s_above;
    ///@brief Second derivative of the numerator below the pseudo-threshold.
    ///
    /// For the D-wave one additionally needs the second derivative of numerator evaluated
    /// at pseudo-threshold. This derivative is constant for 's' below and above
    /// pseudo-threshold and is automatically determined by the matching Process.
    Complex numerator_2_derivative_s_below;
};

class NumeratorEta3Pi: public Numerator{
public:
    NumeratorEta3Pi(double mass_1, double mass_4,
                       FunctionSetOmnes omnes,
                       FunctionSet matched_tilde, phase::PhaseEta3Pi phases,
                       int isospin, double epsilon, double validity, double cutoff,
                       SubtractionConstant sub, int max_subs);    
    
    Complex numerator(double s)const override;
    
private:
    phase::PhaseEta3Pi phases;
    array::ArrayEta3Pi mandelstam_list;

    void build_numerator_real_imag_parts() override;
    void build_integrand_expansions() override;
    void build_splines() override;
};

class NumeratorEtap3Pi: public Numerator{
public:
    NumeratorEtap3Pi(double mass_1, double mass_4,
                       FunctionSetOmnes omnes,
                       FunctionSet matched_tilde, phase::PhaseEtap3Pi phases,
                       int isospin, double epsilon, double validity, double cutoff,
                       SubtractionConstant sub, int max_subs);    
    
    Complex numerator(double s)const override;
    
private:
    phase::PhaseEtap3Pi phases;
    array::ArrayEtap3Pi mandelstam_list;

    void build_numerator_real_imag_parts() override;
    void build_integrand_expansions() override;
    void build_splines() override;
};

class NumeratorEtapEtaPiPi: public Numerator{
public:
    NumeratorEtapEtaPiPi(double mass_1, double mass_3, double mass_4,
                       FunctionSetOmnes omnes,
                       FunctionSet matched_tilde, phase::PhaseEtapEtaPiPi phases,
                       int isospin, double epsilon, double validity, double cutoff,
                       SubtractionConstant sub, int max_subs);    
    
    Complex numerator(double s)const override;

    std::vector<Complex> integrand_expansion_list_t;
private:
    phase::PhaseEtapEtaPiPi phases;
    array::ArrayEtapEtaPiPi mandelstam_list;

    /// Array to save the real part of the t-channel numerator function.
    std::vector<double> real_numerator_t;
    /// Array to save the imaginary part of the numerator function.
    std::vector<double> imag_numerator_t;

    gsl::Interpolate real_num_t;
    gsl::Interpolate imag_num_t;
    
    void build_numerator_real_imag_parts() override;
    void build_integrand_expansions() override;
    void build_splines() override;
};

class NumeratorV3Pi: public Numerator{
public:
    NumeratorV3Pi(double mass_1, double mass_4,
                       FunctionSetOmnes omnes,
                       FunctionSet matched_tilde, phase::PhaseV3Pi phases,
                       int isospin, double epsilon, double validity, double cutoff,
                       SubtractionConstant sub, int max_subs);    
    
    Complex numerator(double s)const override;
    
private:
    phase::PhaseV3Pi phases;
    array::ArrayV3Pi mandelstam_list;

    void build_numerator_real_imag_parts() override;
    void build_integrand_expansions() override;
    void build_splines() override;
};

class NumeratorX3Pi: public Numerator{
public:
    NumeratorX3Pi(double mass_1, double mass_4,
                       FunctionSetOmnes omnes,
                       FunctionSet matched_tilde, phase::PhaseX3Pi phases,
                       int isospin, double epsilon, double validity, double cutoff,
                       SubtractionConstant sub, int max_subs);    
    
    Complex numerator(double s)const override;
    
private:
    phase::PhaseX3Pi phases;
    array::ArrayX3Pi mandelstam_list;

    void build_numerator_real_imag_parts() override;
    void build_integrand_expansions() override;
    void build_splines() override;
};

class NumeratorT3Pi: public Numerator{
public:
    NumeratorT3Pi(double mass_1, double mass_4,
                       FunctionSetOmnes omnes,
                       FunctionSet matched_tilde, phase::PhaseT3Pi phases,
                       int isospin, double epsilon, double validity, double cutoff,
                       SubtractionConstant sub, int max_subs);    
    
    Complex numerator(double s)const override;
    
private:
    phase::PhaseT3Pi phases;
    array::ArrayT3Pi mandelstam_list;

    void build_numerator_real_imag_parts() override;
    void build_integrand_expansions() override;
    void build_splines() override;
    void build_numerator_deriv() override;
};



// Interpolate the matched expansions to decrease computation time
// and calculate the dispersion integrals. The results for the latter are
// stored in lists.
/// Partent class to calculate the dispersion integrals. 
class Dispersion{
public:
    Dispersion(double mass_1, double mass_2, double mass_3, double mass_4,
                    FunctionSetOmnes omnes,
                    FunctionSet matched_tilde,
                    int isospin, double epsilon, double validity, double cutoff,
                    SubtractionConstant sub, int max_subs);
    ///< @param mass_1, mass_2, mass_3 masses of the decay products
    ///< @param mass_4 decay mass
    ///<@param omnes Omnès functions, which can be switched to an evaluation above the cut, below the cut or along the egg-like contour.
    /// Note that the imported FunctionSet depends on the isospin I.
    ///<@param matched_tilde angular average multiplied with kinematic factor 'nu(s)' such that no singularities are contained.
    /// Note that the importet FunctionSet depends on the isospin I.
    ///<@param epsilon, validity matching point and range within which an expansion within the dispersion integral is used.
    ///<@param cutoff numeric cutoff of the dispersion integral
    ///<@param max_subs maximum number of subtraction constants (not valid for all of the processes)
    
    ///Function that returns the splined result of the dispersion integral.
    virtual Complex operator()(double s, Setting evaluation)const=0; 

    double mass_1;
    int isospin;
    //thresholds
    double s_I, s_III, s_0;
    double t_I, t_III, t_0;
    double validity, cutoff;
    int max_subs;
    double infinitesimal=0.0000001;
    
    /// List for the dispersion integral evaluated infinitesimally above the cut.
    std::vector<Complex> dispersion_integral_above={};
    /// List for the dispersion integral evaluated infinitesimally below the cut.
    std::vector<Complex> dispersion_integral_below={};
 
    cauchy::Interpolate spline_disp_above;
    cauchy::Interpolate spline_disp_below;

    ///Function to import the matched numerator from the Numerator class.
    virtual Complex numerator(double s)const=0;
    
    /// Analytically calculated R integral of order 1/2 contributing to the dispersion integral.
    Complex analytic_integral_R12(Complex s, double min, double max, Setting evaluation)const;
    ///<@param s Mandelstam s
    ///<@param min lower limit of integration
    ///<@param max upper limit of integration
    ///<@param evaluation choose whether to evaluate "above" cut, "below" cut or along the "egg"-like contour

    ///@brief Analytically calculated R integral of order (2n+1)/2 contributing to the dispersion integral.
    ///
    ///Calculated recursively via R12.
    Complex analytic_integral_R(Complex s, double min, double max, Setting evaluation, int n)const;
    ///<@param s Mandelstam s
    ///<@param min lower limit of integration
    ///<@param max upper limit of integration
    ///<@param evaluation choose whether to evaluate "above" cut, "below" cut or along the "egg"-like contour
    ///<@param n is the order of the R integral is determined by (2n+1)/2
    
    /// Analytically calculated Q integral of order 1/2 contributing to the dispersion integral.
    Complex analytic_integral_Q12(Complex s, double min, double max)const;
    ///<@param s Mandelstam s
    ///<@param min lower limit of integration
    ///<@param max upper limit of integration

    ///@brief Analytically calculated Q integral of order (2n+1)/2 contributing to the dispersion integral.
    ///
    ///Calculated recursively via Q12.
    Complex analytic_integral_Q(Complex s, double min, double max, int n)const;
    ///<@param s Mandelstam s
    ///<@param min lower limit of integration
    ///<@param max upper limit of integration
    ///<@param n is the order of the Q integral is determined by (2n+1)/2

    ///Numerically calculating the integrals with the Cauchy singularity.
    virtual Complex numeric_integral_cauchy(double s, double lower_limit,
                                      double upper_limit, bool adaptive)const=0;
    ///<@param s Mandelstam s
    ///<@param min lower limit of integration
    ///<@param max upper limit of integration
    ///<@param precision use "on" to evaluate the integrals with a very high precision (adaptive routines)
    /// use "off" to decrease the computation time a lot
    
    
    ///Numerically calculating the integrals with the pseudo-threshold singularity.
    virtual Complex numeric_integral_pseudo(Complex s, double lower_limit,
                               double upper_limit, bool adaptive)const=0;
    ///<@param s Mandelstam s
    ///<@param min lower limit of integration
    ///<@param max upper limit of integration
    ///<@param precision use "on" to evaluate the integrals with a very high precision (adaptive routines)
    /// use "off" to decrease the computation time a lot
    
    
    /// Calculated the dispersion integrals. Used to create lists.
    virtual Complex dispersion_integral(Complex s, Setting evaluation)const=0;

    virtual void build_dispersion_integral_lists()=0;
    virtual void build_splines()=0;
    
    /// Integrand to evaluate the integrals with the singularity at pseudo-threshold.
    virtual Complex integrand_pseudo(Complex s, double s_prime)const=0;
    ///<@param s complex-valued Mandelstam s
    ///<@param s_prime real-valued integration variable  
    
    /// Integrand to evaluate the integrals with the Cauchy singularity.
    virtual Complex integrand_cauchy(double s, double s_prime)const=0;
    ///<@param s real-valued Mandelstam s
    ///<@param s_prime real-valued integration variable
        
    /// Function that returns the result of the dispersion integral, where real-valued s is below the pseudo-threshold.    
    virtual Complex disp_integral_below_pseudo(Complex s, bool adaptive, Setting evaluation)const=0;
    /// Function that returns the result of the dispersion integral, where real-valued s is above the pseudo-threshold.  
    virtual Complex disp_integral_above_pseudo(Complex s, bool adaptive, Setting evaluation)const=0;
    /// Function that returns the result of the dispersion integral, where s is below two-pion threshold or in the complex plane.
    virtual Complex disp_integral_trivial(Complex s, bool adaptive)const=0;
};

class DispersionEta3Pi: public Dispersion{
public:
    DispersionEta3Pi(double mass_1, double mass_4,
                    FunctionSetOmnes omnes,
                    FunctionSet matched_tilde,
                    phase::PhaseEta3Pi phase,
                    int isospin, double epsilon, double validity, double cutoff,
                    SubtractionConstant sub, int max_subs);
    
    Complex operator()(double s, Setting evaluation)const override; 
private:
    phase::PhaseEta3Pi phases;
    path::PolarEgg polar_egg;
    array::ArrayEta3Pi mandelstam_list;
    NumeratorEta3Pi numeratorc;
    
    std::vector<Complex> dispersion_integral_egg={};

    cauchy::Interpolate spline_disp_egg;
        
    Complex numerator(double s)const override;
    
    Complex numeric_integral_cauchy(double s, double lower_limit,
                                      double upper_limit, bool adaptive)const override;
    
    Complex numeric_integral_pseudo(Complex s, double lower_limit,
                               double upper_limit, bool adaptive)const override;
    
    Complex dispersion_integral(Complex s, Setting evaluation)const override;

    Complex numerator_1_derivative_s_above;
    Complex numerator_1_derivative_s_below;
    
    void build_dispersion_integral_lists() override;
    void build_splines() override;
    cauchy::Interpolate integrand_expansion_s;
    
    Complex integrand_pseudo(Complex s, double s_prime)const override;

    Complex integrand_cauchy(double s, double s_prime)const override;
    
    double cusp=4.0*pow(mass_kaon(),2.0)/pow(mass_pi(),2.0)*pow(mass_1,2.0);

    
    Complex disp_integral_below_pseudo(Complex s, bool adaptive, Setting evaluation)const override;
    Complex disp_integral_above_pseudo(Complex s, bool adaptive, Setting evaluation)const override;
    Complex disp_integral_trivial(Complex s, bool adaptive)const override;
};

class DispersionEtap3Pi: public Dispersion{
public:
    DispersionEtap3Pi(double mass_1, double mass_4,
                    FunctionSetOmnes omnes,
                    FunctionSet matched_tilde,
                    phase::PhaseEtap3Pi phase,
                    int isospin, double epsilon, double validity, double cutoff,
                    SubtractionConstant sub, int max_subs);
    
    Complex operator()(double s, Setting evaluation)const override; 
private:
    phase::PhaseEtap3Pi phases;
    path::PolarEgg polar_egg;
    array::ArrayEtap3Pi mandelstam_list;
    NumeratorEtap3Pi numeratorc;

    std::vector<Complex> dispersion_integral_egg={};

    cauchy::Interpolate spline_disp_egg;
        
    Complex numerator(double s)const override;
    
    
    Complex numeric_integral_cauchy(double s, double lower_limit,
                                      double upper_limit, bool adaptive)const override;
    
    Complex numeric_integral_pseudo(Complex s, double lower_limit,
                               double upper_limit, bool adaptive)const override;
    
    Complex dispersion_integral(Complex s, Setting evaluation)const override;

    Complex numerator_1_derivative_s_above;
    Complex numerator_1_derivative_s_below;
    
    void build_dispersion_integral_lists() override;
    void build_splines() override;
    cauchy::Interpolate integrand_expansion_s;
    
    Complex integrand_pseudo(Complex s, double s_prime)const override;

    Complex integrand_cauchy(double s, double s_prime)const override;
    
    double cusp=4.0*pow(mass_kaon(),2.0)/pow(mass_pi(),2.0)*pow(mass_1,2.0);

    
    Complex disp_integral_below_pseudo(Complex s, bool adaptive, Setting evaluation)const override;
    Complex disp_integral_above_pseudo(Complex s, bool adaptive, Setting evaluation)const override;
    Complex disp_integral_trivial(Complex s, bool adaptive)const override;
};

class DispersionEtapEtaPiPi: public Dispersion{
public:
    DispersionEtapEtaPiPi(double mass_1, double mass_3, double mass_4,
                    FunctionSetOmnes omnes,
                    FunctionSet matched_tilde,
                    phase::PhaseEtapEtaPiPi phase,
                    int isospin, double epsilon, double validity, double cutoff,
                    SubtractionConstant sub, int max_subs);
    
    Complex operator()(double s, Setting evaluation)const override;

private:
    phase::PhaseEtapEtaPiPi phases;
    path_eta_pi_pi::Path gamma;
    array::ArrayEtapEtaPiPi mandelstam_list;
    NumeratorEtapEtaPiPi numeratorc;

    std::vector<Complex> dispersion_integral_t_channel_above={};
    std::vector<Complex> dispersion_integral_t_channel_below={};
    
    // Evaluation along the curve parameters for integration contour in eta'-> eta pi pi
    std::vector<Complex> dispersion_integral_minus_upper_a;
    std::vector<Complex> dispersion_integral_minus_upper_b;
    std::vector<Complex> dispersion_integral_minus_lower_a;
    std::vector<Complex> dispersion_integral_minus_lower_b;
    std::vector<Complex> dispersion_integral_plus_upper_a;
    std::vector<Complex> dispersion_integral_plus_upper_b;
    std::vector<Complex> dispersion_integral_plus_lower_a;
    std::vector<Complex> dispersion_integral_plus_lower_b;
    std::vector<Complex> dispersion_integral_zero_upper_a;
    std::vector<Complex> dispersion_integral_zero_upper_b;
    std::vector<Complex> dispersion_integral_zero_lower_a;
    std::vector<Complex> dispersion_integral_zero_lower_b;

    // only for eta'-> eta pi pi
    cauchy::Interpolate spline_disp_t_channel_above;
    cauchy::Interpolate spline_disp_t_channel_below;
    cauchy::Interpolate spline_disp_minus_upper_a;
    cauchy::Interpolate spline_disp_minus_upper_b;
    cauchy::Interpolate spline_disp_minus_lower_a;
    cauchy::Interpolate spline_disp_minus_lower_b;
    cauchy::Interpolate spline_disp_plus_upper_a;
    cauchy::Interpolate spline_disp_plus_upper_b;
    cauchy::Interpolate spline_disp_plus_lower_a;
    cauchy::Interpolate spline_disp_plus_lower_b;
    cauchy::Interpolate spline_disp_zero_upper_a;
    cauchy::Interpolate spline_disp_zero_upper_b;
    cauchy::Interpolate spline_disp_zero_lower_a;
    cauchy::Interpolate spline_disp_zero_lower_b;    
    
    Complex numerator(double s)const override; 
    
    Complex numeric_integral_cauchy(double s, double lower_limit,
                                      double upper_limit, bool adaptive)const override;
     
    Complex numeric_integral_pseudo(Complex s, double lower_limit,
                               double upper_limit, bool adaptive)const override;
    
    Complex dispersion_integral(Complex s, Setting evaluation)const override;

    Complex numerator_1_derivative_s_above;
    Complex numerator_1_derivative_s_below;
    
    void build_dispersion_integral_lists() override;
    void build_splines() override;
    cauchy::Interpolate integrand_expansion_s;
    cauchy::Interpolate integrand_expansion_t;
    
    Complex integrand_pseudo(Complex s, double s_prime)const override;
    
    Complex integrand_cauchy(double s, double s_prime)const override;
    
    double cusp=4.0*pow(mass_kaon(),2.0)/pow(mass_pi(),2.0)*pow(mass_1,2.0);
   
    Complex disp_integral_below_pseudo(Complex s, bool adaptive, Setting evaluation)const override;
    Complex disp_integral_above_pseudo(Complex s, bool adaptive, Setting evaluation)const override;
    Complex disp_integral_trivial(Complex s, bool adaptive)const override;
};

class DispersionV3Pi: public Dispersion{
public:
    DispersionV3Pi(double mass_1, double mass_4,
                    FunctionSetOmnes omnes,
                    FunctionSet matched_tilde,
                    phase::PhaseV3Pi phase,
                    int isospin, double epsilon, double validity, double cutoff,
                    SubtractionConstant sub, int max_subs);
    
    Complex operator()(double s, Setting evaluation)const override; 

    phase::PhaseV3Pi phases;
    path::PolarEgg polar_egg;
    array::ArrayV3Pi mandelstam_list;
    NumeratorV3Pi numeratorc;
    
    std::vector<Complex> dispersion_integral_egg={};

    cauchy::Interpolate spline_disp_egg;
        
    Complex numerator(double s)const override;
    
    Complex numeric_integral_cauchy(double s, double lower_limit,
                                      double upper_limit, bool adaptive)const override;
    
    Complex numeric_integral_pseudo(Complex s, double lower_limit,
                               double upper_limit, bool adaptive)const override;
    
    Complex dispersion_integral(Complex s, Setting evaluation)const override;

    Complex numerator_1_derivative_s_above;
    Complex numerator_1_derivative_s_below;
    
    void build_dispersion_integral_lists() override;
    void build_splines() override;
    cauchy::Interpolate integrand_expansion_s;
    
    Complex integrand_pseudo(Complex s, double s_prime)const override;

    Complex integrand_cauchy(double s, double s_prime)const override;
    
    Complex disp_integral_below_pseudo(Complex s, bool adaptive, Setting evaluation)const override;
    Complex disp_integral_above_pseudo(Complex s, bool adaptive, Setting evaluation)const override;
    Complex disp_integral_trivial(Complex s, bool adaptive)const override;
};

class DispersionX3Pi: public Dispersion{
public:
    DispersionX3Pi(double mass_1, double mass_4,
                    FunctionSetOmnes omnes,
                    FunctionSet matched_tilde,
                    phase::PhaseX3Pi phase,
                    int isospin, double epsilon, double validity, double cutoff,
                    SubtractionConstant sub, int max_subs);
    
    Complex operator()(double s, Setting evaluation)const override; 

    phase::PhaseX3Pi phases;
    path::PolarEgg polar_egg;
    array::ArrayX3Pi mandelstam_list;
    NumeratorX3Pi numeratorc;

    std::vector<Complex> dispersion_integral_egg={};

    cauchy::Interpolate spline_disp_egg;
        
    Complex numerator(double s)const override;
    
    Complex numeric_integral_cauchy(double s, double lower_limit,
                                      double upper_limit, bool adaptive)const override;
    
    Complex numeric_integral_pseudo(Complex s, double lower_limit,
                               double upper_limit, bool adaptive)const override;
    
    Complex dispersion_integral(Complex s, Setting evaluation)const override;

    Complex numerator_1_derivative_s_above;
    Complex numerator_1_derivative_s_below;
    
    void build_dispersion_integral_lists() override;
    void build_splines() override;
    cauchy::Interpolate integrand_expansion_s;
    
    Complex integrand_pseudo(Complex s, double s_prime)const override;

    Complex integrand_cauchy(double s, double s_prime)const override;
    
    Complex disp_integral_below_pseudo(Complex s, bool adaptive, Setting evaluation)const override;
    Complex disp_integral_above_pseudo(Complex s, bool adaptive, Setting evaluation)const override;
    Complex disp_integral_trivial(Complex s, bool adaptive)const override;
};

class DispersionT3Pi: public Dispersion{
public:
    DispersionT3Pi(double mass_1, double mass_4,
                    FunctionSetOmnes omnes,
                    FunctionSet matched_tilde,
                    phase::PhaseT3Pi phase,
                    int isospin, double epsilon, double validity, double cutoff,
                    SubtractionConstant sub, int max_subs);
    
    Complex operator()(double s, Setting evaluation)const override; 

    phase::PhaseT3Pi phases;
    path::PolarEgg polar_egg;
    array::ArrayT3Pi mandelstam_list;
    NumeratorT3Pi numeratorc;
    
    std::vector<Complex> dispersion_integral_egg={};

    cauchy::Interpolate spline_disp_egg;
        
    Complex numerator(double s)const override;
    
    Complex numeric_integral_cauchy(double s, double lower_limit,
                                      double upper_limit, bool adaptive)const override;
    
    Complex numeric_integral_pseudo(Complex s, double lower_limit,
                               double upper_limit, bool adaptive)const override;

    Complex dispersion_integral(Complex s, Setting evaluation)const override;

    Complex numerator_1_derivative_s_above;
    Complex numerator_1_derivative_s_below;
    Complex numerator_2_derivative_s_above;
    Complex numerator_2_derivative_s_below;
    
    void build_dispersion_integral_lists() override;
    void build_splines() override;
    cauchy::Interpolate integrand_expansion_s;
    
    Complex integrand_pseudo(Complex s, double s_prime)const override;

    Complex integrand_cauchy(double s, double s_prime)const override;
    
    Complex disp_integral_below_pseudo(Complex s, bool adaptive, Setting evaluation)const override;
    Complex disp_integral_above_pseudo(Complex s, bool adaptive, Setting evaluation)const override;
    Complex disp_integral_trivial(Complex s, bool adaptive)const override;
};

} // disp

#endif // DISPERSION_INTEGRAL_H
