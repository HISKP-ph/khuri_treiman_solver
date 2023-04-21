// Evaluate the angular averages along a complex integration contour
// that avoids the cut.
// Additionally match an expansion of the angular average close to scattering thresholds.



#ifndef ANGULAR_AVERAGE_H
#define ANGULAR_AVERAGE_H

#include "path.h"
#include "path_eta_pi_pi.h"
#include "matching.h"
#include "array.h"
#include "grid.h"
#include "cauchy.h"
#include "helpers.h"
#include "type_aliases.h"
#include "facilities.h"
#include "gsl_interface.h"
#include "splined_path.h"
#include "enums.h"

#include <complex>
#include <cmath>
#include <vector>
#include <string>


namespace angular_average{

using enums::Setting;
using enums::SubtractionConstant;
using gsl::GaussLegendre;
using match::Matching;
using constants::pi;
using type_aliases::Complex;
using type_aliases::CFunction;
using array::Array;
using namespace std::complex_literals;
using type_aliases::Complex;
using type_aliases::CFunction;
using type_aliases::Curve;
using cauchy::complex_integration;
using namespace constants;
using namespace enums;


using Function=const std::function<Complex(double)>&;
using FunctionSet=const std::function<Complex(double, int, SubtractionConstant, Setting)>&;
// If one sets the string to either "above" or "below" one evaluates the function for Mandelstam s (in units of squared pion mass)
// infinitesimally above or below the cut. If the string is chosen to be "egg", the corresponding function will be a function in terms of
// a real-valued curve-paramter (i.e. an angle in radians) and not in terms of 's' anymore.

///Parent Class to calculate the angular averages
class AngularAverage{
public:
    AngularAverage(double mass_1, double mass_2, double mass_3, double mass_4,
                        FunctionSet amplitude, int isospin,
                        double cutoff, int iteration, double epsilon, double validity, SubtractionConstant sub);
        ///<@param mass_1, mass_2, mass_3 masses of the decay products
        ///<@param mass_4 decay mass
        ///<@param amplitude basis functions
        ///<@param cutoff numerical cutoff for the dispersion integral
        ///<@param iteration if iteration=0, use the fact that the Omnès function is perfectly Schwartz
        ///<@param epsilon, validity matching point and range for the matching procedure of the angular average
        ///<@param sub specify the subtraction constant whose basis functions shall be evaluated (c.f. "enums.h")

    ///Function that returns the splined matched angular averages.
    virtual Complex operator()(double s)const=0;

    double mass_1, mass_2, mass_3, mass_4;
    // thresholds
    double s_I, s_II, s_III, s_IV, s_0;
    double t_I, t_II, t_III, t_IV, t_0;
    double u_II;
    double s_x, t_x;
    double validity;
    double epsilon;
 
    int isospin;
    int iteration;
    SubtractionConstant sub;
    FunctionSet amplitude;
       
    ///Källén function, the last to argumetns are masses and not their squares!
    Complex kallen(double s, double mass1, double mass2)const;
    ///Function nu depending on kinematics
    Complex nu(double s)const;
    ///@brief Specific integrand depending on the reconstruction theorem of the process.
    ///
    ///Note that the order in kappa must be the same for all terms. Therefore some terms need to multiplied by some order of kappe explicitly.
    virtual Complex integrand(double s, double x, Setting evaluation, TypeOfAverage aver=TypeOfAverage::zero)const=0;
    ///Returns the splined angular averages without the matching.
    virtual Complex angular_average(double s)const=0;
    ///Calculates the integral via the function integral_3em or integral_2em.
    virtual Complex eval_int(double s, TypeOfAverage aver=TypeOfAverage::zero)const=0;
    ///Builds the angular averages via a list of mandelstam s from array::Array.
    virtual void build_tilde_functions()=0;
    ///Separates two list from the result of build_tilde_functions.
    virtual void build_matching_list()=0; 
    ///Build the matched list using the lists from build_matching_list and the matching procedure from match::Matching.
    virtual void build_matched_tilde()=0;
    ///Splines the tilde function from the lists of build_tilde_functions.
    virtual void make_splines()=0;
    ///Splines the matched tilde function from the lists of build_matched_tilde.
    virtual void make_splines_matched()=0;
    ///Calculates the angular average integral for a process with 3 equal masses in the final state.
    Complex integral_3em(double s, const splined_path::SplinedPath &polar_egg)const;
    ///Calculates the angular average integral for a process with 2 different masses in the final state.
    Complex integral_2em(double s, TypeOfAverage aver, const path_eta_pi_pi::Path &contour)const;
    ///Helper function for the integral_2em function.
    Complex helper(double x, TypeOfAverage aver ,Curve integrand_lower_a, Curve integrand_upper_a, Curve integrand_lower_b, Curve integrand_upper_b, const path_eta_pi_pi::Path &contour)const;

    match::Matching matched;

    ///Complex list to store the angular average.  
    std::vector<Complex> angular_average_list={};
    ///Complex list to store the matched angular average.
    std::vector<Complex> matched_angular_average;
    
    /// List for the real part of the angular averages. This is needed to match an expansion close to the thresholds.
    std::vector<double>real_angular_average_list;
    /// List for the imaginary part of the angular averages. This is needed to match an expansion close to the thresholds.
    std::vector<double>imag_angular_average_list;

    gsl::Interpolate real_angular_average;
    gsl::Interpolate imag_angular_average;
    
    cauchy::Interpolate spline_angular_average;

    cauchy::Interpolate spline_matched_angular_average;
};


class AngularAverageEta3Pi: public AngularAverage{
public:
    AngularAverageEta3Pi(double mass_1, double mass_4,
                        FunctionSet amplitude, int isospin,
                        double cutoff, int iteration, double epsilon, double validity, SubtractionConstant sub);

    Complex operator()(double s) const override;

private:    
    Complex integrand(double s, double x, Setting evaluation, TypeOfAverage aver=TypeOfAverage::zero)const override;
    Complex angular_average(double s)const override;
    Complex eval_int(double s, TypeOfAverage aver=TypeOfAverage::zero)const override;

    array::ArrayEta3Pi mandelstam_list;

    splined_path::SplinedPathEta3Pi polar_egg;
    splined_path::SplinedPathDerivativeEta3Pi egg_deriv; 
    
    void build_tilde_functions() override;
    void build_matching_list() override; 
    void build_matched_tilde() override;
    void make_splines() override;
    void make_splines_matched() override;
};

class AngularAverageEta3PiC: public AngularAverage{
public:
    AngularAverageEta3PiC(double mass_1, double mass_4,
                        FunctionSet amplitude, int isospin,
                        double cutoff, int iteration, double epsilon, double validity, SubtractionConstant sub);

    Complex operator()(double s) const override;
    
private:    
    Complex integrand(double s, double x, Setting evaluation, TypeOfAverage aver=TypeOfAverage::zero)const override;
    Complex angular_average(double s)const override;
    Complex eval_int(double s, TypeOfAverage aver=TypeOfAverage::zero)const override;

    array::ArrayEta3Pi mandelstam_list;

    splined_path::SplinedPathEta3Pi polar_egg;
    splined_path::SplinedPathDerivativeEta3Pi egg_deriv; 
    
    void build_tilde_functions() override;
    void build_matching_list() override; 
    void build_matched_tilde() override;
    void make_splines() override;
    void make_splines_matched() override;
};

class AngularAverageEtap3Pi: public AngularAverage{
public:
    AngularAverageEtap3Pi(double mass_1, double mass_4,
                        FunctionSet amplitude, int isospin,
                        double cutoff, int iteration, double epsilon, double validity, SubtractionConstant sub);

    Complex operator()(double s) const override;
    
private:    
    Complex integrand(double s, double x, Setting evaluation, TypeOfAverage aver=TypeOfAverage::zero)const override;
    Complex angular_average(double s)const override;
    Complex eval_int(double s, TypeOfAverage aver=TypeOfAverage::zero)const override;

    array::ArrayEtap3Pi mandelstam_list;

    splined_path::SplinedPathEtap3Pi polar_egg;
    splined_path::SplinedPathDerivativeEtap3Pi egg_deriv; 
    
    void build_tilde_functions() override;
    void build_matching_list() override; 
    void build_matched_tilde() override;
    void make_splines() override;
    void make_splines_matched() override;
};

class AngularAverageEtap3PiC: public AngularAverage{
public:
    AngularAverageEtap3PiC(double mass_1, double mass_4,
                        FunctionSet amplitude, int isospin,
                        double cutoff, int iteration, double epsilon, double validity, SubtractionConstant sub);

    Complex operator()(double s) const override;
    
private:    
    Complex integrand(double s, double x, Setting evaluation, TypeOfAverage aver=TypeOfAverage::zero)const override;
    Complex angular_average(double s)const override;
    Complex eval_int(double s, TypeOfAverage aver=TypeOfAverage::zero)const override;

    array::ArrayEtap3Pi mandelstam_list;

    splined_path::SplinedPathEtap3Pi polar_egg;
    splined_path::SplinedPathDerivativeEtap3Pi egg_deriv; 
    
    void build_tilde_functions() override;
    void build_matching_list() override; 
    void build_matched_tilde() override;
    void make_splines() override;
    void make_splines_matched() override;
};

class AngularAverageEtapEtaPiPi: public AngularAverage{
public:
    AngularAverageEtapEtaPiPi(double mass_1, double mass_3, double mass_4,
                        FunctionSet amplitude, int isospin,
                        double cutoff, int iteration, double epsilon, double validity, SubtractionConstant sub);

    Complex operator()(double s) const override;
private:    
    Complex integrand(double s, double x, Setting evaluation, TypeOfAverage aver=TypeOfAverage::zero)const override;
    Complex angular_average(double s)const override;
    Complex eval_int(double s, TypeOfAverage aver=TypeOfAverage::zero)const override;

    array::ArrayEtapEtaPiPi mandelstam_list;
    path_eta_pi_pi::Path contour;
    
    void build_tilde_functions() override;
    void build_matching_list() override;
    void build_matched_tilde() override;
    void make_splines() override;
    void make_splines_matched() override;
    
    ///Complex list to store the t-channel angular average.
    std::vector<Complex> angular_average_intermediate_list={};
    ///Complex list to store the matched t-channel angular average.
    std::vector<Complex> matched_angular_average_intermediate;

    /// List for the real part of the t-channel angular averages. This is needed to match an expansion close to the thresholds.
    std::vector<double>real_angular_average_intermediate_list;
    /// List for the imaginary part of the t-channel angular averages. This is needed to match an expansion close to the thresholds.
    std::vector<double>imag_angular_average_intermediate_list;
    
    match::Matching matched_t_channel;

    Interpolate real_angular_average_intermediate;
    Interpolate imag_angular_average_intermediate;
    cauchy::Interpolate spline_angular_average_intermediate;
    cauchy::Interpolate spline_matched_angular_average_intermediate;
};

class AngularAverageEtapEtaPiPiC: public AngularAverage{
public:
    AngularAverageEtapEtaPiPiC(double mass_1, double mass_3, double mass_4,
                        FunctionSet amplitude, int isospin,
                        double cutoff, int iteration, double epsilon, double validity, SubtractionConstant sub);

    Complex operator()(double s) const override;
private:    
    Complex integrand(double s, double x, Setting evaluation, TypeOfAverage aver=TypeOfAverage::zero)const override;
    Complex angular_average(double s)const override;
    Complex eval_int(double s, TypeOfAverage aver=TypeOfAverage::zero)const override;

    array::ArrayEtapEtaPiPi mandelstam_list;
    path_eta_pi_pi::Path contour;
    
    void build_tilde_functions() override;
    void build_matching_list() override;
    void build_matched_tilde() override;
    void make_splines() override;
    void make_splines_matched() override;

    std::vector<Complex> angular_average_intermediate_list={};
    std::vector<Complex> matched_angular_average_intermediate;

    std::vector<double>real_angular_average_intermediate_list;
    std::vector<double>imag_angular_average_intermediate_list;
    
    match::Matching matched_t_channel;

    Interpolate real_angular_average_intermediate;
    Interpolate imag_angular_average_intermediate;
    cauchy::Interpolate spline_angular_average_intermediate;
    cauchy::Interpolate spline_matched_angular_average_intermediate;
};

class AngularAverageV3Pi: public AngularAverage{
public:
    AngularAverageV3Pi(double mass_1, double mass_4,
                        FunctionSet amplitude, int isospin,
                        double cutoff, int iteration, double epsilon, double validity, SubtractionConstant sub);

    Complex operator()(double s) const override;

private:    
    Complex integrand(double s, double x, Setting evaluation, TypeOfAverage aver=TypeOfAverage::zero)const override;
    Complex angular_average(double s)const override;
    Complex eval_int(double s, TypeOfAverage aver=TypeOfAverage::zero)const override;

    array::ArrayV3Pi mandelstam_list;

    splined_path::SplinedPathV3Pi polar_egg;
    splined_path::SplinedPathDerivativeV3Pi egg_deriv; 
    
    void build_tilde_functions() override;
    void build_matching_list() override; 
    void build_matched_tilde() override;
    void make_splines() override;
    void make_splines_matched() override;
};

class AngularAverageX3Pi: public AngularAverage{
public:
    AngularAverageX3Pi(double mass_1, double mass_4,
                        FunctionSet amplitude, int isospin,
                        double cutoff, int iteration, double epsilon, double validity, SubtractionConstant sub);

    Complex operator()(double s) const override;

private:    
    Complex integrand(double s, double x, Setting evaluation, TypeOfAverage aver=TypeOfAverage::zero)const override;
    Complex angular_average(double s)const override;
    Complex eval_int(double s, TypeOfAverage aver=TypeOfAverage::zero)const override;

    array::ArrayX3Pi mandelstam_list;

    splined_path::SplinedPathX3Pi polar_egg;
    splined_path::SplinedPathDerivativeX3Pi egg_deriv; 
    
    void build_tilde_functions() override;
    void build_matching_list() override; 
    void build_matched_tilde() override;
    void make_splines() override;
    void make_splines_matched() override;
};

class AngularAverageT3Pi: public AngularAverage{
public:
    AngularAverageT3Pi(double mass_1, double mass_4,
                        FunctionSet amplitude, int isospin,
                        double cutoff, int iteration, double epsilon, double validity, SubtractionConstant sub);

    Complex operator()(double s) const override;

private:    
    Complex integrand(double s, double x, Setting evaluation, TypeOfAverage aver=TypeOfAverage::zero)const override;
    Complex angular_average(double s)const override;
    Complex eval_int(double s, TypeOfAverage aver=TypeOfAverage::zero)const override;

    array::ArrayT3Pi mandelstam_list;

    splined_path::SplinedPathT3Pi polar_egg;
    splined_path::SplinedPathDerivativeT3Pi egg_deriv; 
    
    void build_tilde_functions() override;
    void build_matching_list() override; 
    void build_matched_tilde() override;
    void make_splines() override;
    void make_splines_matched() override;
};


} // angular_average

#endif // ANGULAR_AVERAGE_H
