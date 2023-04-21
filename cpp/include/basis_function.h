#ifndef BASIS_FUNCTION_H
#define BASIS_FUNCTION_H


#include "dispersion_integral.h"
#include "path.h"
#include "type_aliases.h"
#include "facilities.h"
#include "angular_average.h"
#include "gsl_interface.h"
#include "path_eta_pi_pi.h"
#include "enums.h"
#include "splined_omnes.h"

#include <complex>
#include <cmath>


///Calculating the basis functions. Used as an intermediate step to discriminate all the different cases.
namespace basis{

using type_aliases::Complex;
using constants::pi;
using namespace angular_average;
using namespace disp;
using type_aliases::Complex;
using type_aliases::CFunction;
using namespace enums;
using namespace std::complex_literals;


using FunctionSetOmnes=const std::function<Complex(double, int, Setting)>;
// for Omnès function
using FunctionSet=const std::function<Complex(double, int, SubtractionConstant)>&;
// for tilde function
using LongFunctionSet=const std::function<Complex(double, int, SubtractionConstant, Setting)>&;
// for basis function


/// Parent Class to initialize the interation procedure. The homogeneous solution is solely given in terms of the Omnès function.
class HomogeneousSolution{
public:
    HomogeneousSolution(double mass_1, double mass_2, double mass_3, double mass_4,
                        FunctionSetOmnes omnes);
    ///<@param mass_1, mass_2, mass_3 masses of the decay products
    ///<@param mass_4 decay mass
    ///<@param omnes Omnès function
    
    /// Initial 'guesses' for the basis amplitudes depending on the subtraction constant.
    virtual Complex initial_guess(double s, int isospin, SubtractionConstant sub, Setting evaluation)const=0;  
    FunctionSetOmnes omnes;   
};

class HomogeneousSolutionEta3Pi: public HomogeneousSolution{
public:
    HomogeneousSolutionEta3Pi(double mass_1, double mass_4,
                        FunctionSetOmnes omnes);

    Complex initial_guess(double s, int isospin, SubtractionConstant sub, Setting evaluation)const override;   
    
private:
    path::PolarEgg polar_egg;   
};

class HomogeneousSolutionEtap3Pi: public HomogeneousSolution{
public:
    HomogeneousSolutionEtap3Pi(double mass_1, double mass_4,
                        FunctionSetOmnes omnes);

    Complex initial_guess(double s, int isospin, SubtractionConstant sub, Setting evaluation)const override;   
    
private:
    path::PolarEgg polar_egg;  
};

class HomogeneousSolutionEtapEtaPiPi: public HomogeneousSolution{
public:
    HomogeneousSolutionEtapEtaPiPi(double mass_1, double mass_3, double mass_4,
                        FunctionSetOmnes omnes);

    Complex initial_guess(double s, int isospin, SubtractionConstant sub, Setting evaluation)const override;   
    
private:
    path_eta_pi_pi::Path gamma;    
};

class HomogeneousSolutionV3Pi: public HomogeneousSolution{
public:
    HomogeneousSolutionV3Pi(double mass_1, double mass_4,
                        FunctionSetOmnes omnes);

    Complex initial_guess(double s, int isospin, SubtractionConstant sub, Setting evaluation)const override;   
    
private:
    path::PolarEgg polar_egg;   
};

class HomogeneousSolutionX3Pi: public HomogeneousSolution{
public:
    HomogeneousSolutionX3Pi(double mass_1, double mass_4,
                        FunctionSetOmnes omnes);

    Complex initial_guess(double s, int isospin, SubtractionConstant sub, Setting evaluation)const override;   
    
private:
    path::PolarEgg polar_egg;   
};

class HomogeneousSolutionT3Pi: public HomogeneousSolution{
public:
    HomogeneousSolutionT3Pi(double mass_1, double mass_4,
                        FunctionSetOmnes omnes);

    Complex initial_guess(double s, int isospin, SubtractionConstant sub, Setting evaluation)const override;   
    
private:
    path::PolarEgg polar_egg;   
};


///Parent Class to setup all the different angular averages.
class BasisTilde{    
public:
    BasisTilde() = default;    
    
    ///Returns the matched angular average from angular_average::AngularAverage.
    virtual Complex matched_tilde(double s, int isospin, SubtractionConstant sub)const =0;
};

class BasisTildeEta3Pi: public BasisTilde{
private:
    //Needed to access angular averages for each set of basis amplitudes and each isospin
    const AngularAverageEta3Pi ang_a0_0;
    const AngularAverageEta3Pi ang_a0_1;
    const AngularAverageEta3Pi ang_a0_2;
    const AngularAverageEta3Pi ang_b0_0;
    const AngularAverageEta3Pi ang_b0_1;
    const AngularAverageEta3Pi ang_b0_2;
    const AngularAverageEta3Pi ang_a1_0;
    const AngularAverageEta3Pi ang_a1_1;
    const AngularAverageEta3Pi ang_a1_2;
    //C-violating
    const AngularAverageEta3PiC ang_g1_1;
    const AngularAverageEta3PiC ang_h1_1;
    const AngularAverageEta3PiC ang_h1_2;

    
public:
    BasisTildeEta3Pi(double mass_1, double mass_4,
               LongFunctionSet amplitude,
               double cutoff, int iteration, double epsilon, double validity);
        ///<@param mass_1 mass of the decay products
        ///<@param mass_4 decay mass
        ///<@param amplitude basis functions
        ///<@param cutoff numerical cutoff for the dispersion integral
        ///<@param iteration if iteration=0, use the fact that the Omnès function is perfectly Schwartz
        ///<@param epsilon, validity matching point and range for the matching procedure of the angular average

    Complex matched_tilde(double s, int isospin, SubtractionConstant sub)const override ;
};

class BasisTildeEtap3Pi: public BasisTilde{
private:
    //Needed to access angular averages for each set of basis amplitudes and each isospin
    const AngularAverageEtap3Pi ang_a0_0;
    const AngularAverageEtap3Pi ang_a0_1;
    const AngularAverageEtap3Pi ang_a0_2;
    const AngularAverageEtap3Pi ang_b0_0;
    const AngularAverageEtap3Pi ang_b0_1;
    const AngularAverageEtap3Pi ang_b0_2;
    const AngularAverageEtap3Pi ang_a1_0;
    const AngularAverageEtap3Pi ang_a1_1;
    const AngularAverageEtap3Pi ang_a1_2;
    //C-violating
    const AngularAverageEtap3PiC ang_g1_1;
    const AngularAverageEtap3PiC ang_h1_1;
    const AngularAverageEtap3PiC ang_h1_2;

    
public:
    BasisTildeEtap3Pi(double mass_1, double mass_4,
               LongFunctionSet amplitude,
               double cutoff, int iteration, double epsilon, double validity);
        ///<@param mass_1 mass of the decay products
        ///<@param mass_4 decay mass
        ///<@param amplitude basis functions
        ///<@param cutoff numerical cutoff for the dispersion integral
        ///<@param iteration if iteration=0, use the fact that the Omnès function is perfectly Schwartz
        ///<@param epsilon, validity matching point and range for the matching procedure of the angular average

    Complex matched_tilde(double s, int isospin, SubtractionConstant sub)const override;
};

// Analogously for eta'-> eta pi pi
class BasisTildeEtapEtaPiPi: public BasisTilde{
private:
    //Needed to access angular averages for each set of basis amplitudes and each isospin
    AngularAverageEtapEtaPiPi ang_a0_000;
    AngularAverageEtapEtaPiPi ang_a0_010;
    AngularAverageEtapEtaPiPi ang_b0_000;
    AngularAverageEtapEtaPiPi ang_b0_010;
    AngularAverageEtapEtaPiPi ang_c0_000;
    AngularAverageEtapEtaPiPi ang_c0_010;
    AngularAverageEtapEtaPiPi ang_d0_000;
    AngularAverageEtapEtaPiPi ang_d0_010;
    //C-violating
    AngularAverageEtapEtaPiPiC ang_a1_111;
    AngularAverageEtapEtaPiPiC ang_a1_110;
    AngularAverageEtapEtaPiPiC ang_b1_111;
    AngularAverageEtapEtaPiPiC ang_b1_110;


    
public:
    BasisTildeEtapEtaPiPi(double mass_1, double mass_3, double mass_4,
               LongFunctionSet amplitude,
               double cutoff, int iteration, double epsilon, double validity);
        ///<@param mass_1, mass_3 masses of the decay products 
        ///<@param mass_4 decay mass
        ///<@param amplitude basis functions
        ///<@param cutoff numerical cutoff for the dispersion integral
        ///<@param iteration if iteration=0, use the fact that the Omnès function is perfectly Schwartz
        ///<@param epsilon, validity matching point and range for the matching procedure of the angular average

    Complex matched_tilde(double s, int isospin, SubtractionConstant sub)const override;
};

class BasisTildeV3Pi: public BasisTilde{
private:
    //Needed to access angular averages for each set of basis amplitudes and each isospin
    const AngularAverageV3Pi ang_a1;
    const AngularAverageV3Pi ang_b1;

    
public:
    BasisTildeV3Pi(double mass_1, double mass_4,
               LongFunctionSet amplitude,
               double cutoff, int iteration, double epsilon, double validity);
        ///<@param mass_1 mass of the decay products
        ///<@param mass_4 decay mass
        ///<@param amplitude basis functions
        ///<@param cutoff numerical cutoff for the dispersion integral
        ///<@param iteration if iteration=0, use the fact that the Omnès function is perfectly Schwartz
        ///<@param epsilon, validity matching point and range for the matching procedure of the angular average

    Complex matched_tilde(double s, int isospin, SubtractionConstant sub)const override ;
};

class BasisTildeX3Pi: public BasisTilde{
private:
    //Needed to access angular averages for each set of basis amplitudes and each isospin
    const AngularAverageX3Pi ang_a1;
    const AngularAverageX3Pi ang_b1;

    
public:
    BasisTildeX3Pi(double mass_1, double mass_4,
               LongFunctionSet amplitude,
               double cutoff, int iteration, double epsilon, double validity);
        ///<@param mass_1 mass of the decay products
        ///<@param mass_4 decay mass
        ///<@param amplitude basis functions
        ///<@param cutoff numerical cutoff for the dispersion integral
        ///<@param iteration if iteration=0, use the fact that the Omnès function is perfectly Schwartz
        ///<@param epsilon, validity matching point and range for the matching procedure of the angular average

    Complex matched_tilde(double s, int isospin, SubtractionConstant sub)const override ;
};

class BasisTildeT3Pi: public BasisTilde{
private:
    //Needed to access angular averages for each set of basis amplitudes and each isospin
    const AngularAverageT3Pi ang_a1;
    const AngularAverageT3Pi ang_b1;

    
public:
    BasisTildeT3Pi(double mass_1, double mass_4,
               LongFunctionSet amplitude,
               double cutoff, int iteration, double epsilon, double validity);
        ///<@param mass_1 mass of the decay products
        ///<@param mass_4 decay mass
        ///<@param amplitude basis functions
        ///<@param cutoff numerical cutoff for the dispersion integral
        ///<@param iteration if iteration=0, use the fact that the Omnès function is perfectly Schwartz
        ///<@param epsilon, validity matching point and range for the matching procedure of the angular average

    Complex matched_tilde(double s, int isospin, SubtractionConstant sub)const override ;
};



/// Parent class to calculate the basis amplitudes with the angular averages from BasisTilde.
class BasisAmplitude{
public:
    BasisAmplitude(double mass_1, double mass_2, double mass_3, double mass_4,
                   FunctionSetOmnes omnes,
                   FunctionSet matched_tilde,
                   double cutoff, double epsilon, double validity, int max_subs);
        ///<@param mass_1, mass_2, mass_3 masses of the decay products
        ///<@param mass_4 decay mass
        ///<@param omnes Omnès function
        ///<@param matched_tilde matched angular average from BasisTilde
        ///<@param cutoff numerical cutoff for the dispersion integral
        ///<@param epsilon, validity matching point and range for the matching procedure of the angular average
        ///<@param max_subs maximum number of subtraction constants (not valid for all of the processes)
    
    
    
    /// Function that returns the basis amplitudes.
    virtual Complex amplitude(double s, int isospin, SubtractionConstant sub,
                      Setting evaluation)const=0;
    
    FunctionSetOmnes omnes;
    int max_subs;
private:
    ///Function that calculates the dispersion integral via disp::Dispersion.
    virtual Complex dispersion_integral(double s, int isospin,
                                SubtractionConstant sub, Setting evaluation)const=0;
    /// Helper function that adds the polynomial and the dispersion integral for the function amplitude.
    virtual Complex helper(double s, int isospin, SubtractionConstant sub, Setting evaluation,
                   double number_subtractions, double power_of_s=0.,
                   double polynomial_on_off=0.)const=0;
};

class BasisAmplitudeEta3Pi: public BasisAmplitude{
public:
    BasisAmplitudeEta3Pi(double mass_1, double mass_4,
                   FunctionSetOmnes omnes,
                   FunctionSet matched_tilde,
                   phase::PhaseEta3Pi phases,
                   double cutoff,double epsilon, double validity, int max_subs);   
    
    Complex amplitude(double s, int isospin, SubtractionConstant sub,
                      Setting evaluation)const override;
    
    
private:
    path::PolarEgg polar_egg;

    //Needed to access dispersion integrals for isospin 0, 1 and 2
    DispersionEta3Pi disp_a0_0;
    DispersionEta3Pi disp_a0_1;
    DispersionEta3Pi disp_a0_2;
    DispersionEta3Pi disp_b0_0;
    DispersionEta3Pi disp_b0_1;
    DispersionEta3Pi disp_b0_2;
    DispersionEta3Pi disp_a1_0;
    DispersionEta3Pi disp_a1_1;
    DispersionEta3Pi disp_a1_2;
    //C-violating
    DispersionEta3Pi disp_g1_1;
    DispersionEta3Pi disp_h1_1;
    DispersionEta3Pi disp_h1_2;
    
    
    Complex dispersion_integral(double s, int isospin,
                                SubtractionConstant sub, Setting evaluation)const override;
    Complex helper(double s, int isospin, SubtractionConstant sub, Setting evaluation,
                   double number_subtractions, double power_of_s=0.,
                   double polynomial_on_off=0.)const override;
};

class BasisAmplitudeEtap3Pi: public BasisAmplitude{
public:
    BasisAmplitudeEtap3Pi(double mass_1, double mass_4,
                   FunctionSetOmnes omnes,
                   FunctionSet matched_tilde,
                   phase::PhaseEtap3Pi phases,
                   double cutoff,double epsilon, double validity, int max_subs);   
    
    Complex amplitude(double s, int isospin, SubtractionConstant sub,
                      Setting evaluation)const override;
    
    
private:
    path::PolarEgg polar_egg;

    //Needed to access dispersion integrals for isospin 0, 1 and 2
    DispersionEtap3Pi disp_a0_0;
    DispersionEtap3Pi disp_a0_1;
    DispersionEtap3Pi disp_a0_2;
    DispersionEtap3Pi disp_b0_0;
    DispersionEtap3Pi disp_b0_1;
    DispersionEtap3Pi disp_b0_2;
    DispersionEtap3Pi disp_a1_0;
    DispersionEtap3Pi disp_a1_1;
    DispersionEtap3Pi disp_a1_2;
    //C-violating
    DispersionEtap3Pi disp_g1_1;
    DispersionEtap3Pi disp_h1_1;
    DispersionEtap3Pi disp_h1_2;
    
    
    Complex dispersion_integral(double s, int isospin,
                                SubtractionConstant sub, Setting evaluation)const override;
    Complex helper(double s, int isospin, SubtractionConstant sub, Setting evaluation,
                   double number_subtractions, double power_of_s=0.,
                   double polynomial_on_off=0.)const override;
};

class BasisAmplitudeEtapEtaPiPi: public BasisAmplitude{
public:
    BasisAmplitudeEtapEtaPiPi(double mass_1, double mass_3, double mass_4,
                   FunctionSetOmnes omnes,
                   FunctionSet matched_tilde,
                   phase::PhaseEtapEtaPiPi phases,
                   double cutoff, double epsilon, double validity, int max_subs);   
    
    Complex amplitude(double s, int isospin, SubtractionConstant sub,
                      Setting evaluation)const override;
    
    
private:
    path_eta_pi_pi::Path gamma;
    
    //Needed to access dispersion integrals for isospin 0, 1 and 2
    DispersionEtapEtaPiPi disp_a0_000;
    DispersionEtapEtaPiPi disp_a0_010;
    DispersionEtapEtaPiPi disp_b0_000;
    DispersionEtapEtaPiPi disp_b0_010;
    DispersionEtapEtaPiPi disp_c0_000;
    DispersionEtapEtaPiPi disp_c0_010;
    DispersionEtapEtaPiPi disp_d0_000;
    DispersionEtapEtaPiPi disp_d0_010;
    //C-violating
    DispersionEtapEtaPiPi disp_a1_111;
    DispersionEtapEtaPiPi disp_a1_110;
    DispersionEtapEtaPiPi disp_b1_111;
    DispersionEtapEtaPiPi disp_b1_110;
       
    Complex dispersion_integral(double s, int isospin,
                                SubtractionConstant sub, Setting evaluation)const override;
    
    Complex helper(double s, int isospin,
                                  SubtractionConstant sub, Setting evaluation,
                                  double number_subtractions, double power_of_s=0.,
                                  double polynomial_on_off=0)const override;

};

class BasisAmplitudeV3Pi: public BasisAmplitude{
public:
    BasisAmplitudeV3Pi(double mass_1, double mass_4,
                   FunctionSetOmnes omnes,
                   FunctionSet matched_tilde,
                   phase::PhaseV3Pi phases,
                   double cutoff,double epsilon, double validity, int max_subs);   
    
    Complex amplitude(double s, int isospin, SubtractionConstant sub,
                      Setting evaluation)const override;
    
private:
    path::PolarEgg polar_egg;

    //Needed to access dispersion integrals for isospin 0, 1 and 2
    DispersionV3Pi disp_a1;
    DispersionV3Pi disp_b1;    
    
    Complex dispersion_integral(double s, int isospin,
                                SubtractionConstant sub, Setting evaluation)const override;
    Complex helper(double s, int isospin, SubtractionConstant sub, Setting evaluation,
                   double number_subtractions, double power_of_s=0.,
                   double polynomial_on_off=0.)const override;
};

class BasisAmplitudeX3Pi: public BasisAmplitude{
public:
    BasisAmplitudeX3Pi(double mass_1, double mass_4,
                   FunctionSetOmnes omnes,
                   FunctionSet matched_tilde,
                   phase::PhaseX3Pi phases,
                   double cutoff,double epsilon, double validity, int max_subs);   
    
    Complex amplitude(double s, int isospin, SubtractionConstant sub,
                      Setting evaluation)const override;
    
private:
    path::PolarEgg polar_egg;

    //Needed to access dispersion integrals for isospin 0, 1 and 2
    DispersionX3Pi disp_a1;
    DispersionX3Pi disp_b1;    
    
    Complex dispersion_integral(double s, int isospin,
                                SubtractionConstant sub, Setting evaluation)const override;
    Complex helper(double s, int isospin, SubtractionConstant sub, Setting evaluation,
                   double number_subtractions, double power_of_s=0.,
                   double polynomial_on_off=0.)const override;
};

class BasisAmplitudeT3Pi: public BasisAmplitude{
public:
    BasisAmplitudeT3Pi(double mass_1, double mass_4,
                   FunctionSetOmnes omnes,
                   FunctionSet matched_tilde,
                   phase::PhaseT3Pi phases,
                   double cutoff,double epsilon, double validity, int max_subs);   
    
    Complex amplitude(double s, int isospin, SubtractionConstant sub,
                      Setting evaluation)const override;
    
private:
    path::PolarEgg polar_egg;

    //Needed to access dispersion integrals for isospin 0, 1 and 2
    DispersionT3Pi disp_a1;
    DispersionT3Pi disp_b1;    
    
    Complex dispersion_integral(double s, int isospin,
                                SubtractionConstant sub, Setting evaluation)const override;
    Complex helper(double s, int isospin, SubtractionConstant sub, Setting evaluation,
                   double number_subtractions, double power_of_s=0.,
                   double polynomial_on_off=0.)const override;
};



} // basis

#endif // BASIS_FUNCTION_H

