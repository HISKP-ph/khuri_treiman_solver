// Interpolate the omnes functions for isospin I=0,1,2 in order to improve the performance

#ifndef SPLINED_OMNES_H
#define SPLINED_OMNES_H

#include "cauchy.h"
#include "phase.h"
#include "path.h"
#include "path_eta_pi_pi.h"
#include "omnes.h"
#include "constants.h"
#include "facilities.h"
#include "mandelstam.h"
#include "type_aliases.h"
#include "array.h"
#include "enums.h"

#include <array>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
using type_aliases::Complex;
using constants::pi;
using namespace path_eta_pi_pi;
using namespace constants;
using namespace enums;

namespace splined_omnes {

/// Parent class to spline the Omnès functions.
class OmnesSpline{
public:
    OmnesSpline(double mass_1, double mass_2, double mass_3, double mass_4, double cutoff);
    ///<@param mass_1, mass_2, mass_3 masses of the decay products
    ///<@param mass_4 decay mass
    ///<@param cutoff numerical cutoff used in the dispersion integral
    ///<@param phase_I the complete path name of the phase shift with isospin I that shall be imported
    

    
    /// @brief Interpolated Omnès functions for isospin I.
    ///
    /// The parameter `Setting` can be either set to `above`, `below` or `egg`
    /// and decides whether to evaluate the Omnès function infinitesimally above the cut,
    /// below the cut or along the egg-like part. For the latter case the Omnès function
    /// is expressed in terms of a real-valued curve-parameter describing the complex contour.
    virtual Complex operator()(double curve_parameter, int isospin, Setting evaluation=Setting::above)const=0;

    /// scattering threshold in s
    double s_I;
    /// scattering threshold in t
    double t_I;
    double infinitesimal=0.00001;

    path::PolarEgg polar_egg;
    path_eta_pi_pi::Path gamma;

   

private:
    /// Function that builds the Omnès functions.
    virtual void build_omnes()=0;
    /// Function that splines the Omnès functions.
    virtual void spline_omnes()=0;
};

class OmnesSplineEta3Pi: public OmnesSpline{
public:
    OmnesSplineEta3Pi(double mass_1, double mass_4, double cutoff, 
                phase::PhaseEta3Pi phases);
    Complex operator()(double curve_parameter, int isospin, Setting evaluation=Setting::above)const;
private:
    void build_omnes() override;
    void spline_omnes() override;
    array::ArrayEta3Pi grid;

    phase::PhaseEta3Pi ph;
        
    /// Evaluate Omnes function of each isospin 0 for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_0_above;
    /// Evaluate Omnes function of each isospin 1 for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_1_above;
    /// Evaluate Omnes function of each isospin 2 for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_2_above;
    
    /// Evaluate Omnes function of each isospin 0 for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_0_below;
    /// Evaluate Omnes function of each isospin 1 for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_1_below;
    /// Evaluate Omnes function of each isospin 2 for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_2_below;
    
    /// Evaluate Omnes function of each isospin 0 for complex arguments on egg-like contour
    std::vector<Complex>omnes_0_egg;
    /// Evaluate Omnes function of each isospin 1 for complex arguments on egg-like contour
    std::vector<Complex>omnes_1_egg;
    /// Evaluate Omnes function of each isospin 2 for complex arguments on egg-like contour
    std::vector<Complex>omnes_2_egg;

    cauchy::Interpolate spline_omnes_0_above;
    cauchy::Interpolate spline_omnes_1_above;
    cauchy::Interpolate spline_omnes_2_above;
    
    cauchy::Interpolate spline_omnes_0_below;
    cauchy::Interpolate spline_omnes_1_below;
    cauchy::Interpolate spline_omnes_2_below;
    
    cauchy::Interpolate spline_omnes_0_egg;
    cauchy::Interpolate spline_omnes_1_egg;
    cauchy::Interpolate spline_omnes_2_egg;
};

class OmnesSplineEtap3Pi: public OmnesSpline{
public:
    OmnesSplineEtap3Pi(double mass_1, double mass_4, double cutoff, 
                phase::PhaseEtap3Pi phases);
    Complex operator()(double curve_parameter, int isospin, Setting evaluation=Setting::above)const;
private:
    void build_omnes() override;
    void spline_omnes() override;
    array::ArrayEtap3Pi grid;

    phase::PhaseEtap3Pi ph;
        
    /// Evaluate Omnes function of isospin 0 for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_0_above;
    /// Evaluate Omnes function of isospin 1 for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_1_above;
    /// Evaluate Omnes function of isospin 2 for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_2_above;
    
    /// Evaluate Omnes function of isospin 0 for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_0_below;
    /// Evaluate Omnes function of isospin 1 for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_1_below;
    /// Evaluate Omnes function of isospin 2 for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_2_below;
    
    /// Evaluate Omnes function of isospin 0 for complex arguments on egg-like contour
    std::vector<Complex>omnes_0_egg;
    /// Evaluate Omnes function of isospin 1 for complex arguments on egg-like contour
    std::vector<Complex>omnes_1_egg;
    /// Evaluate Omnes function of isospin 2 for complex arguments on egg-like contour
    std::vector<Complex>omnes_2_egg;

    cauchy::Interpolate spline_omnes_0_above;
    cauchy::Interpolate spline_omnes_1_above;
    cauchy::Interpolate spline_omnes_2_above;
    
    cauchy::Interpolate spline_omnes_0_below;
    cauchy::Interpolate spline_omnes_1_below;
    cauchy::Interpolate spline_omnes_2_below;
    
    cauchy::Interpolate spline_omnes_0_egg;
    cauchy::Interpolate spline_omnes_1_egg;
    cauchy::Interpolate spline_omnes_2_egg;
};

class OmnesSplineEtapEtaPiPi: public OmnesSpline{
public:
    OmnesSplineEtapEtaPiPi(double mass_1, double mass_3, double mass_4, double cutoff, 
                phase::PhaseEtapEtaPiPi phases);
    Complex operator()(double curve_parameter, int isospin, Setting evaluation=Setting::above)const;
private:
    void build_omnes() override;
    void spline_omnes() override;
    array::ArrayEtapEtaPiPi grid;

    phase::PhaseEtapEtaPiPi ph;
        
    /// Evaluate Omnes function of isospin 0 for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_0_above;
    /// Evaluate Omnes function of isospin 1 for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_1_above;
    /// Evaluate Omnes function of isospin 2 for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_2_above;
    
    /// Evaluate Omnes function of isospin 0 for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_0_below;
    /// Evaluate Omnes function of isospin 1 for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_1_below;
    /// Evaluate Omnes function of isospin 2 for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_2_below;

    cauchy::Interpolate spline_omnes_0_above;
    cauchy::Interpolate spline_omnes_1_above;
    cauchy::Interpolate spline_omnes_2_above;
    
    cauchy::Interpolate spline_omnes_0_below;
    cauchy::Interpolate spline_omnes_1_below;
    cauchy::Interpolate spline_omnes_2_below;
    
    /// Evaluate Omnes function of eta pi for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_eta_pi_above;
    /// Evaluate Omnes function of eta pi for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_eta_pi_below;
    
    /// Evaluate Omnes function of eta pi for complex arguments on egg-like contour
    std::vector<Complex>omnes_eta_pi_plus_upper_a;
    /// Evaluate Omnes function of eta pi for complex arguments on egg like contour
    std::vector<Complex>omnes_eta_pi_plus_upper_b;
    /// Evaluate Omnes function of eta pi for complex arguments on egg like contour
    std::vector<Complex>omnes_eta_pi_plus_lower_a;
    /// Evaluate Omnes function of eta pi for complex arguments on egg like contour
    std::vector<Complex>omnes_eta_pi_plus_lower_b;
    
    /// Evaluate Omnes function of eta pi for complex arguments on egg like contour
    std::vector<Complex>omnes_eta_pi_zero_upper_a;
    /// Evaluate Omnes function of eta pi for complex arguments on egg like contour
    std::vector<Complex>omnes_eta_pi_zero_upper_b;
    /// Evaluate Omnes function of eta pi for complex arguments on egg like contour
    std::vector<Complex>omnes_eta_pi_zero_lower_a;
    /// Evaluate Omnes function of eta pi for complex arguments on egg like contour
    std::vector<Complex>omnes_eta_pi_zero_lower_b;
    
    /// Evaluate Omnes function of isospin 0 for complex arguments on egg like contour
    std::vector<Complex>omnes_0_minus_upper_a;
    /// Evaluate Omnes function of isospin 0 for complex arguments on egg like contour
    std::vector<Complex>omnes_0_minus_upper_b;
    /// Evaluate Omnes function of isospin 0 for complex arguments on egg like contour
    std::vector<Complex>omnes_0_minus_lower_a;
    /// Evaluate Omnes function of isospin 0 for complex arguments on egg like contour
    std::vector<Complex>omnes_0_minus_lower_b;

    /// Evaluate Omnes function of isospin 1 for complex arguments on egg like contour (only necessary for C-violation)
    std::vector<Complex>omnes_1_minus_upper_a;
    /// Evaluate Omnes function of isospin 1 for complex arguments on egg like contour (only necessary for C-violation)
    std::vector<Complex>omnes_1_minus_upper_b;
    /// Evaluate Omnes function of isospin 1 for complex arguments on egg like contour (only necessary for C-violation)
    std::vector<Complex>omnes_1_minus_lower_a;
    /// Evaluate Omnes function of isospin 1 for complex arguments on egg like contour (only necessary for C-violation)
    std::vector<Complex>omnes_1_minus_lower_b;
    
    
    cauchy::Interpolate spline_omnes_eta_pi_above;
    cauchy::Interpolate spline_omnes_eta_pi_below;
    
    cauchy::Interpolate spline_omnes_eta_pi_plus_upper_a;
    cauchy::Interpolate spline_omnes_eta_pi_plus_upper_b;
    cauchy::Interpolate spline_omnes_eta_pi_plus_lower_a;
    cauchy::Interpolate spline_omnes_eta_pi_plus_lower_b;
    
    cauchy::Interpolate spline_omnes_eta_pi_zero_upper_a;
    cauchy::Interpolate spline_omnes_eta_pi_zero_upper_b;
    cauchy::Interpolate spline_omnes_eta_pi_zero_lower_a;
    cauchy::Interpolate spline_omnes_eta_pi_zero_lower_b;
    
    cauchy::Interpolate spline_omnes_0_minus_upper_a;
    cauchy::Interpolate spline_omnes_0_minus_upper_b;
    cauchy::Interpolate spline_omnes_0_minus_lower_a;
    cauchy::Interpolate spline_omnes_0_minus_lower_b;

    cauchy::Interpolate spline_omnes_1_minus_upper_a;
    cauchy::Interpolate spline_omnes_1_minus_upper_b;
    cauchy::Interpolate spline_omnes_1_minus_lower_a;
    cauchy::Interpolate spline_omnes_1_minus_lower_b;

};

class OmnesSplineV3Pi: public OmnesSpline{
public:
    OmnesSplineV3Pi(double mass_1, double mass_4, double cutoff, 
                phase::PhaseV3Pi phases);
    Complex operator()(double curve_parameter, int isospin, Setting evaluation=Setting::above)const;
private:
    void build_omnes() override;
    void spline_omnes() override;
    array::ArrayV3Pi grid;

    phase::PhaseV3Pi ph;

    /// Evaluate Omnes function for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_above;
    /// Evaluate Omnes function for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_below;
    /// Evaluate Omnes function for complex arguments on egg-like contour
    std::vector<Complex>omnes_egg;


    cauchy::Interpolate spline_omnes_above;
    
    cauchy::Interpolate spline_omnes_below;
    
    cauchy::Interpolate spline_omnes_egg;
};

class OmnesSplineX3Pi: public OmnesSpline{
public:
    OmnesSplineX3Pi(double mass_1, double mass_4, double cutoff, 
                phase::PhaseX3Pi phases);
    Complex operator()(double curve_parameter, int isospin, Setting evaluation=Setting::above)const;
private:
    void build_omnes() override;
    void spline_omnes() override;
    array::ArrayX3Pi grid;

    phase::PhaseX3Pi ph;

    /// Evaluate Omnes function for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_above;
    /// Evaluate Omnes function for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_below;
    /// Evaluate Omnes function for complex arguments on egg-like contour
    std::vector<Complex>omnes_egg;


    cauchy::Interpolate spline_omnes_above;
    
    cauchy::Interpolate spline_omnes_below;
    
    cauchy::Interpolate spline_omnes_egg;
};

class OmnesSplineT3Pi: public OmnesSpline{
public:
    OmnesSplineT3Pi(double mass_1, double mass_4, double cutoff, 
                phase::PhaseT3Pi phases);
    Complex operator()(double curve_parameter, int isospin, Setting evaluation=Setting::above)const;
private:
    void build_omnes() override;
    void spline_omnes() override;
    array::ArrayT3Pi grid;

    phase::PhaseT3Pi ph;

    /// Evaluate Omnes function for complex arguments infinitesimally above the cut
    std::vector<Complex>omnes_above;
    /// Evaluate Omnes function for complex arguments infinitesimally below the cut
    std::vector<Complex>omnes_below;
    /// Evaluate Omnes function for complex arguments on egg-like contour
    std::vector<Complex>omnes_egg;


    cauchy::Interpolate spline_omnes_above;
    
    cauchy::Interpolate spline_omnes_below;
    
    cauchy::Interpolate spline_omnes_egg;
};


}

#endif // SPLINED_OMNES_H
