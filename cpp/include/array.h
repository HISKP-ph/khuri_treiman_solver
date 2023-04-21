// Define several arrays of Mandelstam s for which several parts of the amplitude are evaluated

#ifndef ARRAY_H
#define ARRAY_H

#include "path.h"
#include "enums.h"
#include <vector>
#include <cmath>
using namespace enums;
using namespace constants;



namespace array{

/// Parent class for generating arrays for the evaluation.
class Array{
public:
    Array(double mass_1, double mass_2, double mass_3, double mass_4, double cutoff);
    ///< @param mass_1, mass_2, mass_3 masses of the decay products
    ///< @param mass_4 decay mass
    ///< @param cutoff numerical cutoff used in the dispersion integral
 
    /// @brief Array in Mandelstam s starting at threshold for pipi scattering. more points close to cusps, use for I=1,2.
    /// Use this to evaluate the angular average.
    std::vector<double> interval_s_th;
    

    /// @brief Array in Mandelstam s use for I=1,2 amplitudes with more points close to thresholds.
    /// Use this to evaluate the omnes function and the basis amplitude.
    std::vector<double> interval_s;
    /// @brief Array in Mandelstam s. use for I=0 amplitudes.
    /// Additionally more points close to the cusp in the S-wave phase shift.
    /// Not used in the current version, but could be included later.
    std::vector<double> interval_s_0;
    
    
    /// @brief Additional array to evaluate the dispersion integral.
    /// Here points very close to the pseudo threshold are left out.
    std::vector<double> interval_disp;
    /// @brief List used to expand the interval. Only evaluate the expansion along these short lists to save computation time.
    std::vector<double> interval_disp_expansion;
    
    /// Used to expand the integrand of the dispersion  integral in a small range around the pseudo threshold.
    double eps=0.00001;
    /// scattering threshold (lower limit)
    double s_I;
    /// transition to region II, where s- hits sth
    double s_II;
    /// pseudo threshold
    double s_III;
    /// scattering threshold (upper limit)
    double s_IV;
    double t_I, t_III, t_IV;
    
    double decay_mass;
    double cutoff;

    /// function that generates list with constant step_size.
    void constant_spacing(std::vector<double>& list, double min, double max, double step_size);
    /// function that generates list with decreasing step_size towards min,max.
    void asymptotic_spacing(std::vector<double>& list, double min, double max);

private:
    /// Function that creats the intervals.    
    virtual void create_interval()=0;
};

class ArrayEta3Pi: public Array{
public:
    ArrayEta3Pi(double mass_1, double mass_4, double cutoff);    
    
    /// Array for the egg-like part of the complex integration contour for the curve parameter from 0 to 2pi.
    std::vector<double> interval_phi;
    /// Array which only includes the values of Mandelstam s at which the complex contour starts and ends.
    std::vector<double> interval_egg;
private:    
    void create_interval() override;
};

class ArrayEtap3Pi: public Array{
public:
    ArrayEtap3Pi(double mass_1, double mass_4, double cutoff);

    /// Array for the egg-like part of the complex integration contour for the curve parameter from 0 to 2pi.
    std::vector<double> interval_phi;
    /// Array which only includes the values of Mandelstam s at which the complex contour starts and ends.
    std::vector<double> interval_egg;
private:    
    void create_interval() override;
};

class ArrayEtapEtaPiPi: public Array{
public:
    ArrayEtapEtaPiPi(double mass_1, double mass_3, double mass_4, double cutoff);
    
    /// Analogous to the s lists, but here for Mandelstam t.
    std::vector<double> interval_t;
    /// Analogous to the s lists, but here for Mandelstam t.
    std::vector<double> interval_t_th;
    
    /// Analogous to the s lists, but here for Mandelstam t.
    std::vector<double> interval_t_disp;
    /// Analogous to the s lists, but here for Mandelstam t.
    std::vector<double> interval_t_disp_expansion;
    
    /// Array for the egg-like part in eta'-> eta pi pi.
    std::vector<double> interval_y;
    /// Array only for t in region III, used to spline path and its derivative.
    std::vector<double> interval_t_path_spline;    
private:    
    void create_interval() override;
};

class ArrayV3Pi: public Array{
public:
    ArrayV3Pi(double mass_1, double mass_4, double cutoff);    
    
    /// Array for the egg-like part of the complex integration contour for the curve parameter from 0 to 2pi.
    std::vector<double> interval_phi;
    /// Array which only includes the values of Mandelstam s at which the complex contour starts and ends.
    std::vector<double> interval_egg;
private:    
    void create_interval() override;
};

class ArrayX3Pi: public Array{
public:
    ArrayX3Pi(double mass_1, double mass_4, double cutoff);    
    
    /// Array for the egg-like part of the complex integration contour for the curve parameter from 0 to 2pi.
    std::vector<double> interval_phi;
    /// Array which only includes the values of Mandelstam s at which the complex contour starts and ends.
    std::vector<double> interval_egg;
private:    
    void create_interval() override;
};

class ArrayT3Pi: public Array{
public:
    ArrayT3Pi(double mass_1, double mass_4, double cutoff);    
    
    /// Array for the egg-like part of the complex integration contour for the curve parameter from 0 to 2pi.
    std::vector<double> interval_phi;
    /// Array which only includes the values of Mandelstam s at which the complex contour starts and ends.
    std::vector<double> interval_egg;
private:    
    void create_interval() override;
};


} // array

#endif // ARRAY_H
