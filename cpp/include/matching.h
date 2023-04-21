#ifndef MATCHING_H
#define MATCHING_H


#include "helpers.h"
#include "type_aliases.h"
#include "facilities.h"
#include "gsl_interface.h"
#include <cmath>
#include <complex>
#include <vector>

namespace match {

using type_aliases::Complex;
using gsl::Interpolate;
using namespace std::complex_literals;

/// @brief Class for the matching procedure around the various thresholds.
class Matching {
public:
    Matching(double mass_1, double mass_2, double mass_3, double mass_4,
             double epsilon, double validity);
        ///< @param mass_1, mass_2, mass_3 masses of the decay products
        ///< @param mass_4 decay mass
        ///< @param epsilon the value at which an expansion within the dispersion integral is
        /// matched to the exact value
        ///< @param validity the range of validity within which the expansion is used
    
    /// Function for the matching of the tilde function for S-waves.
    double s_match(double s, const Interpolate& amplitude);
    /// Function for the matching of the tilde function for P-waves.
    double p_match(double s, const Interpolate& amplitude);
    /// Function for the matching of the tilde function for D-waves.
    double d_match(double s, const Interpolate& amplitude);
    
    
    
    
    /// Gives the value of the derivative needed in the matching of the dispersion integral.
    double p_wave_parameter(double s, const Interpolate& amplitude);
    /// Gives the value of the first and second derivative needed in the matching of the dispersion integral.
    std::tuple<double,double> d_wave_parameters(double s, const Interpolate& amplitude);

    /// Matched expansion for the dispersion relation around the pseudo threshold.
    Complex s_wave_integrand_expansion(double s, const Interpolate& amplitude);
    /// Matched expansion for the dispersion relation around the pseudo threshold.
    Complex p_wave_integrand_expansion(double s, const Interpolate& amplitude);
    /// Matched expansion for the dispersion relation around the pseudo threshold.
    Complex d_wave_integrand_expansion(double s, const Interpolate& amplitude);
    
    
    
private:
    double mass_1, mass_2, mass_3, mass_4;
    double epsilon, validity;
    double threshold_12, threshold_34, pseudo_threshold_34;    
    
    /// Shortens the notation for the matched coefficients, parameters do not have a special meaning.
    double parameter_simplification(const Interpolate& amplitude, double s, double i,
                                    double j, double k, double l, double m, double s1, double n, double p, double q);
    
    /// Contains exceptions based on wrong input.
    void errors_matching();
};



} // match


#endif // MATCHING_H
