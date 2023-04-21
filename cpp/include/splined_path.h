// Interpolate the egg-like part of the complex integration contour to improve the performance

#ifndef SPLINED_PATH_H
#define SPLINED_PATH_H

#include "cauchy.h"
#include "phase.h"
#include "path.h"
#include "omnes.h"
#include "constants.h"
#include "facilities.h"
#include "mandelstam.h"
#include "type_aliases.h"
#include "array.h"

#include <array>
#include <cmath>
#include <complex>
#include <vector>
using type_aliases::Complex;

namespace splined_path {
/// @brief Parent Class to interpolate the path, additionally define some
/// integration limits for a comfortable use in angular_averages.h.
class SplinedPath{   
public:
    SplinedPath(double mass_1, double mass_2, double mass_3,
                double mass_4, double cutoff);
    ///< @param mass_1, mass_2, mass_3 masses of the decay products
    ///< @param mass_4 decay mass
    ///< @param cutoff numerical cutoff for the dispersion integral
    
    double s_I, s_III, s_IV, s_0;

    /// list for the complex contour
    std::vector<Complex>egg_list;
    /// list for the lower limits of the curve parameter    
    std::vector<double>phi_lower_list;
    /// list for the upper limits of the curve parameter
    std::vector<double>phi_upper_list;

    path::PolarEgg egg;
    
    /// interpolation of the egg-like complex contour
    cauchy::Interpolate spline_egg;
    /// interpolation of the lower limits of the curve parameter
    gsl::Interpolate spline_phi_lower;
    /// interpolation of the upper limits of the curve parameter
    gsl::Interpolate spline_phi_upper;
    /// Function that builds the path.
    virtual void build_egg()=0;
    /// Function that splines the path.
    virtual void spline_path()=0;

    double s_upper(double s)const;
    double s_lower(double s)const;

private:
    /// Kinematic function to determine the integration limits for the polar-egg.
    Complex kappa(double s)const;

};

class SplinedPathEta3Pi: public SplinedPath{   
public:
    SplinedPathEta3Pi(double mass_1, double mass_4, double cutoff);
    ///< @param mass_1 mass of the decay products
    ///< @param mass_4 decay mass
    ///< @param cutoff numerical cutoff for the dispersion integral
    
    array::ArrayEta3Pi s_list;
private:
    void build_egg() override;
    void spline_path() override;
};

class SplinedPathEtap3Pi: public SplinedPath{   
public:
    SplinedPathEtap3Pi(double mass_1, double mass_4, double cutoff);
    ///< @param mass_1 mass of the decay products
    ///< @param mass_4 decay mass
    ///< @param cutoff numerical cutoff for the dispersion integral
    
    array::ArrayEtap3Pi s_list;
private:
    void build_egg() override;
    void spline_path() override;
};

class SplinedPathV3Pi: public SplinedPath{   
public:
    SplinedPathV3Pi(double mass_1, double mass_4, double cutoff);
    ///< @param mass_1 mass of the decay products
    ///< @param mass_4 decay mass
    ///< @param cutoff numerical cutoff for the dispersion integral
    
    array::ArrayV3Pi s_list;
private:
    void build_egg() override;
    void spline_path() override;
};

class SplinedPathX3Pi: public SplinedPath{   
public:
    SplinedPathX3Pi(double mass_1, double mass_4, double cutoff);
    ///< @param mass_1 mass of the decay products
    ///< @param mass_4 decay mass
    ///< @param cutoff numerical cutoff for the dispersion integral
    
    array::ArrayX3Pi s_list;
private:
    void build_egg() override;
    void spline_path() override;
};

class SplinedPathT3Pi: public SplinedPath{   
public:
    SplinedPathT3Pi(double mass_1, double mass_4, double cutoff);
    ///< @param mass_1 mass of the decay products
    ///< @param mass_4 decay mass
    ///< @param cutoff numerical cutoff for the dispersion integral
    
    array::ArrayT3Pi s_list;
private:
    void build_egg() override;
    void spline_path() override;
};


/// Parent Class to interpolate the derivative of the path.
class SplinedPathDerivative{
public:
    SplinedPathDerivative() = default;

    /// Function that builds the derivative.
    virtual void build_derivative()=0;
    /// Function that splines the derivative.
    virtual void spline_derivative()=0;

    /// List of the derivative of the complex contour.
    std::vector<Complex>egg_derivative;
    /// Interpolation of the derivative of the egg-like complex contour.
    cauchy::Interpolate spline_egg_derivative;
};

class SplinedPathDerivativeEta3Pi: public SplinedPathDerivative{
public:
    SplinedPathDerivativeEta3Pi(double mass_1, double mass_4, double cutoff);
    ///< @param mass_1 mass of the decay products
    ///< @param mass_4 decay mass
    ///< @param cutoff numerical cutoff for the dispersion integral
private:
    array::ArrayEta3Pi s_list;
    SplinedPathEta3Pi egg;
    void build_derivative() override;
    void spline_derivative() override;
};

class SplinedPathDerivativeEtap3Pi: public SplinedPathDerivative{
public:
    SplinedPathDerivativeEtap3Pi(double mass_1, double mass_4, double cutoff);
    ///< @param mass_1 mass of the decay products
    ///< @param mass_4 decay mass
    ///< @param cutoff numerical cutoff for the dispersion integral
private:
    array::ArrayEtap3Pi s_list;
    SplinedPathEtap3Pi egg;
    void build_derivative() override;
    void spline_derivative() override;
};

class SplinedPathDerivativeV3Pi: public SplinedPathDerivative{
public:
    SplinedPathDerivativeV3Pi(double mass_1, double mass_4, double cutoff);
    ///< @param mass_1 mass of the decay products
    ///< @param mass_4 decay mass
    ///< @param cutoff numerical cutoff for the dispersion integral
private:
    array::ArrayV3Pi s_list;
    SplinedPathV3Pi egg;
    void build_derivative() override;
    void spline_derivative() override;
};

class SplinedPathDerivativeX3Pi: public SplinedPathDerivative{
public:
    SplinedPathDerivativeX3Pi(double mass_1, double mass_4, double cutoff);
    ///< @param mass_1 mass of the decay products
    ///< @param mass_4 decay mass
    ///< @param cutoff numerical cutoff for the dispersion integral
private:
    array::ArrayX3Pi s_list;
    SplinedPathX3Pi egg;
    void build_derivative() override;
    void spline_derivative() override;
};

class SplinedPathDerivativeT3Pi: public SplinedPathDerivative{
public:
    SplinedPathDerivativeT3Pi(double mass_1, double mass_4, double cutoff);
    ///< @param mass_1 mass of the decay products
    ///< @param mass_4 decay mass
    ///< @param cutoff numerical cutoff for the dispersion integral
private:
    array::ArrayT3Pi s_list;
    SplinedPathT3Pi egg;
    void build_derivative() override;
    void spline_derivative() override;
};

}

#endif // SPLINED_PATH_H
