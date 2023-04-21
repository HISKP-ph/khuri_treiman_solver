// Evaluate the angular averages along a complex integration contour
// that avoids the cut.
// Additionally match an expansion of the angular average close to scattering thresholds.



#ifndef PATH_ETA_PI_PI_H
#define PATH_ETA_PI_PI_H

#include "array.h"
#include "type_aliases.h"
#include "cauchy.h"

#include <complex>
#include <cmath>


namespace path_eta_pi_pi{

using type_aliases::Complex;
using namespace std::complex_literals;
using namespace constants;
using namespace enums;

/// Class for Path deformation for a final state with two different masses.
class Path{
public:
    Path(double mass_1, double mass_2, double mass_3, double mass_4);
        ///< @param mass_1, mass_2 charged pion masses
        ///< @param mass_3 eta mass
        ///< @param mass_4 eta' mass
    

    
    /// Contour in terms of a curve-parameter y with simplified interface.
    Complex operator()(double y, Setting evaluation)const;
    
    /// Derivative of the contour,
    Complex derivative(double y, Setting evaluation)const;
    
    /// Functions needed to determine the integration limits.
    double integration_limit(double s, TypeOfAverage average,
                              PinocchioEndpoints endpoints)const;
    
    
    /// Express mandelstam t in terms of curve parameter ya.
    double t_to_ya(double ya)const;
    /// Express mandelstam t in terms of curve parameter yb.
    double t_to_yb(double yb)const;
    
    
    /// This function is used to evaluate the limit for the integration in region III.
    double ya(double t)const;
    /// This function is used to evaluate the limit for the integration in region III.
    double yb(double t)const;
    
    
    double delta;
    
private:
    double mass_1, mass_2, mass_3, mass_4;
    double s_I, s_III, s_IV, s_0;
    double t_I, t_III, t_IV, t_0;
    
    /// @brief Use endpoints of Pinocchio-contour to parameterize the path in terms of real Mandelstam s.
    ///
    /// Caution: Cannot directly integrate over it, due to end-point-singularities.
    /// Therefore `contour` will be used.
    Complex s_plus_minus(double s, TypeOfAverage average, PinocchioEndpoints endpoints)const;
    
    /// Contour in terms of a curve-parameter y.
    Complex contour(double y, TypeOfAverage average,
                       PinocchioEndpoints endpoints, IntegrationRegion region)const;
    
    /// Kinematic function
    double kappa_pi(double s)const;
    /// Kinematic function
    double kappa_eta(double t)const;
    /// Källen function
    Complex kallen(double s, double mass1, double mass2)const;//Källén function

};

} // path_eta_pi_pi

#endif // PATH_ETA_PI_PI_H
