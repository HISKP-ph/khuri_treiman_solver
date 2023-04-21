#ifndef MANDELSTAM_HEADER_GUARD
#define MANDELSTAM_HEADER_GUARD

#include "phase_space.h"
#include "facilities.h"
#include "type_aliases.h"

#include <cmath>
#include <complex>

/// @brief Provide Mandelstam variables for a general four-particle process as
/// well as simplified computation in several special cases.
namespace mandelstam {
using phase_space::sigma;
using phase_space::rho;
using type_aliases::Complex;
using namespace std::complex_literals;
using facilities::square;

struct DivisionByZero : public std::exception {
    DivisionByZero() {}
    const char* what() const noexcept {return "s==0 not allowed";}
};

inline Complex kaellen(const Complex& a, const Complex& b, const Complex& c)
    /// the Kaellen function
{
    return a*a + b*b + c*c - 2.0 * (a*b + a*c + b*c);
}

inline Complex t(const Complex& s, double z, double squared_1,
        double squared_2, double squared_3, double squared_4)
    /// the Mandelstam variable t in the CMS
{
    if (s==0.0)
        throw DivisionByZero{};
    double sum{squared_1+squared_2+squared_3+squared_4};
    double delta_1{squared_1-squared_2};
    double delta_2{squared_3-squared_4};
    Complex kaellen_1{kaellen(s,squared_1,squared_2)};
    Complex kaellen_2{kaellen(s,squared_3,squared_4)};
    return (sum - s - (delta_1*delta_2 - z*std::sqrt(kaellen_1*kaellen_2))/s)
        / 2.0;
}

inline Complex u(const Complex& s, double z, double squared_1,
        double squared_2, double squared_3, double squared_4)
    /// the Mandelstam variable u in the CMS
{
    return t(s,-z,squared_1,squared_2,squared_4,squared_3);
}

inline double s_greater(double pion_mass, double virtuality)
    /// @brief the upper bound of the region in which t is complex for
    /// photon+pion->pion+pion
{
    if (virtuality<0.0)
        throw std::invalid_argument{"virtuality needs to be non-negative"};
    double temp{std::sqrt(virtuality) + pion_mass};
    return temp*temp;
}

inline double s_smaller(double pion_mass, double virtuality)
    /// @brief the upper bound of the region in which t is complex for
    /// photon+pion->pion+pion
{
    return s_greater(-pion_mass,virtuality);
}

inline Complex a_photon_pion(Complex s, double pion_mass,
        double virtuality)
{
    return (3.0*pion_mass*pion_mass + virtuality - s) / 2.0;
}

inline Complex b_photon_pion(Complex s, double pion_mass,
        double virtuality, bool schwarz=false)
{
    const auto phase_space{schwarz ? sigma(pion_mass, s) : rho(pion_mass, s)};
    if (virtuality <= 0.0)
        return 0.5 * phase_space
            * std::sqrt(kaellen(s,virtuality,pion_mass*pion_mass));

    Complex sqrt_1{std::sqrt(s - s_greater(pion_mass,virtuality))};
    Complex sqrt_2{std::sqrt(s - s_smaller(pion_mass,virtuality))};
    return 0.5*phase_space*sqrt_1*sqrt_2;
}

inline Complex t_photon_pion(const Complex& s, double z,
        double pion_mass, double virtuality, bool schwarz=false)
    /// the Mandelstam variable t for photon+pion->pion+pion in the CMS
{
    return a_photon_pion(s,pion_mass,virtuality)
        + z*b_photon_pion(s,pion_mass,virtuality,schwarz);
}

inline Complex t_photon_pion_min(const Complex&s, double pion_mass,
        double virtuality)
    /// @brief the Mandelstam variable t for photon+pion->pion+pion in the CMS
    /// evaluated at z=-1.0
{
    return t_photon_pion(s,-1.0,pion_mass,virtuality);
}

inline Complex t_photon_pion_max(const Complex&s, double pion_mass,
        double virtuality)
    /// @brief the Mandelstam variable t for photon+pion->pion+pion in the CMS
    /// evaluated at z=1.0
{
    return t_photon_pion(s,1.0,pion_mass,virtuality);
}

/// @brief The characteristics of the singular region where Mandelstam t hits
/// the branch point at the two-pion threshold.
///
/// The region is contained in a square in the complex plane whose left boundary
/// is a vertical line at `left`, the right one is a vertical line at `right`,
/// the upper/lower one a horizontal line at +/- `imaginary_radius`.
class Critical {
public:
    constexpr Critical(double pion_mass, double virtuality)
        : pion_mass{pion_mass}, virtuality{virtuality} {}

    double imaginary_radius() const
        /// @brief Return an upper bound for the maximal imaginary value of the
        /// region.
    {
        return std::abs(virtuality - 8.0*square(pion_mass)) / 3.0;
    }

    constexpr double left() const
    {
        return 0.5 * (virtuality - square(pion_mass));
    }

    constexpr double right() const
    {
        return virtuality - 5.0*square(pion_mass);
    }
private:
    double pion_mass;
    double virtuality;
};
} // mandelstam

#endif // MANDELSTAM_HEADER_GUARD
