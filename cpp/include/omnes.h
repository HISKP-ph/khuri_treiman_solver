#ifndef OMNES_FUNCTION_H
#define OMNES_FUNCTION_H

#include "cauchy.h"
#include "constants.h"
#include "gsl_interface.h"
#include "helpers.h"
#include "phase_space.h"
#include "type_aliases.h"

#include <cmath>
#include <complex>

/// The Omnes function of an arbitrary phase.

/// The way of thinking about the Omnes function is that it provides a
/// different function for different phases.
/// Hence, a class `Omnes` is provided. Each instance of the class provides
/// a different Omnes function, i.e. one for a different specific phase.
namespace omnes {
using namespace std::complex_literals;
using type_aliases::Complex;
using type_aliases::CFunction;
using helpers::hits_threshold;

/// The Omnes function for arbitrary phases and thresholds.
template<typename Integrate=gsl::Cquad>
class Omnes;

using OmnesF = Omnes<>;

template<typename Integrate>
class Omnes {
public:
    Omnes(const gsl::Function& phase, double threshold,
            double minimal_distance, gsl::Settings config=gsl::Settings{});
        ///< @param phase The phase of the Omnes function above its branch cut.
        ///< @param threshold The start of the branch cut of the Omnes function.
        ///< @param minimal_distance Half the width of a band
        ///< around the cut. For arguments of the Omnes function in this band,
        ///< a different prescription is used for the evaluation of the Omnes
        ///< function to take care of the singularity in the integral.
        ///< @param config The settings for the integration routine.

    Omnes(const gsl::Function& phase, double threshold, double constant,
            double cut, double minimal_distance,
            gsl::Settings config=gsl::Settings{});
        ///< @param phase The phase of the Omnes function in
        ///< [`threshold`,`cut`].
        ///< @param threshold The start of the branch cut of the Omnes function.
        ///< @param constant The phase of the Omnes function is set equal to
        ///< `constant` along the real line above `cut`.
        ///< @param cut Cf. `phase` and `constant`
        ///< @param minimal_distance Half the width of a band
        ///< around the cut. For arguments of the Omnes function in this band,
        ///< a different prescription is used for the evaluation of the Omnes
        ///< function to take care of the singularity in the integral.
        ///< @param config The settings for the integration routine.

    Complex operator()(Complex s) const;
        ///< Evaluate the Omnes function at `s`.

    double derivative_at_zero() const noexcept {return derivative;}
        ///< Return the derivative of the Omnes function at the origin.

    double branch_point() const noexcept {return threshold;}
        ///< Return the branch point.
private:
    const gsl::Function phase_below; // phase below `cut`
    const double constant;
    const double threshold;
    const double cut;
    const double minimal_distance;
    const Integrate integrate;
    const double derivative;

    Complex upper(const Complex& s) const;
        // Evaluate the Omnes function in the upper half of the complex plane
    bool hits_cut(const Complex& s) const;
        // Return true if `s` is in the region around the branch cut, false
        // otherwise.
    Complex threshold_presciption(double s) const;
        // Calculate the Omnes function if `s` is close to `threshold`.
    Complex
            ordinary_prescription(const Complex& s) const;
        // Calculate the Omnes function if `s` is not close to the branch
        // cut.
    Complex cut_prescription(double s) const;
        // Calculate the Omnes function if `s` is close to the branch cut.
    double phase(double s) const;
        // Calculate the phase of the Omnes function along the branch cut.
    double abs_cut(double s) const;
        // Calculate the absolute value of the Omnes function along the branch
        // cut.
};

inline double derivative_0(const gsl::Function& phase, double threshold,
       double cut, double constant, const gsl::Integration& integrate)
    // Return the derivative of the Omnes function at s=0. The parameters are
    // the same as the ones with the same name in the constructor of
    // `class Omnes`.
{
    double first{integrate([&phase](double x){return phase(x)/(x*x);},
            threshold,cut).first};
    double second{constant/cut};
    return (first + second)/constants::pi();
}

template<typename T>
Omnes<T>::Omnes(const gsl::Function& phase, double threshold,
        double minimal_distance, gsl::Settings config)
: phase_below{phase},
    constant{0.0}, // value of the `constant` is irrelevant if `cut` is infinity
    threshold{threshold}, cut{std::numeric_limits<double>::infinity()},
    minimal_distance{minimal_distance},
    integrate{config},
    derivative{derivative_0(phase,threshold,cut,constant,integrate)}
{
}

template<typename T>
Omnes<T>::Omnes(const gsl::Function& phase, double threshold, double constant,
        double cut, double minimal_distance, gsl::Settings config)
: phase_below{phase}, constant{constant}, threshold{threshold}, cut{cut},
    minimal_distance{minimal_distance},
    integrate{config},
    derivative{derivative_0(phase,threshold,cut,constant,integrate)}
{
}

template<typename T>
Complex Omnes<T>::operator()(Complex s) const
{
    // Apply the Schwartz reflection principle.
    if (s.imag()<0)
        return std::conj(upper(std::conj(s)));
    else
        return upper(s);
}

template<typename T>
Complex Omnes<T>::upper(const Complex& s) const
{
    if (hits_threshold(threshold, s, minimal_distance))
        return threshold_presciption(s.real());
    if (hits_cut(s))
        return cut_prescription(s.real());
    else
        return ordinary_prescription(s);
}

template<typename T>
bool Omnes<T>::hits_cut(const Complex& s) const
{
    return s.real()>=threshold && std::abs(s.imag())<=minimal_distance;
}

template<typename T>
Complex Omnes<T>::threshold_presciption(double s) const
{
    // average
    return (cut_prescription(threshold+minimal_distance)
            +ordinary_prescription(threshold-minimal_distance)) / 2.0;
}

template<typename T>
Complex Omnes<T>::ordinary_prescription(
        const Complex& s) const
{
    Complex above_cut{std::log(1.0-s/cut)};
    auto integral{std::get<0>(cauchy::c_integrate(
                [&s,this](double z){return phase_below(z)/(z*(z-s));},
                threshold,cut,integrate))};
    return std::exp((s*integral-constant*above_cut)/constants::pi());
}

template<typename T>
Complex Omnes<T>::cut_prescription(double s) const
{
    const Complex arg{phase(s) * 1.0i};
    return abs_cut(s) * std::exp(arg);
}

template<typename T>
double Omnes<T>::phase(double s) const
{
    if (s<cut)
        return phase_below(s);
    else
        return constant;
}

inline double abs_helper(double s, double value)
    // Simplify the calculation of the absolute value of the Omnes
    // function along the cut.
{
    double temp{1.0 - s/value};
    return std::log(std::abs(1.0/temp));
}

template<typename T>
double Omnes<T>::abs_cut(double s) const
{
    double phase_at_s{phase(s)};
    auto integral{integrate(
                [&s,&phase_at_s,this](double z)
                    {return (phase_below(z)-phase_at_s)/(z*(z-s));},
                threshold,cut).first};
    double a{s<cut ? constant-phase_at_s : 0.0};
    return std::exp((s*integral + a*abs_helper(s,cut)
            + phase_at_s*abs_helper(s,threshold))/constants::pi());
}

template<typename T>
Complex second_sheet(const Omnes<T>& o, const CFunction& amplitude,
        const Complex& s)
    /// Evaluate the Omnes function on the second Riemann sheet.
    /// @param o The Omnes function.
    /// @param amplitude The two-particle scattering amplitude associated with
    /// the phase of the Omnes function.
    /// @param s The value of Mandelstam s.
{
    const double pion_mass{std::sqrt(o.branch_point() / 4.0)};
    const Complex z{2.0i};
    return o(s) / (1.0 + z * phase_space::rho(pion_mass, s) * amplitude(s));
}
} // omnes

#endif // OMNES_FUNCTION_H
