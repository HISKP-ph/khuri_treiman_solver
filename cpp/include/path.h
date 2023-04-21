#ifndef PATH_H
#define PATH_H

#include "cauchy.h"
#include "constants.h"
#include "facilities.h"
#include "mandelstam.h"
#include "type_aliases.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <exception>
#include <limits>
#include <vector>

/// The parametrization of the Pinocchio path.

/// A general interface of curves in the complex plane with a paremetrization
/// that can be piecewise targeted at parametrizations of the Pinocchio path
/// as well as an explicit parametrization of this path are provided.
/// The general interface is the class `Path`, the explicit parametrization
/// the class `Pinocchio`. Currently, the latter is restricted to decays
/// into three particles of equal mass.
namespace path {
using facilities::square;
using type_aliases::Complex;
using type_aliases::Curve;

/// An interval along the real line.
struct Interval {
    double lower;
    double upper;

    bool contains(double x)
        /// Check if `x` is contained in the interval.
    {
        return lower <= x && x <= upper;
    }
};

/// @brief A general interface of a curve in the complex plane with a mapping
/// from Mandelstam s to the curve parameter.
class Path {
public:
    virtual Complex operator()(double curve_parameter) const = 0;
        ///< Evaluate the curve.

    Complex derivative(double curve_parameter, double step_size=1e-8) const;
        ///< Evaluate the derivative of the curve.

    virtual Interval curve_parameter(double mandelstam_s) const = 0;
        ///< @brief Lower and upper values corersponding to `mandelstam_s`.
        ///<
        ///< In this way, each value of `mandelstam_s` gets mapped to a
        ///< subset of the path, connecting
        ///< `Path::operator()(curve_parameter(mandelstam_s).lower)`
        ///< with
        ///< `Path::operator()(curve_parameter(mandelstam_s).upper)`.

    virtual Interval domain() const = 0;
        ///< The domain of definition.

    virtual ~Path() = default;
};

/// @brief An explicit parametrization of the egg-like part of the Pinocchio
/// path.
///
/// Usually, this does not need to be used directly, instead `Pinocchio`
/// should be preferred.
class PolarEgg {
public:
    PolarEgg(double decay_mass, double pion_mass);

    Complex operator()(double phi) const;

    Interval phi(double mandelstam_variable) const;
private:
    double radius(double cosine) const;
    double sum() const;
    double squared_diff() const;
    double radius_helper(double cosine, int flag) const;
    double angle(double cosine, int flag) const;
private:
    double decay_mass;
    double pion_mass;
};

/// An explicit (piecewise) parametrization of the Pinocchio path. 
class Pinocchio: public Path {
public:
    Pinocchio(double decay_mass, double pion_mass,
              double cut=std::numeric_limits<double>::infinity(),
              double epsilon=0.0);
        ///< @param decay_mass the mass of the decaying particle
        ///< @param pion_mass the mass of the pions (more generally speaking,
        ///< the mass of a single particle in the final state)
        ///< @param cut the parametrization is valid for value of Mandelstam s
        ///< below this cutoff
        ///< @param epsilon the integration path is shifted away from the cut
        ///< by this amount
    Complex operator()(double curve_parameter) const override;
    Interval curve_parameter(double mandelstam_s) const override;
    Interval domain() const override;
    std::vector<double> parameter_boundaries() const;
        ///< @brief The values of the curve parameter at which one segment of
        ///< the piecewise parametrization ends and another starts.
    std::vector<double> mandelstam_boundaries() const;
        ///< @brief The mandelstam values corresponding to
        ///< `parameter_boundaries`.
        ///<
        ///< That is, `curve_parameter` evaluated at an element of
        ///< `mandelstam_boundaries` equals the corresponding element of
        ///< `parameter_boundaries`.
private:
    double one() const { return 4.0 * square(pion_mass); }
    double two() const { return 0.5 * (square(decay_mass)-square(pion_mass)); }
    double three() const { return square(decay_mass - pion_mass); }
    double four() const { return square(decay_mass + pion_mass); }

    double lower(double mandelstam_s) const;
        /// The curve parameter corresponding to the lower end.
    double upper(double mandelstam_s) const;
        /// The curve parameter corresponding to the upper end.
    double shift() const { return 2.0 * one(); }
        /// In the upper half of the path, the curve parameter is shifted by
        /// this value.
private:
    double decay_mass;
    double pion_mass;
    double cut;
    double epsilon;
    PolarEgg egg;
};

/// This error is thrown whenever an explicit parametrization of a path fails.
class ParametrizationError: public std::exception {
public:
    explicit ParametrizationError(const char* message) : message(message) {}
    explicit ParametrizationError(const std::string& message)
        : message{message} {}

    const char* what() const noexcept override
    {
        return message.c_str();
    }
private:
    std::string message;
};
} // path

#endif // PATH_H
