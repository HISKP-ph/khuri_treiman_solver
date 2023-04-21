#ifndef ASYMPTOTIC
#define ASYMPTOTIC

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <utility>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include "constants.h"
#include "gsl_interface.h"

namespace asymptotic{
using constants::pi;
using gsl::Interpolate;
using gsl::Function;

/// @brief Return phase that approaches a constant at high energies.
/// 
/// The decorated phase is used below `matching_point`, above which `limit`
/// is approached eventually. The transition is one time differentiable.
class Asymptotic1{
private:
    Function phase;
    /// Point where the phase is matched.
    double matching_point;
    /// Limit that is approached.
    double limit;
    /// Phase at the matching point.
    double value;
    /// Derivative of the phase at the matching point.
    double d_value;
    /// Parameter1 from matching value and d_value at the matching point.
    double param1;
    /// Parameter2 from matching value and d_value at the matching point.
    double param2;
public:
    Asymptotic1(Function phasein, double matching_pointin, double limitin=pi());

    double continued_phase(double s)
    {
        if (s < matching_point){
            return phase(s);
        }
        else{
            return limit - param1 / (param2 + s / matching_point);
        }
    }

    double operator()(double s)
    {
        return continued_phase(s);
    }
};

/// @brief Return phase that approaches a constant at high energies for splines.
///
/// The decorated phase is used below `matching_point`, above which `limit`
/// is approached eventually. The transition is one times differentiable.
class Asymptotic1s{
private:
    Interpolate phase;
    /// Point where the phase is matched.
    double matching_point;
    /// Limit that is approached.
    double limit;
    /// Phase at the matching point.
    double value;
    /// Derivative of the phase at the matching point.
    double d_value;
    /// Parameter1 from matching value and d_value at the matching point.
    double param1;
    /// Parameter2 from matching value and d_value at the matching point.
    double param2;
public:
    Asymptotic1s() = default;
    Asymptotic1s(Interpolate phasein, double matching_pointin, double limitin=pi());

    double continued_phase(double s) const
    {
        if (s < matching_point){
            return phase(s);
        }
        else{
            return limit - param1 / (param2 + s / matching_point);
        }
    }

    double operator()(double s) const
    {
        return continued_phase(s);
    }
};

/// @brief Return phase that equals `phase` below `lower` and `limit` above `upper`.
/// Between `lower` and `upper` there is a in [`lower`, `upper`] one
/// times continuously differentiable smooth connection.
///
/// The matching is achieved using a third-order polynomial, which
/// requires to use also the phase as a part of the function in
/// [`lower`, `upper`].
class Asymptotic2{
private:
    Function phase;
    /// matching point lower
    double lower;
    /// matching point upper
    double upper;
    /// Limit that is approached.
    double limit;
public:
    Asymptotic2(Function phasein, double lowerin, double upperin, double limitin=pi());

    double connect(double s);

    double continued_phase(double s)
    {
        if (s < lower){
            return phase(s);
        }
        else if (s < upper){
            double c = connect(s);
            return (1.0 - c) * phase(s) + c * limit;
        }
        else{
            return limit;
        }
    }

    double operator()(double s)
    {
        return continued_phase(s);
    }
};

}

#endif // ASYMPTOTIC