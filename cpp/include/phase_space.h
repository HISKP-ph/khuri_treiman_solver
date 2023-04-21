#ifndef PHASE_SPACE_H
#define PHASE_SPACE_H

#include "facilities.h"
#include "type_aliases.h"

#include <complex>

/// @brief Different versions of the two-particle phase space.
namespace phase_space {
using type_aliases::Complex;
using facilities::square;

inline int signum_im(const Complex& x)
    /// Signum of the imaginary part of a number.
{
    return imag(x) < 0.0 ? -1 : 1;
}

inline Complex alt_sqrt(const Complex& x)
    /// @brief Square root function with cut on positive real axis.
    /// (std::sqrt has cut on negative real axis.)
{
    return static_cast<double>(signum_im(x))*std::sqrt(x);
}

inline Complex rho(double mass, const Complex& s)
    /// The two-body phase space (cuts along [4mass^2,\infty) and (-\infty,0]).

    /// This one does not obey the Schwarz reflection principle.
{
    return alt_sqrt(1.0 - 4.0*square(mass)/s);
}

inline Complex sigma(double mass, const Complex& s)
    /// The two-body phase space (cut along [0,4mass^2]).

    /// This one obeys the Schwarz reflection principle.
{
    return std::sqrt(1.0 - 4.0*square(mass)/s);
}
} // phase_space

#endif // PHASE_SPACE_H
