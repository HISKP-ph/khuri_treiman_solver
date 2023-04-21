#ifndef CHPT_AMPLITUDES_H
#define CHPT_AMPLITUDES_H

#include "constants.h"
#include "facilities.h"
#include "phase_space.h"
#include "type_aliases.h"

#include <complex>

/// @brief \f$I=J=1\f$ \f$\pi\pi\to\pi\pi\f$ ChPT-partial-wave amplitudes up to
/// NLO in terms of the pion decay constant in the chiral limt.
namespace chpt {
using constants::pi;
using facilities::square;
using phase_space::rho;
using phase_space::sigma;
using type_aliases::Complex;
using std::pow;
using std::sqrt;
using std::log;

inline Complex t2(double mass, Complex s, double pion_decay)
    /// The CHPT LO amplitude.

    /// @param mass pion mass in physical units
    /// @param s Mandelstam s in physical units
    /// @param pion_decay pion decay constant (either in chiral limit or not)
    /// in physical units
{
    return (s-4.0*square(mass)) / (96.0*square(pion_decay)*pi());
}

Complex t4(double mass, Complex s, double pion_decay, double l_diff);
    ///< The CHPT NLO amplitude.

    ///< @param mass pion mass in physical units
    ///< @param s Mandelstam s in physical units
    ///< @param pion_decay pion decay constant in the chiral limit
    ///< in physical units
    ///< @param l_diff linear combination of LECs:
    ///< l_diff := 48\pi^2 (l_2 - 2l_1)

} // chpt

#endif // CHPT_AMPLITUDES_H
