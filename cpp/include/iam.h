#ifndef IAM_PION_DECAY_CHIRAL_H
#define IAM_PION_DECAY_CHIRAL_H

#include "chpt.h"
#include "facilities.h"

namespace iam {
using chpt::Complex;
using chpt::t2;
using chpt::t4;
using facilities::square;

Complex iam_nlo(double mass, Complex s, double pion_decay, double l_diff);
    ///< The IAM amplitude up to NLO on the first Riemann sheet.

    ///< @param mass pion mass in physical units
    ///< @param s Mandelstam s in physical units
    ///< @param pion_decay pion decay constant in the chiral limit
    ///< in physical units
    ///< @param l_diff linear combination of LECs: l_diff := 48\pi^2 (l_2 - 2l_1)
} // iam

#endif // IAM_PION_DECAY_CHIRAL_H
