#include "iam.h"

namespace iam {
Complex iam_nlo(double mass, Complex s, double pion_decay, double l_diff)
{
    const Complex lo = t2(mass,s,pion_decay);
    const Complex nlo = t4(mass,s,pion_decay,l_diff);
    return square(lo) / (lo - nlo);
}
} // iam
