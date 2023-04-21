#ifndef HELPERS_HEADER_GUARD
#define HELPERS_HEADER_GUARD

#include "facilities.h"
#include "type_aliases.h"

/// Useful helpers that are needed in several places and that are specific
/// to Khuri Treiman equations/dispersion relations.
namespace helpers {
using type_aliases::Complex;
using facilities::square;

inline constexpr double threshold(double pion_mass)
    /// Return the two-pion threshold.
{
    return 4.0*square(pion_mass);
}

inline bool hits_threshold(double threshold, const Complex& s,
                           double minimal_distance)
    /// Return true if `s` hits the two-pion threshold.
{
    return std::abs(s-threshold) < minimal_distance;
}

inline bool hits_threshold_m(double pion_mass, const Complex& s,
                             double minimal_distance)
    /// Return true if `s` hits the two-pion threshold.
{
    return hits_threshold(threshold(pion_mass), s, minimal_distance);
}
} // helpers
#endif // HELPERS_HEADER_GUARD
