#include "constants.h"
#include "partial_wave.h"
#include "phase_space.h"

#include <sstream>
#include <stdexcept>

namespace wave {

using mandelstam::a_photon_pion;
using mandelstam::b_photon_pion;

Complex legendre(Complex z)
{
    if (z.imag() == 0.0) {
        // ensure that the branch cut is approached from above
        z += Complex{0.0, 1e-20};
    }
    return 0.5 * (std::log((1.0 + z) / (1.0 - z))
                  - Complex{0.0, constants::pi() * phase_space::signum_im(z)});
}

Complex kernel_core(Complex mandelstam_s, double x, double pion_mass,
        double virtuality)
{
    const auto a{a_photon_pion(mandelstam_s, pion_mass, virtuality)};
    const auto b{b_photon_pion(mandelstam_s, pion_mass, virtuality, true)};
    const auto xi{(x - a) / b};
    return 3.0 * ((1.0 - xi * xi) * legendre(xi) + xi) / b;
}

Complex kernel_core_asymptotic(Complex mandelstam_s, double x,
        double pion_mass, double virtuality)
{
    const auto a{a_photon_pion(mandelstam_s, pion_mass, virtuality)};
    const auto b{b_photon_pion(mandelstam_s, pion_mass, virtuality, true)};
    const auto diff{x - a};
    const auto factor{b / diff};
    const auto factor2{factor * factor};
    Complex res{0.0};
    Complex power{1.0};
    for (int i{0}; std::abs(power) > 1e-16; ++i) {
       res += power * (2.0 / (2.0 * i + 1.0) / (2.0 * i + 3.0));
       power *= factor2;
    }
    return 3.0 * res / diff;
}

bool converges(Complex mandelstam_s, double x, double pion_mass,
        double virtuality)
    // check if values are inside radius of convergence of asymptotic kernel
{
    const auto a{a_photon_pion(mandelstam_s, pion_mass, virtuality)};
    const auto b{b_photon_pion(mandelstam_s, pion_mass, virtuality, true)};
    return std::abs(x - a) > std::abs(b);
}

bool use_asymptotic(Complex mandelstam_s, double x, double pion_mass,
        double virtuality, std::optional<bool> asymptotic)
    // decide if to use asymptotic prescription of kernel or ordinary one
{
    return asymptotic
        ? *asymptotic
        : converges(mandelstam_s, x, pion_mass, virtuality);
}

Complex kernel(Complex mandelstam_s, double x, double pion_mass,
        double virtuality, int subtractions, std::optional<bool> asymptotic)
{

    const auto core{
        use_asymptotic(mandelstam_s, x, pion_mass, virtuality, asymptotic)
        ?  kernel_core_asymptotic(mandelstam_s, x, pion_mass, virtuality)
        : kernel_core(mandelstam_s, x, pion_mass, virtuality)};
    switch (subtractions) {
        case 1:
            return core - 2.0 / x;
        case 2: {
            const auto a{a_photon_pion(mandelstam_s, pion_mass, virtuality)};
            return core - 2.0 / x - 2.0 * a / x / x;
            }
        default:
            std::ostringstream message;
            message << "kernel not implemented for " << subtractions
                << " subtractions";
            throw std::runtime_error{message.str()};
    }
}

std::vector<Complex> projected_polynomial(Complex mandelstam_s,
        double pion_mass, double virtuality, int subtractions,
        double omnes_derivative)
{
    switch (subtractions) {
        case 1:
            return {2.0};
        case 2: {
            const auto a{a_photon_pion(mandelstam_s, pion_mass, virtuality)};
            return {2.0 + omnes_derivative * 2.0 * a, 2.0 * a};
            }
        default:
            std::ostringstream message;
            message << "projected polynomial not implemented for "
                << subtractions << " subtractions";
            throw std::runtime_error{message.str()};
    }
}

OmnesBasis::OmnesBasis(OmnesF  omn, int subtractions,
                       double pion_mass, double virtuality)
        : omn{std::move(omn)}, subtractions{subtractions}, pion_mass{pion_mass},
          virt{virtuality}
{
    if (subtractions < 1) {
        throw std::runtime_error{"Number of subtractions smaller than 1."};
    }
}

Complex OmnesBasis::operator()(std::size_t i, Complex s) const
{
    return std::pow(s, i) * omn(s);
}

std::size_t OmnesBasis::size() const
{
    // The class invariant assures that `subtractions` is non-negative.
    using V = std::make_unsigned_t<decltype(subtractions)>;
    return static_cast<V>(subtractions);
}
} // wave
