#ifndef PARTIAL_WAVE_HEADER
#define PARTIAL_WAVE_HEADER

#include "cauchy.h"
#include "constants.h"
#include "gsl_interface.h"
#include "mandelstam.h"
#include "omnes.h"
#include "type_aliases.h"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

/// @brief Partial waves

/// The partial P-wave of the scattering of an I=0, J=1, P=C=-1 particle and a
/// pion into two pions. Currently, this expression is not valid in the
/// decay region.
namespace wave {
using type_aliases::Complex;
using omnes::OmnesF;

// -- Kernel method -----------------------------------------------------------

Complex legendre(Complex z);
    ///< The lowest legendre function of the second kind.

    ///< The expression is valid in the entire complex plane with a cut
    ///< along [-1, 1].

Complex kernel(Complex mandelstam_s, double x, double pion_mass,
        double virtuality, int subtractions,
        std::optional<bool> asymptotic=std::nullopt);
    ///< The kernel used for the computation of the hat function.

    ///< @param mandelstam_s the value of Mandelstam s
    ///< @param x the value of the integration variable in the Mandelstam s
    ///< plane
    ///< @param pion_mass the mass of the pion
    ///< @param virtuality the four-momentum squared of the photon
    ///< @param subtractions the number of subtractions in the spectral
    ///< representation of the single variable function
    ///< @param asymptotic if specified and true, the kernel is computed
    ///< using an series expansion valid for certain values of `x` and
    ///< `mandelstam_s`, otherwise the ordinary expansion is used.
    ///< If unspecified, it is automatically decided which method of
    ///< computation works best.

std::vector<Complex> projected_polynomial(Complex mandelstam_s,
        double pion_mass, double virtuality, int subtractions,
        double omnes_derivative);
    ///< The projected subtraction polynomial that appears in the hat function.

/// The P-Wave.
/// @tparam Basis: should provide the same public methods as `kernel::Basis`,
/// a valid parameter is, e.g., `OmnesBasis`
template<typename Basis>
class PWaveKernel {
public:
    PWaveKernel(const Basis& basis, const gsl::Interval& sites,
            int spectral_subtractions);
        /// @param basis a set of basis functions, provided e.g. via
        /// `kernel::Basis`
        /// @param sites the sites between which the basis function is
        /// interpolated
        /// @param spectral_subtractions the number of subtractions in the
        /// spectral representation
    Complex operator()(std::size_t i, Complex mandelstam_s) const;
        /// Evaluate the partial wave basis.
    Complex evaluate(Complex mandelstam_s,
            std::vector<Complex> subtraction_constants) const;
        /// Evaluate the partial wave.
    Complex hat(std::size_t i, Complex mandelstam_s) const;
        /// Evaluate the `i`-th hat function.
private:
    Basis basis;
    int spectral_subtractions;
    gsl::Cquad integrate;
    std::vector<gsl::Interpolate> imaginary_parts;
};

template<typename Basis>
PWaveKernel<Basis>::PWaveKernel(const Basis& basis,
        const cauchy::Interval& sites, int spectral_subtractions)
    : basis{basis}, spectral_subtractions{spectral_subtractions},
    integrate{gsl::Cquad{}}, imaginary_parts{}
{
    const auto size{basis.size()};
    imaginary_parts.reserve(size);
    for (std::size_t i{0}; i < size; ++i) {
        gsl::Interval imag;
        imag.reserve(sites.size());
        std::transform(sites.cbegin(), sites.cend(), std::back_inserter(imag),
                [&basis, i](double x){ return basis(i, x).imag(); });
        imaginary_parts.push_back(gsl::Interpolate(sites, imag,
                    gsl::InterpolationMethod::linear));
    }
}

template<typename Basis>
Complex PWaveKernel<Basis>::hat(std::size_t i, Complex mandelstam_s) const
{
    const auto& inter{imaginary_parts.at(i)};
    const Complex integral{std::get<0>(cauchy::c_integrate(
            [&inter, this, mandelstam_s](double x)
            {
                const auto ker{kernel(mandelstam_s, x, basis.mass(),
                                      basis.virtuality(),
                                      spectral_subtractions)};
                return ker * inter(x);
            },
            inter.front(), inter.back(), integrate))};
    const auto poly{projected_polynomial(mandelstam_s, basis.mass(),
            basis.virtuality(), spectral_subtractions,
            basis.omnes().derivative_at_zero())};
    return poly.at(i) + integral / constants::pi();
}

template<typename Basis>
Complex PWaveKernel<Basis>::operator()(std::size_t i, Complex mandelstam_s) const
{
    return basis(i, mandelstam_s) + hat(i, mandelstam_s);
}

template<typename Basis>
Complex PWaveKernel<Basis>::evaluate(Complex mandelstam_s,
        std::vector<Complex> subtraction_constants) const
{
    if (subtraction_constants.size() != basis.size()) {
        throw std::runtime_error{"number of subtraction constants does not"
                                 " match size of basis"};
    }
    Complex result{0.0};
    for (std::size_t i = 0; i < basis.size(); ++i) {
        result += this->operator()(i, mandelstam_s) * subtraction_constants[i];
    }
    return result;
}


// -- Direct brute-force integration ------------------------------------------

/// The P-wave.
/// @tparam Basis: should provide the same public methods as `kernel::Basis`,
/// a valid parameter is, e.g., `OmnesBasis`
template<typename Basis>
class PWave {
public:
    PWave(const Basis& basis): basis{basis}, integrate{gsl::Cquad{}} {};
        /// @param basis a set of basis functions, provided e.g. via
        /// `kernel::Basis`
    Complex operator()(std::size_t i, Complex mandelstam_s) const;
        /// Evaluate the partial wave basis.
    Complex evaluate(Complex mandelstam_s,
            std::vector<Complex> subtraction_constants) const;
        /// Evaluate the partial wave.
    Complex hat(std::size_t i, Complex mandelstam_s) const;
        /// Evaluate the `i`-th hat function.
private:
    Basis basis;
    gsl::Cquad integrate;
};


template<typename Basis>
Complex PWave<Basis>::operator()(std::size_t i, Complex mandelstam_s) const
{
    return basis(i, mandelstam_s) + hat(i, mandelstam_s);
}

template<typename Basis>
Complex PWave<Basis>::evaluate(Complex mandelstam_s,
        std::vector<Complex> subtraction_constants) const
{
    if (subtraction_constants.size() != basis.size()) {
        throw std::runtime_error{"number of subtraction constants does not"
                                 " match size of basis"};
    }
    Complex result{0.0};
    for (std::size_t i = 0; i < basis.size(); ++i) {
        result += this->operator()(i, mandelstam_s) * subtraction_constants[i];
    }
    return result;
}

template<typename Basis>
Complex PWave<Basis>::hat(std::size_t i, Complex mandelstam_s) const
{
    const auto function{
        [mandelstam_s, i, this](double z)
        {
            const auto t{mandelstam::t_photon_pion(mandelstam_s, z,
                                                   basis.mass(),
                                                   basis.virtuality())};
            return (1.0 - z*z) * basis(i, t);
        }
    };
    const auto integral{cauchy::c_integrate(function, -1.0, 1.0, integrate)};
    return std::get<0>(integral) * 1.5;
}


// -- Omnes basis for partial wave projection ----------------------------------

/// As `kernel::Basis`, but with plain Omnes functions instead of full
/// Khuri-Treiman solutions.
class OmnesBasis {
public:
    OmnesBasis(OmnesF  omn, int subtractions,
               double pion_mass, double virtuality);
    ///< @param o the Omnes function
    ///< @param subtraction the number of subtractions
    ///< @param pion_mass the pion mass
    ///< @param virtuality the 'mass' squared of the I=0, J=1, P=C=-1
    ///< particle. Might take on arbitrary values (i.e. 0 and negative
    ///< values are allowed, too).
    Complex operator()(std::size_t i, Complex s) const;
    ///< @brief Evaluate the basis function with subtraction polynomial
    ///< s^`i` at `s`.
    [[nodiscard]] std::size_t size() const;
    ///< The number of basis functions (i.e., of subtractions).
    [[nodiscard]] double virtuality() const noexcept { return virt; }
    [[nodiscard]] double mass() const noexcept { return pion_mass; }
    [[nodiscard]] const omnes::OmnesF& omnes() const { return omn; }
private:
    OmnesF omn;
    int subtractions;
    double pion_mass;
    double virt;
};
} // wave
#endif // PARTIAL_WAVE_HEADER
