"""
=======================================================
Parametrization of pi+pi -> pi+pi scattering amplitudes
=======================================================

Parametrizations of the pi+pi -> pi+pi scattering amplitudes and phases of the
s0-wave and the p-wave from

    <arXiv:1907.13162 by J. R. Peláez, A. Rodas, J. Ruiz de Elvira>.

Along the real axis the parametrizations are valid up to 2 GeV.
The p-wave amplitude below 1.12 GeV is holomorphic and hence valid in the
complex plane, the same holds true for the s-wave amplitude below 1.4 GeV.
However, above the 2-Kaon-threshold the s-wave has additional sheets, so the
first and second are not longer continuously connected.

The notation
from the aforementioned paper is used in the implementation. The values for
the non-dimensionless fit parameters/default arguments are given in GeV, so
Mandelstam s should accordingly be given in units of (GeV)^2.
"""
import numpy as np
import numpy.polynomial.chebyshev as cheb
from numpy.polynomial.polynomial import Polynomial
from numpy.lib.scimath import sqrt as csqrt, log as clog
from scipy.misc import derivative

from khuri.amplitude import from_cot, second_sheet
from khuri.phase_space import rho


# general / helpers ###########################################################


PION_MASS = 0.13957  # in GeV
KAON_MASS = 0.496  # in GeV
S_MATCHING = 1.4**2  # in GeV


def arccot2(x):
    """The arccot for real arguments."""
    return np.arctan2(1.0, np.real(x))


@np.vectorize
def arctan2_alt(x, y):
    """The arctan2 with values in [0, 2pi]."""
    arc = np.arctan2(x, y)
    if arc < 0:
        arc += 2.0 * np.pi
    return arc


def chebyshev_derivative(number, x):
    """Evaluate the derivative of `number`-th Chebyshev polynomial at `x`."""
    coefficients = np.zeros(number + 1)
    coefficients[number] = 1
    return cheb.chebval(x, cheb.chebder(coefficients))


def momentum(s):
    return csqrt(s) * rho(PION_MASS, s) / 2.0


def omega(s, s_0, alpha):
    """A conformal variable.

    Parameters
    ----------
    s
        the Mandelstam variable `s`
    s_0
        usually non-negative. Above `s_0`, the conformal variable is complex.
    alpha
        the center of the conformal expansion
    """
    sqrt_s = csqrt(s)
    diff = alpha * csqrt(s_0 - s)
    return (sqrt_s - diff) / (sqrt_s + diff)


def omega_2(s, s_m=S_MATCHING, constant=2.0):
    """Variable for Chebyshev polynomials used in the high-energy region."""
    sqrt_m = csqrt(s_m)
    return 2.0 * (csqrt(s) - sqrt_m) / (constant - sqrt_m) - 1.0


# s-wave ######################################################################


S_WAVE_F_0_POLE = (0.996 - 0.025j)**2

S_WAVE_CONFORMAL_COEFFICIENTS_1 = (12.2, -0.9, 15.9, -5.7, -22.5, 6.9)
S_WAVE_F_0_COEFFICIENTS_1 = (1.29, -1.08, -0.043, -0.068)
S_WAVE_ADLER_1 = 0.137

S_WAVE_CONFORMAL_COEFFICIENTS_2 = (12.2, -1.2, 15.5, -6.0, -21.4, 6.3)
S_WAVE_F_0_COEFFICIENTS_2 = (1.22, -1.16, -0.01, -0.075)
S_WAVE_ADLER_2 = 0.135


def omega_1(s):
    return omega_2(s, s_m=4.0*KAON_MASS**2, constant=1.5)


def loop_function(s, mass):
    phase_space = rho(mass, s)
    return ((2.0 + phase_space * clog((phase_space - 1) / (phase_space + 1)))
            / np.pi)


def generate_conformal_amplitude(conformal_coeff, adler):
    conformal_polynomial = Polynomial(conformal_coeff)
    s_0 = 4.0 * KAON_MASS**2
    a2 = adler**2

    def conformal_amplitude(s):
        conf = conformal_polynomial(omega(s, s_0=s_0, alpha=1))
        prefactor = PION_MASS**2 / (s - 0.5 * a2)
        psi = prefactor * (a2 / PION_MASS / csqrt(s) + conf)
        return 1.0 / (psi - 1j * rho(PION_MASS, s))

    return conformal_amplitude


def generate_f_0_amplitude(coefficients, pole):
    chebyshev = cheb.chebval(omega_1(pole), coefficients)
    loop_kaon = loop_function(pole, KAON_MASS)
    loop_pion = loop_function(pole, PION_MASS)
    phase_space = rho(PION_MASS, pole)
    denominator = (loop_pion * pole).imag - 2.0 * (phase_space * pole).real

    cheb_kaon = chebyshev * loop_kaon

    g_parameter = -(pole + cheb_kaon).imag / denominator
    m_parameter = (cheb_kaon.real
                   + ((loop_pion.imag - 2.0 * phase_space.real) * abs(pole)**2
                      - cheb_kaon.imag
                      * ((pole * loop_pion).real
                         + 2.0 * (pole * phase_space).imag))
                   / denominator)

    def amplitude(s):
        chebyshev = cheb.chebval(omega_1(s), coefficients)
        loop_1 = loop_function(s, PION_MASS)
        loop_2 = loop_function(s, KAON_MASS)
        return (s * g_parameter
                / (m_parameter - s - loop_1 * s * g_parameter
                   - loop_2 * chebyshev))

    return amplitude


def end_of_increase(values):
    """The position of the first value that is smaller than its predecessor."""
    for i, (previous, current) in enumerate(zip(values, values[1:])):
        if previous > current:
            return i
    return None


def generate_s_wave(conformal_coeff, adler, f_0_coefficients, f_0_pole):
    conformal_amplitude = generate_conformal_amplitude(conformal_coeff, adler)
    f_0_amplitude = generate_f_0_amplitude(f_0_coefficients, f_0_pole)

    def amplitude(s):
        conf = conformal_amplitude(s)
        f_0 = f_0_amplitude(s)
        return conf + f_0 + 2j * rho(PION_MASS, s) * conf * f_0

    def phase(s):
        amp = amplitude(s)
        space = rho(PION_MASS, s)
        numerator = 2.0 * space * amp.real
        denominator = 1.0 - 2.0 * space * amp.imag
        return 0.5 * arctan2_alt(numerator, denominator)

    energies_squared = np.linspace(0.98**2, 1.0**2, 10000)
    phases = phase(energies_squared)
    boundary = energies_squared[end_of_increase(phases)]

    if boundary is None:
        def increasing_phase(s):
            s = s.real
            return phase(s)
    else:
        def increasing_phase(s):
            s = np.asarray(s).real
            shift = np.pi
            return np.piecewise(s,
                                [s < boundary, s >= boundary],
                                [phase, lambda x: phase(x) + shift])

    def inelasticity(s):
        amp = amplitude(s)
        space = rho(PION_MASS, s)
        return np.sqrt((2.0 * space * amp.real)**2
                       + (1.0 - 2.0 * space * amp.imag)**2)

    return (increasing_phase,
            amplitude,
            second_sheet(PION_MASS)(amplitude),
            inelasticity)


(s_wave_phase,
 s_wave,
 s_wave_2,
 s_wave_inelasticity) = generate_s_wave(S_WAVE_CONFORMAL_COEFFICIENTS_2,
                                        S_WAVE_ADLER_2,
                                        S_WAVE_F_0_COEFFICIENTS_2,
                                        S_WAVE_F_0_POLE)


# p-wave ######################################################################


P_WAVE_RHO_MASS = 0.7752  # in GeV

P_WAVE_CONFORMAL_COEFFICIENTS_1 = (0.97, 0.12, -0.18, 0.40, 1.65)
P_WAVE_D_0_1 = 11.1
P_WAVE_D_1_1 = 5.6
P_WAVE_K_0_1 = 0.35
P_WAVE_E_0_1 = -0.146
P_WAVE_E_1_1 = 0.337
P_WAVE_E_2_1 = 0.198

P_WAVE_CONFORMAL_COEFFICIENTS_2 = (0.97, 0.11, -0.13, 0.47, 1.36)
P_WAVE_D_0_2 = 3.6
P_WAVE_D_1_2 = 2.4
P_WAVE_K_0_2 = 0.17
P_WAVE_E_0_2 = -0.16
P_WAVE_E_1_2 = 0.043
P_WAVE_E_2_2 = 0.06


def p_wave_generate_cot_phase(conformal_coeff, rho_mass=P_WAVE_RHO_MASS):
    conformal_polynomial = Polynomial(conformal_coeff)
    rho2 = rho_mass**2
    pion = PION_MASS**3
    s_0 = 1.43**2

    def cot_phase(s):
        conf = conformal_polynomial(omega(s, s_0=s_0, alpha=0.3))
        sqrt_s = csqrt(s)
        prefactor = sqrt_s * (rho2 - s) / 2 / momentum(s)**3
        remove_spurious = 2 * pion / rho2 / sqrt_s
        return prefactor * (remove_spurious + conf)

    return cot_phase


def p_wave_extend_phase(d_0, d_1, s_m=S_MATCHING):
    """Extend phase to be valid up to 2 GeV."""
    def wrapper(low_energy_phase):
        value = low_energy_phase(s_m)
        d_value = derivative(low_energy_phase, s_m, dx=1e-6)

        delta = ((d_value - d_0 * chebyshev_derivative(2, -1.0)
                  - d_1 * chebyshev_derivative(3, -1.0))
                 / chebyshev_derivative(1, -1.0))

        def high_energy_phase(s):
            chebyshev = cheb.chebval(omega_2(s, s_m=s_m), [0, delta, d_0, d_1])
            return value + chebyshev + delta - d_0 + d_1

        def phase(s):
            s = np.asarray(s)
            return np.piecewise(s,
                                [s <= s_m, s > s_m],
                                [low_energy_phase, high_energy_phase])

        return phase

    return wrapper


def generate_p_wave(conformal_coeff, d_0, d_1,
                    rho_mass=P_WAVE_RHO_MASS,
                    s_m=S_MATCHING):
    cot_phase = p_wave_generate_cot_phase(conformal_coeff, rho_mass=rho_mass)

    # @p_wave_extend_phase(d_0, d_1, s_m=s_m)
    def phase(s):
        return arccot2(cot_phase(s))

    amplitude = from_cot(PION_MASS)(cot_phase)
    amplitude_2 = second_sheet(PION_MASS)(amplitude)

    return phase, amplitude, amplitude_2


def p_wave_inelasticity(s, k_0, s_e=1.12**2):
    s = np.asarray(s)

    def above(x):
        return 1.0 - k_0 * (1.0 - x / s_e)**2

    return np.piecewise(s, [s < s_e, s >= s_e], [1, above])


p_wave_phase, p_wave, p_wave_2 = generate_p_wave(
                    conformal_coeff=P_WAVE_CONFORMAL_COEFFICIENTS_2,
                    d_0=P_WAVE_D_0_2, d_1=P_WAVE_D_1_2)
