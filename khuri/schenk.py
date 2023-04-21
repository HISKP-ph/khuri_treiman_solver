"""
===============================================================
The Schenk parametrization of pion+pion -> pion+pion scattering
===============================================================

The original form of the parametrization is given in
<A. Schenk, Absorption and dispersion of pions at finite temperature.
Nucl. Phys. B363, 97 (1991)>.

More recent values for the parameters are given in
<G. Colangelo, J. Gasser, H.Leutwyler, Nuclear Physics B 603 (2001) 125–179>.
"""
import functools

import numpy as np

from khuri.amplitude import from_cot
from khuri.phase_space import rho


def tan_phase(mandelstam_s, isospin, pion_mass, peak, coefficients):
    """Return the tangent of the phase.

    Parameters
    ----------
    mandelstam_s: float or array_like
        the value of the Mandelstam variable s, at which the phase is evaluated
    isospin: int
        the value of the isospin of the partial wave (one of 0, 1, 2)
    pion_mass: float
        the value of the pion mass
    peak: float
        the value of Mandelstam s at which the phase equals 90 degree
    coefficients: iterable
        the (dimensionless) coefficients of the polynomial in the
        parametrization
    """
    mandelstam_s = np.asarray(mandelstam_s)
    return np.piecewise(mandelstam_s,
                        [mandelstam_s != peak],
                        [lambda x: _tan_phase(x, isospin, pion_mass, peak,
                                              coefficients),
                         lambda _: np.inf])


def to_spin(isospin):
    """Return the spin corresponding to a given isospin."""
    spin = {
        0: 0,
        1: 1,
        2: 0,
    }
    return spin[isospin]


def _tan_phase(mandelstam_s, isospin, pion_mass, peak, coefficients):
    threshold = 4.0 * pion_mass**2
    param = mandelstam_s / threshold - 1.0
    phase_space = rho(pion_mass, mandelstam_s)
    polynomial = sum(coeff * param**i for i, coeff in enumerate(coefficients))
    peak_factor = (threshold - peak) / (mandelstam_s - peak)
    return phase_space * param**to_spin(isospin) * polynomial * peak_factor


def partial_wave(mandelstam_s, isospin, pion_mass, peak, coefficients):
    """Return the partial wave amplitude.

    Parameters
    ----------
    mandelstam_s: float or array_like
        the value of the Mandelstam variable s, at which the phase is evaluated
    isospin: int
        the value of the isospin of the partial wave (one of 0, 1, 2)
    pion_mass: float
        the value of the pion mass
    peak: float
        the value of Mandelstam s at which the phase equals 90 degree
    coefficients: iterable
        the (dimensionless) coefficients of the polynomial in the
        parametrization
    """
    @from_cot(pion_mass, spin=to_spin(isospin))
    def amp(s_value):
        tan = tan_phase(s_value, isospin, pion_mass, peak, coefficients)
        return inverse(tan)
    return amp(mandelstam_s)


def inverse(values):
    return np.piecewise(values, [values != 0.0],
                        [lambda x: 1.0 / x, lambda _: np.inf])


def literature(func):
    @functools.wraps(func)
    def wrapper(mandelstam_s, isospin, pion_mass):
        params = literature_values(isospin, pion_mass)
        return func(mandelstam_s, isospin, pion_mass, *params)
    return wrapper


@literature
def tan_phase_lit(*args, **kwargs):
    """The tangent of the phase w/ parameters set to literature values.

    Parameters
    ----------
    mandelstam_s: float or array_like
        the value of the Mandelstam variable s, at which the phase is evaluated
    isospin: int
        the value of the isospin of the partial wave (one of 0, 1, 2)
    pion_mass: float
        the value of the pion mass
    """
    return tan_phase(*args, **kwargs)


@literature
def partial_wave_lit(*args, **kwargs):
    """The partial wave amplitude w/ parameters set to literature values.

    Parameters
    ----------
    mandelstam_s: float or array_like
        the value of the Mandelstam variable s, at which the phase is evaluated
    isospin: int
        the value of the isospin of the partial wave (one of 0, 1, 2)
    pion_mass: float
        the value of the pion mass
    """
    return partial_wave(*args, **kwargs)


def literature_values(isospin, pion_mass):
    """Return literature values of parameters.

    Parameters
    ----------
    isospin: int
        the value of the isospin of the partial wave (one of 0, 1, 2)
    pion_mass: float
        the value of the pion mass
    """
    if isospin == 0:
        return 36.77 * pion_mass**2, (0.22, 0.268, -0.139e-1, -0.139e-2)
    if isospin == 1:
        return 30.72 * pion_mass**2, (0.379e-1, 0.14e-4, -0.673e-4, 0.163e-7)
    if isospin == 2:
        return -21.62 * pion_mass**2, (-0.444e-1, -0.857e-1, -0.221e-2,
                                       -0.129e-3)
    raise ValueError(f'There is no partial wave with isospin {isospin}.')
