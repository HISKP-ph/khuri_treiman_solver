"""
========================================
Singularities in Khuri Treiman equations
========================================
"""
import numpy as np
from numpy.lib.scimath import sqrt  # handles negative reals correctly


CUBE_ROOT = (-1.0 + sqrt(3) * 1j) / 2.0


@np.vectorize
def cube_root(value):
    """Return the cube root."""
    return complex(value)**(1/3)


def cubic_equation(coeff, number: int):
    """Solve a cubic equation.

    Parameters
    ----------
    coeff:
        the coefficients of the equation. That is, the cubic equation is of
        the form:
        coeff[0] + coeff[1]*x + coeff[2]*x**2 + coeff[3]*x**3
    number: int
        pick one of three solutions of the equation via setting `number` to
        0, 1, or 2.
    """
    c0, c1, c2, c3 = coeff
    delta0 = c2**2 - 3.0 * c3 * c1
    delta1 = 2.0 * c2**3 - 9.0 * c3 * c2 * c1 + 27.0 * c3**2 * c0
    sqr = sqrt(delta1**2 - 4.0 * delta0**3)
    delta = cube_root((delta1 + sqr) / 2.0)
    threshold = 1e-20
    if np.any(abs(delta) < threshold):
        delta = cube_root((delta1 - sqr) / 2.0)
    delta *= CUBE_ROOT**number
    return -1.0 * (c2 + delta + delta0 / delta) / 3.0 / c3


def _coefficient_3(cos2):
    return cos2 - 1


def _coefficient_2(cos2, decay_mass_2, mass_2):
    return 2.0 * ((decay_mass_2 - 5.0 * mass_2)
                  - cos2 * (decay_mass_2 + 3.0 * mass_2))


def _coefficient_1(cos2, decay_mass_2, mass_2):
    return (cos2 * (decay_mass_2 + 3.0 * mass_2)**2
            - (decay_mass_2 - 5.0 * mass_2)**2)


def _coefficient_0(cos2, decay_mass_2, mass_2):
    return - 4.0 * cos2 * mass_2 * (decay_mass_2 - mass_2)**2


def coefficients(cos2, decay_mass_2, mass_2):
    """The coefficients of the cubic equation determining the singularities."""
    return (_coefficient_0(cos2, decay_mass_2, mass_2),
            _coefficient_1(cos2, decay_mass_2, mass_2),
            _coefficient_2(cos2, decay_mass_2, mass_2),
            _coefficient_3(cos2))


def singularity(cos2, decay_mass_2, mass_2, number: int):
    """Return the singularity.

    Parameters
    ----------
    cos2: float or array_like of floats
        the squared cosine of the scattering angle
    decay_mass_2: float
        the squared mass of the decaying particle
    mass_2: float
        the squared mass of one of the equal-mass particles in the final state
    number: int
        pick one of three singularities via setting `number` to 0, 1, or 2.
    """
    coeff = coefficients(cos2, decay_mass_2, mass_2)
    return cubic_equation(coeff, number)


def singularity_curve(decay_mass_2, mass_2, number, points=200, eps=1e-6):
    """Return the values of the singularity curve as parametrized by the angle.

    Parameters
    ----------
    decay_mass_2: float
        the squared mass of the decaying particle
    mass_2: float
        the squared mass of one of the equal-mass particles in the final state
    number: int
        pick one of three singularities via setting `number` to 0, 1, or 2.
    points: int, optional
        the number of sites
    eps: float, optional
        the endpoints of the curve correspond to the values `0 + eps`
        and `1 - eps` of the squared cosine of the scattering angle.
        (At 0 and 1 the cubic equation reduces to a quadratic one, so the
        general solution of the first stops working.)
    """
    cos2 = np.linspace(0 + eps, 1 - eps, points)
    return singularity(cos2, decay_mass_2, mass_2, number)


def singularity_curves(decay_mass_2, mass_2, points=200, eps=1e-6):
    """Return the values of the three singularity curves.

    Parameters
    ----------
    decay_mass_2: float
        the squared mass of the decaying particle
    mass_2: float
        the squared mass of one of the equal-mass particles in the final state
    points: int, optional
        the number of sites
    eps: float, optional
        the endpoints of the curve correspond to the values `0 + eps`
        and `1 - eps` of the squared cosine of the scattering angle.
        (At 0 and 1 the cubic equation reduces to a quadratic one, so the
        general solution of the first stops working.)

    Returns
    -------
    result: dict
        the three entries with keys 0, 1, and 2 contain each the values of one
        singularity curve
    """
    return {k: singularity_curve(decay_mass_2, mass_2, k, points, eps)
            for k in range(3)}
