"""
============================================================
General two particle into two particle scattering amplitudes
============================================================

It is assumed that all particles have the same mass.
"""
from typing import Optional
from functools import wraps

import numpy as np

from khuri.phase_space import rho


def from_cot(mass: float, spin: Optional[int] = None):
    """Obtain amplitude from cotangent of the phase.

    Parameters
    ----------
    mass: float
        the mass of an individual particle
    spin: int, optional
        partial waves with non-zero angular momentum vanish at threshold.
        If `spin` is specified and non-zero, this behaviour is enforced.
    """
    def wrapper(cot_phase):
        def amplitude(*args):
            return 1 / (cot_phase(*args) - 1j) / rho(mass, args[-1])
        if not spin:
            return wraps(cot_phase)(amplitude)

        @wraps(cot_phase)
        def modified_amplitude(*args):
            threshold = 4.0 * mass**2
            mandelstam_s = np.asarray(args[-1], dtype=np.complex128)
            return np.piecewise(mandelstam_s,
                                [mandelstam_s != threshold],
                                [lambda x: amplitude(*args[:-1], x),
                                 lambda _: 0.0])
        return modified_amplitude
    return wrapper


def from_phase(mass, inelasticity=None):
    """Obtain amplitude from the phase."""
    if inelasticity is None:
        def inelasticity(_):
            return 1.0

    def wrapper(phase):
        @wraps(phase)
        def amplitude(*args):
            phase_value = phase(*args)
            inelastic = inelasticity(args[-1])
            numerator = inelastic * np.exp(2j * phase_value) - 1
            return numerator / rho(mass, args[-1]) / 2j
        return amplitude
    return wrapper


def second_sheet(mass):
    """Obtain the amplitude on the second sheet from the one on the first."""
    def wrapper(amplitude):
        @wraps(amplitude)
        def second(*args):
            first = amplitude(*args)
            return first / (1 + 2j * rho(mass, args[-1]) * first)

        return second
    return wrapper


def first_sheet(mass):
    """Obtain the amplitude on the first sheet from the one on the second."""
    def wrapper(amplitude):
        @wraps(amplitude)
        def first(*args):
            second = amplitude(*args)
            return second / (1 - 2j * rho(mass, args[-1]) * second)

        return first
    return wrapper
