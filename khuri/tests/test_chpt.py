import numpy as np
import pytest

from khuri import chpt
from khuri.phase_space import rho
from khuri.tests.helpers import schwarz


PION_MASS = 0.139
PION_DECAY_CHIRAL = 0.08


def lo(mandelstam_s):
    return chpt.lo(PION_MASS, mandelstam_s, PION_DECAY_CHIRAL)


def nlo(mandelstam_s):
    lec = 6.0
    return chpt.nlo(PION_MASS, mandelstam_s, PION_DECAY_CHIRAL, lec)


AMPLITUDES = (lo, nlo)


@pytest.mark.parametrize('function', AMPLITUDES)
def test_schwarz(function):
    real_part = 0.760**2
    imaginary_parts = np.linspace(-1e4, 1e4, 20)
    mandelstam_s = real_part + 1j * imaginary_parts
    schwarz(function, mandelstam_s)


def test_perturbative_unitarity():
    mandelstam_s = np.linspace(5.0, 40.0, 50) * PION_MASS**2
    one = np.imag(nlo(mandelstam_s))
    two = rho(PION_MASS, mandelstam_s) * np.abs(lo(mandelstam_s))**2
    assert np.allclose(one, two)
