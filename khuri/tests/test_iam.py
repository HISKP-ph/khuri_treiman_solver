import numpy as np

from khuri import iam
from khuri.phase_space import rho
from khuri.tests.helpers import schwarz


PION_MASS = 0.139


def nlo(mandelstam_s):
    lec = 6.0
    pion_decay_chiral = 0.08
    return iam.nlo(PION_MASS, mandelstam_s, pion_decay_chiral, lec)


def test_schwarz():
    real_part = 0.760**2
    imaginary_parts = np.linspace(-1e4, 1e4, 20)
    mandelstam_s = real_part + 1j * imaginary_parts
    schwarz(nlo, mandelstam_s)


def test_unitarity():
    mandelstam_s = np.linspace(5.0, 40.0, 50) * PION_MASS**2
    one = np.imag(nlo(mandelstam_s))
    two = rho(PION_MASS, mandelstam_s) * np.abs(nlo(mandelstam_s))**2
    assert np.allclose(one, two)
