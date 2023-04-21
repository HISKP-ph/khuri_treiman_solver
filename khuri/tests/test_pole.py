import numpy as np
import pytest

from khuri.tests.test_iam import PION_MASS, nlo
import khuri.pole as kp


# guesses for rho resonance characteristics
RHO = 0.77 - 0.14j
COUPLING = 6.0


@pytest.fixture
def pole_position():
    return np.sqrt(kp.pole(nlo, PION_MASS, RHO**2))


@pytest.fixture
def coupling_value():
    return np.sqrt(kp.coupling(nlo, PION_MASS, RHO**2))


def test_pole(pole_position):
    assert abs(pole_position - RHO) < 0.2
    assert pole_position.imag < 0


def test_coupling(coupling_value, pole_position):
    assert abs(coupling_value - COUPLING) < 1.0
    coupling2 = np.sqrt(kp.coupling(nlo, PION_MASS, pole_position**2,
                                    compute_pole=False))
    assert coupling_value == pytest.approx(coupling2)


def test_pole_and_coupling(pole_position, coupling_value):
    pole, coupling_squared = kp.pole_and_coupling(nlo, PION_MASS, RHO**2)
    assert pole == pytest.approx(pole_position**2)
    assert coupling_squared == pytest.approx(coupling_value**2)


def test_failure():
    with pytest.raises(kp.PoleError):
        kp.pole(lambda _: 0.0, PION_MASS, RHO**2)
