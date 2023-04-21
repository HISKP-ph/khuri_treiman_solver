import numpy as np
import pytest

from khuri import madrid_global as mg
from khuri.tests.helpers import schwarz, connected


THRESHOLD = 4.0 * mg.PION_MASS**2


@pytest.fixture
def vertical_line():
    """Vertical line in the complex s-plane."""
    real_part = 0.76**2
    imaginary_parts = np.linspace(-1e4, 1e4, 100)
    return real_part + 1j * imaginary_parts


@pytest.mark.parametrize('function',
                         [mg.p_wave, mg.p_wave_2, mg.s_wave, mg.s_wave_2],
                         ids=['p_wave', 'p_wave_2', 's_wave', 's_wave_2'])
def test_schwarz(function, vertical_line):
    schwarz(function, vertical_line)


@pytest.mark.parametrize('first, second, upper',
                         [(mg.p_wave, mg.p_wave_2, 1.4),
                          (mg.s_wave, mg.s_wave_2, 2.0 * mg.KAON_MASS)],
                         ids=['p_wave', 's_wave'])
def test_connected(first, second, upper):
    real_axis = np.linspace(THRESHOLD, upper**2, 100)
    connected(first, second, real_axis)
