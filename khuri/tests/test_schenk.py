import numpy as np
import pytest

from khuri.schenk import partial_wave_lit, tan_phase_lit, literature_values
from khuri.phase_space import rho


ISOSPINS = (0, 1, 2)
ISOSPIN_IDS = tuple(f'isospin {i}' for i in ISOSPINS)


@pytest.mark.parametrize('isospin', ISOSPINS, ids=ISOSPIN_IDS)
def test(isospin):
    mandelstam_s = np.linspace(0.3, 0.5)**2
    pion_mass = 0.135
    tan = tan_phase_lit(mandelstam_s, isospin, pion_mass)
    partial_wave = partial_wave_lit(mandelstam_s, isospin, pion_mass)
    tan2 = np.tan(np.angle(partial_wave))
    assert np.allclose(tan, tan2)


@pytest.mark.parametrize('func', (partial_wave_lit, tan_phase_lit),
                         ids=('partial wave', 'tangent'))
def test_isospin_4(func):
    with pytest.raises(ValueError):
        func(0.3**2, 4, 0.135)


@pytest.mark.parametrize('isospin', ISOSPINS, ids=ISOSPIN_IDS)
def test_peak(isospin):
    peak = literature_values(isospin, pion_mass=1.0)[0]
    tan = tan_phase_lit(peak, isospin, pion_mass=1.0)
    partial_wave = partial_wave_lit(peak, isospin, pion_mass=1.0)
    assert tan == np.inf
    assert partial_wave == pytest.approx(1j / rho(1.0, peak))


@pytest.mark.parametrize('isospin', ISOSPINS, ids=ISOSPIN_IDS)
def test_threshold(isospin):
    threshold = 4.0
    tan = tan_phase_lit(threshold, isospin, pion_mass=1.0)
    assert tan == 0.0
    if isospin == 1:
        partial_wave = partial_wave_lit(threshold, isospin, pion_mass=1.0)
        assert partial_wave == 0.0
