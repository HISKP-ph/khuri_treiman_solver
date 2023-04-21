import itertools

import numpy as np
import pytest

from khuri.omnes import generate_omnes, second_sheet
from khuri.gsl import IntegrationRoutine
from khuri import madrid, iam
from khuri.tests.helpers import schwarz, connected
from khuri.tests.test_phases import PHASES, MATCHING_POINTS


# A small value is added to avoid the evaluation of the (non-global)
# parametrization by the Madrid group at threshold
THRESHOLD = (2.0 * madrid.PION_MASS)**2 + 1e-6

PION_MASS = 0.14

SAMPLE_POINTS = 3


def test_unknown_integration():
    with pytest.raises(ValueError):
        generate_omnes(PHASES[0], threshold=THRESHOLD,
                       integration_routine='routine')


def all_omnes_for_phase(phase):
    return (
        generate_omnes(phase, threshold=THRESHOLD, constant=np.pi, cut=1e10),
        generate_omnes(phase, threshold=THRESHOLD),
        generate_omnes(phase, threshold=THRESHOLD, constant=np.pi, cut=5e3,
                       integration_routine=IntegrationRoutine.qag),
        generate_omnes(phase, threshold=THRESHOLD,
                       integration_routine=IntegrationRoutine.qag),
    )


OMNES_LIST = list(itertools.chain.from_iterable(all_omnes_for_phase(phase)
                                                for phase in PHASES))
PHASE_LIST = list(itertools.chain.from_iterable(itertools.repeat(phase, 4)
                                                for phase in PHASES))
OMNES_IDS = [
    'p_wave_cquad_cut',
    'p_wave_cquad',
    'p_wave_qag_cut',
    'p_wave_qag',
    'p_wave_global_cquad_cut',
    'p_wave_global_cquad',
    'p_wave_global_qag_cut',
    'p_wave_global_qag',
    's_wave_global_cquad_cut',
    's_wave_global_cquad',
    's_wave_global_qag_cut',
    's_wave_global_qag',
]


def amplitude(mandelstam_s):
    return iam.nlo(PION_MASS, mandelstam_s, 0.09, 6.0)


SHEET_1 = generate_omnes(lambda x: np.angle(amplitude(x)),
                         threshold=4.0 * PION_MASS**2,
                         constant=np.pi, cut=1e10)


def sheet_2(mandelstam_s):
    return second_sheet(SHEET_1, amplitude, mandelstam_s)


def test_2nd_sheet_connected():
    mandelstam_s = np.linspace(0.7, 1.2, SAMPLE_POINTS)**2
    connected(SHEET_1, sheet_2, mandelstam_s)


def test_2nd_sheet_schwarz():
    schwarz(sheet_2, 1.3 + 4.0j)


@pytest.mark.parametrize('function', OMNES_LIST, ids=OMNES_IDS)
def test_omnes_at_zero(function):
    assert function(0) == pytest.approx(1.0)

@pytest.mark.parametrize('function', OMNES_LIST, ids=OMNES_IDS)
def test_derivative_at_zero(function):
    assert isinstance(function.derivative_at_zero, float)


@pytest.mark.parametrize('function, phase', list(zip(OMNES_LIST, PHASE_LIST)),
                         ids=OMNES_IDS)
def test_phase(function, phase):
    """Check if the phase of the Omnes function agrees with the input phase."""
    # The phase of the s-wave exceeds pi close to kaon threshold, thus
    # for simplicity only the region below is checked.
    mandelstam_s = np.linspace(THRESHOLD, min(MATCHING_POINTS) - 0.02,
                               SAMPLE_POINTS)
    omnes_phase = np.angle(function(mandelstam_s))
    phase_values = phase(mandelstam_s)
    assert np.allclose(omnes_phase[1:], phase_values[1:])


@pytest.mark.parametrize('function', OMNES_LIST, ids=OMNES_IDS)
def test_schwarz(function):
    """Check if the Omnes function fulfills the Schwarz reflection princple."""
    real_part = 0.760**2
    imaginary_parts = [x for x in np.linspace(-1e4, 1e4, SAMPLE_POINTS)
                       if x != 0.0]
    mandelstam_s = real_part + 1j * np.asarray(imaginary_parts)
    schwarz(function, mandelstam_s)
