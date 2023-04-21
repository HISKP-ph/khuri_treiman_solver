import numpy as np
import pytest

import khuri.madrid
import khuri.madrid_global
import khuri.phases


P_WAVE_MATCHING_POINT = 1.2**2
KAON_THRESHOLD = 4.0 * khuri.madrid_global.KAON_MASS**2


@khuri.phases.asymptotic1(matching_point=P_WAVE_MATCHING_POINT)
def p_wave_phase(s):
    return khuri.madrid.phase(s)


@khuri.phases.asymptotic1(matching_point=P_WAVE_MATCHING_POINT)
def p_wave_phase_global(s):
    return khuri.madrid_global.p_wave_phase(s)


@khuri.phases.asymptotic1(matching_point=KAON_THRESHOLD, limit=2.0 * np.pi)
def s_wave_phase_global(s):
    return khuri.madrid_global.s_wave_phase(s)


PHASES = [
    p_wave_phase,
    p_wave_phase_global,
    s_wave_phase_global
]
ORIGINAL = [
    khuri.madrid.phase,
    khuri.madrid_global.p_wave_phase,
    khuri.madrid_global.s_wave_phase
]
LIMITS = [np.pi, np.pi, 2.0 * np.pi]
MATCHING_POINTS = [
    P_WAVE_MATCHING_POINT,
    P_WAVE_MATCHING_POINT,
    KAON_THRESHOLD
]


class TestAsymptotic1:
    @pytest.mark.parametrize('phase, limit', list(zip(PHASES, LIMITS)))
    def test_limit(self, phase, limit):
        """Test if the phase approaches the desired limit."""
        assert phase(1e12) == pytest.approx(limit)

    #@pytest.mark.parametrize('phase, original, matching_point',
    #                         list(zip(PHASES, ORIGINAL, MATCHING_POINTS)))
    #def test_agreement(self, phase, original, matching_point):
    #    """Test agreement with the original phase."""
    #    values = np.linspace(0.3**2, matching_point - 1e-5, 500)
    #    original = original(values)
    #    wrapper = phase(values)
    #    assert np.all(original == wrapper)
