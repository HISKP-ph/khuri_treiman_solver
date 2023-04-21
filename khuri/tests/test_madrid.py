import numpy as np
import pytest

from khuri import madrid


def test_evaluate_below():
    with pytest.raises(ValueError):
        value = madrid.LOWEST / 2.0
        madrid.phase(value**2)


def test_evaluate_above():
    with pytest.raises(ValueError):
        value = madrid.HIGHEST * 2.0
        madrid.phase(value**2)


def test_ascending():
    """Test if the phase of the p wave is ascending around the resonance."""
    energies = np.linspace(madrid.LOWEST + 0.01, 0.8)
    phases = madrid.phase(energies**2)
    assert all(a <= b for a, b in zip(phases, phases[1:]))
