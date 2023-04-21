import numpy as np
import pytest

from khuri import mandelstam


def test_sum():
    mass = 1.0
    virtuality = 30.0
    mandelstam_s = 10.0
    cos = np.linspace(-1.0, 1.0)
    mandelstam_t = mandelstam.t_vector_decay(mandelstam_s, cos, mass,
                                             virtuality)
    mandelstam_u = mandelstam.t_vector_decay(mandelstam_s, -cos, mass,
                                             virtuality)
    total = mandelstam_s + mandelstam_t + mandelstam_u
    expected = 3.0 * mass**2 + virtuality
    assert total == pytest.approx(expected)
