import numpy as np

import khuri.phase_space as kps


MASS = 1.0
THRESHOLD = 4.0 * MASS**2


def rho(s):
    return kps.rho(MASS, s)


def sigma(s):
    return kps.sigma(MASS, s)


def test_compare_along_real_axis():
    """Values above threshold should be independent of the cut structure."""
    mandelstam_s = np.linspace(THRESHOLD, 1000)
    comparision = np.allclose(rho(mandelstam_s), sigma(mandelstam_s))
    assert comparision, test_compare_along_real_axis.__doc__
