import numpy as np
import pytest

from khuri.curved_omnes import CurvedOmnes
from khuri.iam import nlo
from khuri.khuri_treiman import Adaptive, Real, VectorDecay
from khuri.omnes import generate_omnes, second_sheet


def amplitude(mandelstam_s):
    return nlo(1.0, mandelstam_s, 0.09 / 0.14, 6.0)


def phase(mandelstam_s):
    return np.angle(amplitude(mandelstam_s))


def make_curved(curve):
    cut = 1e2
    omnes_function = generate_omnes(phase, 4.0, phase(cut), cut)
    return CurvedOmnes(omnes_function, amplitude, curve)


def test_real():
    curved = make_curved(Real(4.0, 1e5))
    omn = curved.original()
    mandelstam_s = [-10.0, 0.0, 10.0, 1e6, 10.0 - 10.0j, 10.0 + 10.0j]
    assert np.all(curved(mandelstam_s) == omn(mandelstam_s))


@pytest.mark.parametrize('curve',
                         (Adaptive, VectorDecay),
                         ids=('Adaptive', 'VectorDecay'))
def test_rectangle(curve):
    curved = make_curved(curve(1.0, 30.0, 1e5))
    omn = curved.original()

    mandelstam_s_1st_sheet = [-10.0, 0.0, 10.0, 1e6, 10.0 + 10.0j, 10.0 - 1e2j]
    assert np.all(curved(mandelstam_s_1st_sheet)
                  == omn(mandelstam_s_1st_sheet))

    mandelstam_s_2nd_sheet = [5.0 - 1.0j, 8.0 - 1e-2j]
    assert np.all(curved(mandelstam_s_2nd_sheet)
                  == second_sheet(omn, amplitude, mandelstam_s_2nd_sheet))
