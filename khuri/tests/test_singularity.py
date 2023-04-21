import pytest
import numpy as np

from khuri import singularity


def test_cube_root():
    def cube_root(number):
        return np.exp(np.emath.log(number) / 3)
    for i in (-1, 0.5, 2j):
        assert singularity.cube_root(i) == pytest.approx(cube_root(i))


def test_cubic_equation():
    coeff1 = (0, -0.75, 1, 1)
    sol1 = {k: singularity.cubic_equation(coeff1, k) for k in range(3)}
    assert all(v.imag < 1e-10 for v in sol1.values())
    assert ({k: v.real for k, v in sol1.items()}
            == pytest.approx({0: -1.5, 1: 0.5, 2: 0.0}))


def test_singularity_curves_spacelike():
    curves = singularity.singularity_curves(decay_mass_2=-10, mass_2=1)
    for curve in curves.values():
        assert np.all(curve.imag < 1e-8)
        assert np.all(curve.real <= 0.0)


def test_singularity_curves_timelike():
    curves = singularity.singularity_curves(decay_mass_2=30, mass_2=1)
    # one solution extends along the negative real axis
    assert np.all(curves[1].imag < 1e-8)
    assert np.all(curves[1].real <= 0.0)
    # the other solutions are a complex conjugated pair
    inner = slice(1, -1)
    assert curves[0].real[inner] == pytest.approx(curves[2].real[inner])
    assert -curves[0].imag[inner] == pytest.approx(curves[2].imag[inner])
