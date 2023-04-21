import numpy as np
import pytest

from khuri.khuri_treiman import Real

LOWER_VALUE = 4.0
UPPER_VALUE = 100.0


@pytest.fixture
def real():
    """Linear curve along the real axis"""
    return Real(LOWER_VALUE, UPPER_VALUE)


class TestBoundaries:
    def test_lower(self, real):
        assert real.lower() == 0.0

    def test_upper(self, real):
        assert real.upper() == 1.0

    def test_boundaries(self, real):
        assert real.boundaries() == [0.0, 1.0]


class TestHits:
    def test_miss(self, real):
        assert real.hits(1.0j) is None

    def test_hits(self, real):
        assert real.hits(10.0) == (0.0, 1.0)


class TestEvaluate:
    def test_curve(self, real):
        assert real.curve_func(0.0) == LOWER_VALUE
        middle = (UPPER_VALUE + LOWER_VALUE) / 2.0
        assert real.curve_func(0.5) == pytest.approx(middle)
        assert real.curve_func(1.0) == UPPER_VALUE

    def test_derivative(self, real):
        derivative = UPPER_VALUE - LOWER_VALUE
        for x in np.linspace(0.0, 1.0, 20):
            assert real.derivative_func(x) == pytest.approx(derivative)
