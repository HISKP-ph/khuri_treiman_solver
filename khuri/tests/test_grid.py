import pytest

from khuri.khuri_treiman import Real, GridReal

X_SIZE = 20
Z_SIZE = 5


@pytest.fixture
def real_grid():
    curve = Real(4.0, 50.0)
    x_sizes = (X_SIZE, )
    return GridReal(curve, x_sizes, Z_SIZE)


def test_evaluation(real_grid):
    x_index = 8
    z_index = 1
    point = real_grid(x_index, z_index)
    assert real_grid.x(x_index) == point.x
    assert real_grid.derivative(x_index) == point.x_derivative
    assert real_grid.z(z_index) == point.z


def test_x_size(real_grid):
    assert real_grid.x_size() == X_SIZE


def test_z_size(real_grid):
    assert real_grid.z_size() == Z_SIZE


def test_x_parameter(real_grid):
    assert real_grid.x_parameter_lower() == pytest.approx(0.0)
    assert real_grid.x_parameter_upper() == pytest.approx(1.0)
