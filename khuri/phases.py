"""
==================================
High energy continuation of phases
==================================
"""
from functools import wraps

import numpy as np
import numdifftools as nd
# from scipy.misc import derivative


def asymptotic1(matching_point, limit=np.pi):
    """Return phase that approaches a constant at high energies.

    The decorated phase is used below `matching_point`, above which `limit`
    is approached eventually. The transition is one times differentiable.
    """
    def wrapper(phase):
        value = phase(matching_point)
        deriv = nd.Derivative(phase, step=1e-6)
        d_value = deriv(matching_point)
        parameter1 = (limit - value)**2 / matching_point / d_value
        parameter2 = (limit - value) / matching_point / d_value - 1

        if parameter2 <= -1:
            raise ValueError("the continuation of the phase to high energies"
                             " exhibits a pole above the matching point")

        @wraps(phase)
        def continued_phase(s):
            if s < matching_point:
                return phase(s)
            return limit - parameter1 / (parameter2 + s / matching_point)

        return np.vectorize(continued_phase, doc=phase.__doc__)

    return wrapper


def asymptotic2(lower, upper, limit=np.pi):
    """Return phase that equals `phase` below `lower` and `limit` above `upper`.

    Between `lower` and `upper` there is a in [`lower`, `upper`] one
    times continuously differentiable smooth connection.

    Note
    ----
    The matching is achieved using a third-order polynomial, which
    requires to use also the phase as a part of the function in
    [`lower`, `upper`].
    """
    def wrapper(phase):
        @wraps(phase)
        def continued_phase(s):
            if s < lower:
                return phase(s)
            if s < upper:
                c = _connect(lower, upper, s)
                return (1.0 - c) * phase(s) + c * limit
            return limit

        return np.vectorize(continued_phase, doc=phase.__doc__)

    return wrapper


def asymptotic3(lower, upper, limit=np.pi):
    """Return phase that equals `phase` below `lower` and `limit` above `upper`.

    Between `lower` and `upper` there is a in [`lower`, `upper`] one
    times continuously differentiable smooth connection.

    Note
    ----
    The matching is achieved using a polynomial of degree three, which might
    lead to a huge bump in [`lower`, `upper`], depending on the input
    parameters.
    """
    def wrapper(phase):
        value = phase(lower)
        deriv = nd.Derivative(phase, step=1e-6)
        d_value = deriv(matching_point)
        vector = np.array([value, limit, d_value, 0.0])
        matrix = np.array([
            [1.0, lower, lower**2, lower**3],
            [1.0, upper, upper**2, upper**3],
            [0.0, 1.0, 2.0 * lower, 3.0 * lower**2],
            [0.0, 1.0, 2.0 * upper, 3.0 * upper**2],
        ])
        coefficients = np.linalg.inv(matrix) @ vector

        @wraps(phase)
        def continued_phase(s):
            if s < lower:
                return phase(s)
            if s < upper:
                return coefficients @ np.array([1.0, s, s**2, s**3])
            return limit

        return np.vectorize(continued_phase, doc=phase.__doc__)

    return wrapper


def _connect(lower, upper, s):
    a = (s - lower)**2
    b = 3.0*upper - 2.0*s - lower
    c = upper - lower
    return a * b / c**3
