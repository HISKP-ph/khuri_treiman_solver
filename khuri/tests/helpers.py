"""
===========================================
Facilites that are used in different tests.
===========================================
"""
from typing import Any

import numpy as np


def schwarz(func, s, atol=1e-14, rtol=0):
    """Check whether `func` fulfills Schwarz reflection principle."""
    assert np.all(s.imag != 0.0), "tried to check Schwarz at real axis"
    a = func(s)
    b = func(s.conjugate()).conjugate()
    assert np.allclose(a, b, atol=atol, rtol=rtol), compare("a", a, "b", b)


def connected(sheet1, sheet2, s, epsilon=1e-15, rtol=1e-5, atol=1e-8):
    """Check if `sheet1` and `sheet2` are continuously connected."""
    s_plus = s + epsilon * 1j
    s_minus = s_plus.conjugate()

    first_sheet_above = sheet1(s_plus)
    first_sheet_below = sheet1(s_minus)

    second_sheet_above = sheet2(s_plus)
    second_sheet_below = sheet2(s_minus)

    message1 = compare("first_sheet_above", first_sheet_above,
                       "second_sheet_below", second_sheet_below)
    assert np.allclose(first_sheet_above, second_sheet_below, rtol=rtol,
                       atol=atol), message1
    message2 = compare("first_sheet_below", first_sheet_below,
                       "second_sheet_above", second_sheet_above)
    assert np.allclose(first_sheet_below, second_sheet_above, rtol=rtol,
                       atol=atol), message2


def compare(name1: str, val1: Any, name2: str, val2: Any) -> str:
    """Return error message."""
    return f'{name1} = {val1}\n<=>\n{name2} = {val2}'
