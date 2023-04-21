"""
===========================================================
Two-particle phase spaces with different analytic structure
===========================================================
"""
from numpy.lib.scimath import sqrt  # handles negative reals correctly
from numpy import imag, vectorize


@vectorize
def signum_im(x):
    """Signum of the imaginary part of a number."""
    return 1 if imag(x) >= 0.0 else -1


def alt_sqrt(x):
    """Square root function with cut on positive real axis.

    Notes
    -----
    numpy.sqrt/numpy.lib.scimath.sqrt develop a cut on the negative real axis.
    """
    return signum_im(x)*sqrt(x)


def rho(mass, s):
    """The two-body phase space (cuts along [4m**2,infty) and (-infty,0])."""
    return alt_sqrt(1 - 4*mass**2/s)


def sigma(mass, s):
    """The two-body phase space (cut along [0,4m**2])."""
    return sqrt(1 - 4*mass**2/s)
