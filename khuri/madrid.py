"""
==========================================================
The I=J=1 pi pi->pi pi Madrid partial wave parametrization
==========================================================

The notation
from the paper is used in the implementation. The values for
the non-dimensionless fit parameters/default arguments are given in GeV, so
Mandelstam s should accordingly be given in units of (GeV)^2.
"""
import numpy as np


# constraint fit parameters
B0_CFD = 1.043
B1_CFD = 0.19
LAMBDA1_CFD = 1.39
LAMBDA2_CFD = -1.70

# unconstraint fit parameters
B0_UFD = 1.055
B1_UFD = 0.15
LAMBDA1_UFD = 1.57
LAMBDA2_UFD = -1.96

S0 = 1.05**2  # in (GeV)^2

# The masses are parameters that determine the shape of the phase.
# They are not necessarily the physical masses of the particles.
KAON_MASS = 0.496  # in GeV
RHO_MASS = 0.7736  # in GeV
PION_MASS = 0.13957  # in GeV

# The parametrization is not valid for energies below LOWEST
# or above HIGHEST.
LOWEST = 2.0*PION_MASS  # in GeV
HIGHEST = 1.42  # in GeV


def _momentum(s, m):
    """The center of mass momentum.

    Parameters
    ----------
    s
        the squared cm energy
    m
        the pion mass
    """
    return np.sqrt(0.25*s - m**2)


def _w(s, s0=S0):
    """Simplify the calculation of the phase.

    Parameters
    ----------
    s
        the squared cm energy
    s0
    """
    a = np.sqrt(s)
    b = np.sqrt(s0 - s)
    return (a-b) / (a+b)


def phase_low(s, m_pi, m_rho, b0, b1, s0=S0):
    """Calculate the phase for s <= 4m_K^2

    Parameters
    ----------
    s
        the squared cm energy
    m_pi
        the pion mass
    m_rho
        the rho mass
    b0
        a parameter determined by fits
    b1
        a parameter determined by fits
    s0
    """
    en = np.sqrt(s)
    coeff = en * (m_rho**2-s) / (2*_momentum(s, m_pi)**3)
    summ = 2*m_pi**3 / (m_rho**2*en) + b0 + b1*_w(s, s0)
    return np.arctan2(1, coeff*summ)


def phase_high(s, m_k, lambda_0, lambda_1, lambda_2):
    """Calculate the phase for (1420 GeV)^2 > s > 4m_K^2

    Parameters
    ----------
    s
        the squared cm energy
    m_k
        the kaon mass
    lambda_0
        a parameter determined by fits
    lambda_1
        a parameter determined by fits
    lambda_2
        a parameter determined by fits
    """
    temp = np.sqrt(s)/(2*m_k) - 1
    return lambda_0 + lambda_1*temp + lambda_2*temp**2


def lambda0(m_pi, m_k, m_rho, b0, b1, s0):
    s = 4.0*m_k**2
    return phase_low(s, m_pi, m_rho, b0, b1, s0)


@np.vectorize
def phase(
        s,
        m_pi=PION_MASS,
        m_k=KAON_MASS,
        m_rho=RHO_MASS,
        b0=B0_CFD,
        b1=B1_CFD,
        lambda_1=LAMBDA1_CFD,
        lambda_2=LAMBDA2_CFD,
        s0=S0):
    """Calculate the Madrid phase for the I=J=1 pi pi-> pi pi partial wave.

    Parameters
    ----------
    s
        the squared cm energy
    m_pi
        the pion mass
    m_k
        the kaon mass
    m_rho
        the rho mass
    b0
        a parameter determined by fits
    b1
        a parameter determined by fits
    lambda_1
        a parameter determined by fits
    lambda_2
        a parameter determined by fits
    s0
    """
    if s < LOWEST**2:
        raise ValueError(f'{s} is below threshold')
    elif s <= 4*m_k**2:
        return phase_low(s, m_pi, m_rho, b0, b1, s0)
    elif s <= HIGHEST**2:
        l0 = lambda0(m_pi, m_k, m_rho, b0, b1, s0)
        return phase_high(s, m_k, l0, lambda_1, lambda_2)
    else:
        raise ValueError(f'{s} is above region, in which the parametrization'
                         ' of the phase is valid')
