import functools

from khuri.gsl import Settings, IntegrationRoutine
from _khuri_omnes import *
from _khuri_omnes import __doc__ as module_docstring


__doc__ = module_docstring


def _factory(func):
    callable1, callable2 = func()

    @functools.wraps(func)
    def wrapper(*args,
                integration_routine=IntegrationRoutine.cquad,
                **kwargs):
        if integration_routine == IntegrationRoutine.cquad:
            return callable1(*args, **kwargs)
        if integration_routine == IntegrationRoutine.qag:
            return callable2(*args, **kwargs)
        raise ValueError('unknown integration routine')

    return wrapper


@_factory
def generate_omnes():
    """Generate an Omnes function.

    Parameters
    ----------
    phase: callable
        the phase of the Omnes function
    threshold: float
        the threshold in the s-plane
    constant: float, optional
        if specified, the phase is set to constant above `cut`
    cut: float, optional
        if specified, the phase is set to constant above `cut`
    minimal_distance: float, optional
        half the width of a band around the cut. For arguments of the Omnes
        function in this band, a different prescription is used for the
        evaluation of the Omnes function to take care of the singularity in the
        integral.
    conf: Settings, optional
        the settings for the integration routine
    integration_routine: IntegrationRoutine, optional
        specify an adaptive integration routine to be used to compute the
        integral

    Returns
    -------
    A callable that yields the values of the requested Omnes function.
    """
    return OmnesCquad, OmnesQag


def second_sheet(omnes_function, amplitude, mandelstam_s):
    if isinstance(omnes_function, OmnesCquad):
        return second_sheet_cquad(omnes_function, amplitude, mandelstam_s)
    if isinstance(omnes_function, OmnesQag):
        return second_sheet_qag(omnes_function, amplitude, mandelstam_s)
    raise ValueError('unknown type of omnes function')
