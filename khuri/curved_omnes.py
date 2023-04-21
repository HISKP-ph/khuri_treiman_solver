from _khuri_curved_omnes import *
from _khuri_curved_omnes import __doc__ as module_docstring

from khuri.khuri_treiman import Adaptive, VectorDecay, Real


__doc__ = module_docstring


class CurvedOmnes:
    def __init__(self, omnes_function, amplitude, curve):
        args = omnes_function, amplitude, curve
        if isinstance(curve, Adaptive):
            self.func = make_curved_adaptive(*args)
        elif isinstance(curve, VectorDecay):
            self.func = make_curved_vector_decay(*args)
        elif isinstance(curve, Real):
            self.func = make_curved_real(*args)
        else:
            raise ValueError('unknown curve type')

        self.original = self.func.original

    def __call__(self, mandelstam_s):
        return self.func(mandelstam_s)
