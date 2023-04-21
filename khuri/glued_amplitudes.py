"""
===========================================================
Glue between amplitudes and high-energy phase prescriptions
===========================================================
"""
import numpy as np

from khuri.amplitude import from_phase


class GluedAmplitude:
    """Glue together an amplitude and a high-energy phase prescription."""
    def __init__(self, phase_decorator, low_energy_amplitude, transition):
        """
        Parameters
        ----------
        phase_decorator
            the decorator used to continue the phase and amplitude at values
            of Mandelstam s above `transition`
        low_energy_amplitude
            the amplitude used as long as the real part of Mandelstam s stays
            below `transition`
        transition:
            above this point, the high energy prescription sets in
        """
        self.transition = transition
        self.low_energy_amplitude = low_energy_amplitude
        self.phase = phase_decorator(self.low_energy_phase)
        self.amplitude = np.vectorize(self.amplitude_helper)

    def low_energy_phase(self, mandelstam_s):
        return np.angle(self.low_energy_amplitude(mandelstam_s))

    @from_phase(1.0)
    def high_energy_amplitude(self, mandelstam_s):
        return self.phase(mandelstam_s)

    def amplitude_helper(self, mandelstam_s):
        if mandelstam_s.real < self.transition:
            return self.low_energy_amplitude(mandelstam_s)
        assert mandelstam_s.imag < 1e-10
        return self.high_energy_amplitude(mandelstam_s.real)
