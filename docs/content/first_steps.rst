First steps
===========

Try the following inside a python script (or interactive session) to produce a
plot of the Omnes function of the Madrid p-wave:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt

    from khuri import madrid_global, phases, omnes


    THRESHOLD = (2.0 * madrid_global.PION_MASS)**2


    @phases.asymptotic1(matching_point=1.12**2)
    def phase(s):
        return madrid_global.p_wave_phase(s)


    omnes_function = omnes.generate_omnes(phase, threshold=THRESHOLD,
                                          constant=np.pi, cut=1e10)

    energies = np.linspace(0, 1.2, 200)
    omnes_values = omnes_function(energies**2)

    plt.title('The Omnes function of the Madrid p-wave')
    plt.plot(energies, np.real(omnes_values), label='Re')
    plt.plot(energies, np.imag(omnes_values), label='Im')
    plt.xlabel('E/GeV')
    plt.legend()
    plt.show()
