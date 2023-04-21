Khuri Treiman Pinocchio
=======================

The following lines give a step by step explanation regarding the C++ implementation of the iterative Khuri-Treiman solver relying on the Pinocchio-like deformation of the complex integration contour.

A detailed documentation of the C++ code can be done via :ref:`doxygen <doc>`.
This documentation serves the purpose to provide a broader picture of the code structure and explain what each of the relevant files is used for.

All masses, energies and momenta are used in units of the charged pion mass :math:`M_\pi=0.13957018\, \text{GeV}`. Particle masses can be found in ``constants``. 

All interpolations rely on cubic splines. 

Useful algebraic functions/operations (e.g. derivatives, integrals, interpolations) can be found in ``gsl_interface``, while some generalizations for complex functions are included in ``cauchy``.

As we are trying to solve a system of coupled integral equations with an iterative solution algorithm, it is vital for the runtime to interpolate each non-trivial function and integral that are used more than once. For this purpose we often provide first a class evaluating the amplitudes/integrals along predefined lists which are subsequently interpolated in the constructor of a following class. This ensures that the evaluation and interpolation of each integral is only computed once.

In general we define one parent class for each of the following subproblems, while process specific method definitons are relegated to the derived classes.

The naming scheme is as following:

* ``Eta3Pi``: :math:`\eta\to 3\pi` (including the C-violating processes, where specific classes for these are marked by an additional ``C``)
* ``Etap3Pi``: :math:`\eta'\to 3\pi` (including the C-violating processes, where specific classes for these are marked by an additional ``C``)
* ``EtapEtaPiPi``: :math:`\eta'\to \eta\pi\pi` (including the C-violating processes, where specific classes for these are marked by an additional ``C``)
* ``V3Pi``: :math:`V\to 3\pi` :math:`(J^{PC}=1^{--})`
* ``X3Pi``: :math:`X\to 3\pi` :math:`(J^{PC}=1^{-+})`
* ``T3Pi``: :math:`T\to 3\pi` :math:`(J^{PC}=2^{++})`

The final classes from ``iterative_solution`` are exported via pybind11 in ``bindings/khuri_pinocchio_bindings.cpp``. If a new process is implemented the corresponding binding needs to be added there.

Side remark for the process e:math:`\eta'\to \eta\pi\pi`: it is not possible to distinguish all amplitudes unambiguously by their isospin. Therefore replace the isospin in this process by the combination of three quantum numbers, e.g. use 'isospin=110' for a state with total isospin of final state =1, isospin of two-particle state =1 and angular momentum =0. Additionally one has to distinguish between the s-channel (:math:`\eta'\eta\to \pi\pi`) and the t-channel (:math:`\eta'\pi\to \eta\pi`) which have different thresholds.


Lists & Splines
---------------
* ``enums``: Use the 'enum' command to use several switches between the different implemented processes, the subtraction constants defining the basis solutions, the curve parameters describing the complex integration contour, ...
* ``array``: First create a one-dimensional grid in Mandelstam variables along which each amplitude and integral is evaluated and interpolated. Use a narrow spacing within the physical region and more points around scattering thresholds (and thresholds resulting from crossing symmetry). Often one has to fine tune the arrays by hand for each problem to avoid numerical instabilities (for instance by omitting single points causing trouble). 
* ``phase``: This file contains functions to read in scattering phase shifts. Using pre-calculated tables of phases as input, an analytic interpolation of those (with custom continuation to large energies) is returned as output.
* ``asymptotics``: Provides different continuation functions for the phase towards large energies.
* ``splined_omnes``: Calculate and interpolate the Omnès function (defined in ``omnes``) for Mandelstam s infinitesimally above and below the unitarity cut as well as for different curve parameters needed to compute the angular averages.
* ``matching``: The analytic continuation of a :math:`2\to 2` scattering process to the three-body decay introduces poles in the dispersion integral at scattering thresholds (and thresholds resulting from crossing symmetry). To avoid these poles, which scale like :math:`\sqrt{s}` for S-waves, :math:`\sqrt{s}^3` for P-waves and :math:`\sqrt{s}^5` for D-waves, an expansion of the integrand is matched to the exact values. Using these expansions the dispersion integrals can be rewritten such that the poles are shifted out of the integration limits. The matching process has the freedom to adjust the matching point ``epsilon`` and the range of validity ``validity`` within which the expansion is used. Note that there are two types of matching: Matching for the Tilde function, which is denoted by an additional ``_tilde``, and when calculating the disperion integral. These four parameters need to be adjusted for each process and decay mass to resolve instabilities around the thresholds.

Angular Averages 
----------------
* ``path``: Provides a parameterization for the complex integration contour in the angular averages, aka Pinocchio path, for a three-body final state in which all particles have the same mass (e.g. :math:`3\pi` final states). This parameterization is appealing because only one curve parameter is needed to describe the egg-like shape of the curve.
* ``splined_path``: Interpolation of the contour (and its derivative) from ``path``. This speeds up the code.
* ``path_eta_pi_pi``: A more general parameterization than the one given in ``path``, which works for two distinct masses in the final state (e.g. in :math:`\eta\pi\pi` final states).
* ``angular_average``: Define and evaluate the angular averages (given by partial wave projections) along the lists given in ``array``. Subsequently interpolate these lists and finally multiply the angular average with a factor resulting from ``matching`` to directly cancel all poles that cause trouble in the dispersion integral (except for the pseudo-threshold, this will be handled in ``dispersion_integral``). If you would like to change the reconstruction theorem you have to adjust the angular averages accordingly.

Dispersion Integral
-------------------
* ``dispersion_integral``: First apply the matched expressions to the integrand of the dispersion integral to avoid numerical issues close to the singularity at pseudo threshold. If you would like to demand a different asymptotic behavior of phase shifts and/or amplitudes you have to adjust the integrands accordingly. The dispersion integral is decomposed into integrals containing the singularity at pseudo threshold, the cauchy singularity as well as analytical expressions. Thereby the poles of each integral are shifted out of the integration limits. Similar as in ``splined_omnes``, the dispersion integral is evaluated along the predefined lists from ``array`` for Mandelstam s infinitesimally above and below the unitarity cut as well as for different curve parameters needed to compute the angular averages and interpolated afterwards.


Solving for single variable functions
-------------------------------------
* ``basis_function``: Defines 'basis solutions' for the angular averages and the single variable functions. Each subtraction constants has its own set of basis solutions. Additionally basis solutions to initiate the iteration procedure are provided by setting the angular averages to zero (known as homogeneous solution).
* ``iterative_solution``: First initiate Omnès function and homogeneous basis solutions. Then calculate the angular averages, from those the single variable functions and with these again the angular averages. Iterate this procedure until the desired convergence is achieved (vary the number of iterations in the constructor). The results of each iteration will be stored in lists and can be called separately if needed. Furthermore the final solutions are exported as data files.


KT solver
---------
* ``example/KT``: Choose your decay of interest, play around with the parameters to smoothen the amplitude close scattering thresholds ('epsilon_tilde', 'validity_tilde') and pseudo-threshold ('epsilon', 'validity'). Empirically, good values for all of these parameters are between 0.1 and 1. Choose the number of iterations. The number of necessary iterations increases with the phase space and the amount of subtraction constants (iterative solution will fail if one of both is too large; in this case a matrix inversion should be used).


Pybind11
--------
The same thing can be done directly in python using the bindings via pybind11.
An example of to test the functions in a python script is provided below.

.. code-block:: python

		import os.path
		import numpy as np
		import matplotlib.pyplot as plt

		import khuri.khuri_pinocchio as kt_pino


		virtuality = 0.78266 / 0.13957

		directory = './phase/' # put in the directory where the datafile for the phase shift is

		path_phase_1 = os.path.join(directory, 'phase_pipi_1.dat')

		output_directory = os.path.join('./', '1--') # make sure this directory exists

		path_output_a1 = os.path.join(output_directory, 'output_a1_'+str(virtuality)+'.txt')
		path_output_b0 = os.path.join(output_directory, 'output_b0_'+str(virtuality)+'.txt')
				
		kt_pino_config = {
		    'epsilon': 0.1,
		    'validity': 0.05,
		    'epsilon_tilde': 0.7,
		    'validity_tilde': 0.1,
		    'cutoff': 1000,
		    'iterations': 7,
		    'max_subs': 1,
		}
		result = kt_pino.IterationV3Pi(1., virtuality,
		                             path_phase_1, kt_pino_config['epsilon'],
		                             kt_pino_config['validity'], kt_pino_config['epsilon_tilde'],
		                             kt_pino_config['validity_tilde'], kt_pino_config['cutoff'],
		                             kt_pino_config['iterations'], kt_pino_config['max_subs'],
		                             path_output_a1, path_output_b0)

		energies = np.linspace(0, 100, 500)
		kt_result = result(energies, isospin=1, sub=kt_pino.SubtractionConstant.a1,
		            	set=kt_pino.Setting.above, iteration=7)

		plt.title('Example for KT-equations')
		plt.plot(energies, np.real(kt_result), label='Re')
		plt.plot(energies, np.imag(kt_result), label='Im')
		plt.xlabel('E/GeV')
		plt.legend()
		plt.show()