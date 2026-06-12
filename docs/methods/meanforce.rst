Mean Force Method
=================

Overview
--------

The mean-force (PMF) reweighting method reconstructs the free energy along
selected collective variables by computing and integrating the average force
experienced at each umbrella window :cite:`pal2021mftass`.

This approach bypasses explicit histogram construction and WHAM iteration,
providing a direct route to 1D and 2D free energy surfaces from TASS data.

1D free energy
--------------

For each umbrella replica :math:`i` centered at :math:`s_0^{(i)}`, the
instantaneous mean force along the umbrella CV is:

.. math::

   \left\langle \frac{\partial F}{\partial s} \right\rangle_i =
   -\kappa^{(i)} \left\langle s - s_0^{(i)} \right\rangle

where the average is taken over MD frames from ``TMIN`` to ``TMAX``.

When metadynamics is enabled, each frame is additionally reweighted by the
metadynamics bias:

.. math::

   \left\langle \frac{\partial F}{\partial s} \right\rangle_i =
   \frac{\displaystyle\sum_{t} \frac{\partial F}{\partial s}\,
   e^{\beta V_\mathrm{bias}(t)}}
   {\displaystyle\sum_{t} e^{\beta V_\mathrm{bias}(t)}}

The 1D PMF is obtained by integrating the averaged mean force across umbrella
windows:

.. math::

   F(s) = \int_{s_\mathrm{min}}^{s}
   \left\langle \frac{\partial F}{\partial s'} \right\rangle
   \mathrm{d}s'

Output is written to ``free_energy.dat`` and ``av_dfds.dat``.

2D free energy
--------------

For 2D PMF calculation (``PROBABILITY DIMENSION = 2``), the code histograms
configurations along a second CV (specified in ``PROBABILITY CV INDEX``) within
each umbrella window and computes the conditional free energy:

.. math::

   F(s_1, s_2) = -k_B T_0 \ln P(s_1, s_2) + C

where :math:`P(s_1, s_2)` is the reweighted joint probability.  Output is
written to ``free_energy_2D.dat``.

B-spline interpolation
----------------------

When the ``B-SPLINE INTERPOLATION`` flag is set, the 1D PMF is interpolated
onto a finer grid using the bundled `bspline-fortran` library.  The smoothed
profile is saved to ``interp_free_energy.dat``.

Usage
-----

Set in ``input.tass``:

.. code-block:: text

   REWEIGHTING TOOL
   pmf
   PROBABILITY DIMENSION
   2
   PROBABILITY CV INDEX
   1 3

Advantages
----------

* No MPI required — runs entirely in serial
* Direct integration avoids WHAM convergence issues
* Supports both 1D and 2D free energy surfaces
* Compatible with metadynamics bias removal

See also
--------

* :doc:`umbrella` — source of mean force data
* :doc:`wham` — alternative histogram-based approach
* :doc:`../input_reference` — ``REWEIGHTING TOOL = pmf``
