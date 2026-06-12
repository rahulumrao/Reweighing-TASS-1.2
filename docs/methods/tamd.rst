TAMD / d-AFED
=============

Overview
--------

Temperature Accelerated Molecular Dynamics (TAMD) and its adiabatic variant
(d-AFED) accelerate the exploration of extended collective variables by coupling
them to a high-temperature bath while keeping the physical system at
:math:`T_0` :cite:`awasthi2017tass`.

For an extended CV :math:`\mathbf{s} = (s_1, s_2, \ldots)`, the effective
dynamics follow a Langevin equation at an elevated CV temperature :math:`T`:

.. math::

   M_\mathbf{s}\,\ddot{\mathbf{s}} =
   -\nabla_\mathbf{s} F(\mathbf{s}) -
   \gamma_\mathbf{s}\,\dot{\mathbf{s}} +
   \sqrt{2\gamma_\mathbf{s}\,k_B T}\,\boldsymbol{\eta}(t)

where :math:`M_\mathbf{s}` is the fictitious mass, :math:`\gamma_\mathbf{s}` the
friction, and :math:`F(\mathbf{s})` the free energy along the extended CVs.

Temperature reweighting
-----------------------

Because the CVs are sampled at temperature :math:`T` rather than the physical
temperature :math:`T_0`, each configuration must be reweighted by the factor:

.. math::

   w_T = \exp\!\left[
     \left(\frac{1}{k_B T_0} - \frac{1}{k_B T}\right)
     \left(F(\mathbf{s}) - V_\mathrm{bias}(\mathbf{s})\right)
   \right]

In practice, Reweighing-TASS applies the TASS temperature correction alongside
umbrella and metadynamics bias removal.  The relevant input parameters are:

* ``SYSTEM TEMP`` (:math:`T_0`) — physical system temperature
* ``CV TEMP`` (:math:`T`) — extended CV propagation temperature

Extended CVs in TASS
--------------------

In a TASS simulation, one CV is reserved for umbrella sampling (``UCV COLUMN``),
one (optionally) for metadynamics (``MTD CV COLUMN``), and the remaining CVs are
TAMD/d-AFED extended coordinates.  These "TASS CVs" are identified automatically
as all CVs that are neither the umbrella nor the metadynamics index.

The grid parameters for TASS CVs are set in the ``GRIDS`` section and used for
2D histogramming when ``PROBABILITY DIMENSION = 2`` and metadynamics is
disabled (US vs. TASS CV output in ``PROB_2D.dat``).

See also
--------

* :doc:`tass` — full TASS method combining TAMD, US, and MetaD
* :doc:`../input_reference` — ``SYSTEM TEMP`` and ``CV TEMP`` keywords
