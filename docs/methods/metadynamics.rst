Metadynamics
============

Overview
--------

Metadynamics is a non-equilibrium enhanced sampling method that fills the free
energy wells of a collective variable :math:`\xi` by depositing history-dependent
Gaussian hills:

.. math::

   V_\mathrm{MTD}(\xi, t) =
   \sum_{t_i < t}\,
   h \exp\!\left[
     -\frac{(\xi - \xi_i)^2}{2w^2}
   \right]

where :math:`h` is the hill height, :math:`w` the hill width, and
:math:`\xi_i` the CV value at deposition time :math:`t_i`.

Well-tempered metadynamics
--------------------------

Reweighing-TASS supports **well-tempered metadynamics** (WT-MetaD), in which
the hill heights are scaled by a time-dependent factor :math:`c(t)` to ensure
convergence to the correct free energy :cite:`barducci2008well`:

.. math::

   V_\mathrm{WT}(\xi, t) =
   c(t)\,V_\mathrm{MTD}(\xi, t)

The factor :math:`c(t)` depends on the bias factor :math:`\gamma` (keyword
``MTD BIAS FACTOR``) and the extended CV temperature :math:`T`:

.. math::

   c(t) = e^{-\,V_\mathrm{MTD}(\xi, t)\,/\,(\gamma - 1)\,k_B T}

Bias removal in Reweighing-TASS
-------------------------------

When ``MTD ENABLED = y``, the code:

1. Reads hill deposition data from CPMD (``parvar_mtd``, ``colvar_mtd``) or
   PLUMED (``HILLS``) files
2. Reconstructs the time-dependent bias :math:`V_\mathrm{bias}(t)` and the
   WT factor :math:`c(t)`
3. Computes the **unbiased** probability along the MTD CV by reweighting each
   frame with :math:`e^{\beta V_\mathrm{bias}(t)}`

The reweighting uses the CV print frequency (``CV PRINT FREQUENCY``) and hill
deposition frequency (``MTD PRINT FREQUENCY``) to synchronize MD and MTD time
steps.

Input parameters
----------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``MTD ENABLED``
     - ``y`` to enable metadynamics bias removal
   * - ``MTD CV COLUMN``
     - Index of the metadynamics CV
   * - ``MTD BIAS FACTOR``
     - Well-tempered bias factor :math:`\gamma` (K)
   * - ``MTD PRINT FREQUENCY``
     - MD steps between hill depositions
   * - ``CV PRINT FREQUENCY``
     - MD steps between CV output frames

Output
------

When metadynamics is enabled with 2D probability analysis, unbiased
US-vs-MTD distributions are written to ``Pu_2D.dat_<ir>`` for each replica.

See also
--------

* :doc:`tass` — how metadynamics fits into the TASS framework
* :doc:`../input_reference` — keyword reference
