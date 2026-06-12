Umbrella Sampling
=================

Overview
--------

Umbrella sampling (US) enhances the sampling of a rare-event coordinate
:math:`s` by adding a harmonic biasing potential that restrains the system near
a series of predefined centers :math:`\{s_0^{(i)}\}`:

.. math::

   V_\mathrm{US}(s) = \frac{1}{2}\,\kappa\,(s - s_0)^2

where :math:`\kappa` is the force constant (spring constant) and :math:`s_0` is
the umbrella center for a given replica window.

In TASS, umbrella sampling is applied along one selected collective variable
(the **umbrella CV**, index ``UCV COLUMN``).  Each replica is simulated with a
different umbrella center, collectively spanning the relevant region of the free
energy landscape.

Role in Reweighing-TASS
-----------------------

During reweighting, the harmonic bias is removed from the sampled configurations.
For each umbrella window :math:`i`, the code reads:

* The umbrella center :math:`s_0^{(i)}` (``pcons``)
* The force constant :math:`\kappa^{(i)}` (``kcons``)
* The CV trajectory file

The instantaneous mean force along the umbrella coordinate is computed as:

.. math::

   \frac{\partial F}{\partial s}\bigg|_i =
   -\kappa^{(i)}\,\bigl(s - s_0^{(i)}\bigr)

This quantity is averaged over the production segment (``TMIN`` to ``TMAX``) and
integrated across windows to obtain the 1D PMF (see :doc:`meanforce`).

Input parameters
----------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``UCV COLUMN``
     - Index of the umbrella CV in trajectory files
   * - ``NUMBER OF UMBRELLA``
     - Total number of replica windows
   * - ``input.inp``
     - Per-replica ``<s_0>  <kappa>  <cv_file>`` entries

Unit conventions
----------------

Force constants in ``input.inp`` are converted internally:

* **CPMD**: given in atomic units, converted to kcal/mol via
  :math:`\kappa_\mathrm{kcal} = \kappa_\mathrm{a.u.} \times 627.51`
* **PLUMED**: given in kJ/mol, converted to kcal/mol via
  :math:`\kappa_\mathrm{kcal} = \kappa_\mathrm{kJ} \times 0.239006`

See also
--------

* :doc:`meanforce` — PMF integration from umbrella mean forces
* :doc:`wham` — alternative WHAM-based combination of umbrella windows
* :doc:`../input_reference` — full keyword reference
