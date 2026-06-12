WHAM
====

Weighted Histogram Analysis Method
----------------------------------

The Weighted Histogram Analysis Method (WHAM) combines biased probability
histograms from multiple umbrella windows into a single unbiased estimate of
the free energy :cite:`awasthi2017tass`.

In Reweighing-TASS, WHAM is used as the second step after probability
unbiasing: first ``tass_analysis.x`` generates unbiased histograms and a
``whaminput`` control file, then ``wham.x`` (MPI) performs the WHAM iteration.

Workflow
--------

1. Set ``REWEIGHTING TOOL = prob`` in ``input.tass``
2. Run the main analysis:

   .. code-block:: bash

      tass_analysis.x < input.tass

   This produces:

   * ``whaminput`` — WHAM control file
   * ``cv.dat_<ir>`` — per-replica CV histogram data
   * ``PROB.dat_<ir>`` — unbiased 1D probabilities

3. Run WHAM with MPI:

   .. code-block:: bash

      mpirun -np <N_procs> wham.x

WHAM equations
--------------

WHAM iteratively solves for the unbiased free energy :math:`F(\mathbf{q})` on a
discretized grid by minimizing the weighted sum of histogram deviations.  For
umbrella window :math:`i` with bias :math:`V_i(\mathbf{q})`, the unbiased
probability is:

.. math::

   P(\mathbf{q}) = \frac{
     \displaystyle\sum_i N_i \,
     \langle n_i(\mathbf{q}) \rangle \,
     e^{\beta F(\mathbf{q})}\, e^{-\beta V_i(\mathbf{q})}
   }{
     \displaystyle\sum_i N_i \,
     \langle n_i(\mathbf{q}) \rangle \,
     e^{-\beta V_i(\mathbf{q})}
   }

where :math:`N_i` is the number of frames in window :math:`i` and
:math:`\langle n_i(\mathbf{q}) \rangle` is the normalized histogram count.

The ``whaminput`` file format
-----------------------------

Generated automatically by ``tass_analysis.x``:

.. code-block:: text

   <tolerance>  <num_umbrella>
   <dimension>  <temperature_K>
   <gridmin_1>  <gridmax_1>  <griddif_1>
   <gridmin_2>  <gridmax_2>  <griddif_2>   # for 2D
   <s0_1>  <kappa_1>  <num_frames_1>
   <s0_2>  <kappa_2>  <num_frames_2>
   ...

MPI parallelization
-------------------

``wham.x`` is the only MPI-parallel component of Reweighing-TASS.  Build it
with:

.. code-block:: bash

   make wham

Ensure ``mpif90`` (or ``mpiifort`` for Intel) is available and MPI environment
variables are set (see :doc:`../installation`).

Comparison with mean-force method
---------------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 37 38

   * - Feature
     - Mean force (``pmf``)
     - WHAM (``prob`` + ``wham.x``)
   * - Parallelism
     - Serial
     - MPI
   * - Input
     - Direct force integration
     - Histogram iteration
   * - Dimensions
     - 1D and 2D
     - 1D and 2D
   * - MTD support
     - Yes
     - Yes (via unbiasing step)

See also
--------

* :doc:`meanforce` — direct PMF alternative
* :doc:`../usage` — output file descriptions
* :doc:`../installation` — MPI build instructions
