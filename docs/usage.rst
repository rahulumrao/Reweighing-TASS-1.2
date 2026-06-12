Usage
=====

Running the analysis
--------------------

The main executable reads control parameters from **standard input** (typically
redirected from ``input.tass``) and replica file paths from ``input.inp``.

Using the provided script
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   chmod 755 run.sh
   ./run.sh

The ``run.sh`` script creates an ``input.tass`` file and pipes it to
``tass_analysis.x``.

Direct invocation
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   bin/tass_analysis.x < input.tass

Or, if the executable is on ``PATH``:

.. code-block:: bash

   tass_analysis.x < input.tass

Required input files
--------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - File
     - Description
   * - ``input.tass``
     - Control parameters (keywords and values); read from stdin
   * - ``input.inp``
     - Per-replica umbrella parameters and paths to CV/MTD trajectory files
   * - CV trajectory files
     - One file per umbrella replica (paths listed in ``input.inp``)

``input.inp`` format
~~~~~~~~~~~~~~~~~~~~

For each umbrella window ``ir = 1 … N``:

**CPMD** (with metadynamics enabled):

.. code-block:: text

   <umbrella_mean>  <umbrella_kappa>
   <path/to/cvmdck_mtd>
   <path/to/parvar_mtd>
   <path/to/colvar_mtd>

**CPMD** (metadynamics disabled):

.. code-block:: text

   <umbrella_mean>  <umbrella_kappa>
   <path/to/cvmdck_mtd>

**PLUMED** (with metadynamics):

.. code-block:: text

   <umbrella_mean>  <umbrella_kappa>
   <path/to/COLVAR>
   <path/to/HILLS>

**PLUMED** (metadynamics disabled):

.. code-block:: text

   <umbrella_mean>  <umbrella_kappa>
   <path/to/COLVAR>

Units: CPMD umbrella force constants are in atomic units (converted internally
to kcal/mol); PLUMED force constants are in kJ/mol (converted to kcal/mol).

Workflow
--------

.. code-block:: text

   input.tass ──(stdin)──► tass_analysis.x ◄── input.inp
                                  ▲
                                  │ CV trajectory files
                                  │
                    ┌─────────────┴─────────────┐
                    │     REWEIGHTING TOOL      │
                    └─────────────┬─────────────┘
                          pmf   │   prob
                    ┌───────────┴───────────┐
                    ▼                       ▼
         free_energy.dat            whaminput + PROB*.dat
         free_energy_2D.dat                  │
                                             ▼
                                        wham.x (MPI)
                                             │
                                             ▼
                                    WHAM free energy

Output files
------------

Main analysis (``tass_analysis.x``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - Description
   * - ``cv.dat_<ir>``
     - Extracted CV time series for replica ``ir``
   * - ``free_energy.dat``
     - 1D PMF along the umbrella CV (mean-force method)
   * - ``free_energy_2D.dat``
     - 2D free energy surface (mean-force method)
   * - ``av_dfds.dat``
     - Average mean force (:math:`\partial F / \partial s`) at each umbrella center
   * - ``interp_free_energy.dat``
     - B-spline interpolated 1D free energy (when ``B-SPLINE INTERPOLATION`` is set)
   * - ``whaminput``
     - WHAM input file (generated when ``REWEIGHTING TOOL = prob``)
   * - ``PROB.dat_<ir>``
     - Unbiased 1D probability per replica
   * - ``Pu_2D.dat_<ir>``
     - Unbiased 2D probability (US vs. MTD CV)
   * - ``PROB_2D.dat``
     - Unbiased 2D probability (US vs. TASS temperature CV)
   * - ``variance.dat``
     - Block variance data (statistical error analysis)
   * - ``delta_G.dat``
     - Statistical error in free energy differences
   * - ``ct_test.dat_<ir>``
     - Well-tempered metadynamics :math:`c(t)` factor (debug output)
   * - ``vbias_test.dat_<ir>``
     - Accumulated metadynamics bias vs. time (debug output)

WHAM analysis (``wham.x``)
~~~~~~~~~~~~~~~~~~~~~~~~~~

Run after generating ``whaminput`` and ``cv.dat_*`` files:

.. code-block:: bash

   mpirun -np <N> bin/wham.x

WHAM reads ``whaminput`` and the per-replica ``cv.dat_*`` histogram files to
produce a converged multidimensional free energy estimate.

B-spline utilities
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   bin/1d_bspline.x    # interactive 1D interpolation
   bin/2d_bspline.x    # interactive 2D interpolation

Example calculation
-------------------

A complete worked example with bundled CPMD trajectory data, reference output,
and plotting scripts is provided in the ``example/`` directory.  See
:doc:`example` for the step-by-step tutorial.

Common errors
-------------

* **``input.inp doesn't exist``** — create ``input.inp`` with replica paths
  before running.
* **``CV grid out of range``** — a CV value falls outside the ``GRIDS`` bounds;
  adjust ``gridmin``/``gridmax`` or check trajectory data.
* **``NCV is NOT EQUAL TO GRIDS``** — the number of ``GRIDS`` lines must match
  ``NUMBER OF CV``.
* **Case sensitivity** — use uppercase keyword labels as documented in
  :doc:`input_reference`.
