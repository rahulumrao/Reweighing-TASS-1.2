Input reference
===============

Control parameters are supplied via **standard input**, typically from a file
named ``input.tass``.  Each keyword occupies its own line in **uppercase**,
followed by one or more lines containing the corresponding value(s).

.. important::

   Keyword labels are matched case-sensitively.  Use the exact uppercase spellings
   shown below.  Optional boolean flags (e.g. ``B-SPLINE INTERPOLATION``) are
   activated simply by including the keyword line with no value.

Replica-specific data (umbrella centers, force constants, file paths) are read
from the separate file ``input.inp`` — see :doc:`usage`.

Summary table
-------------

.. list-table::
   :header-rows: 1
   :widths: 28 10 12 42

   * - Keyword
     - Type
     - Default
     - Description
   * - ``NUMBER OF CV``
     - integer
     - —
     - Total number of collective variables in the simulation
   * - ``CODE NAME``
     - string
     - —
     - MD package: ``CPMD`` or ``PLUMED``
   * - ``NUMBER OF UMBRELLA``
     - integer
     - —
     - Number of umbrella replica windows
   * - ``UCV COLUMN``
     - integer
     - —
     - Column index of the umbrella (US) CV in trajectory files
   * - ``MTD ENABLED``
     - ``y``/``n``
     - ``n``
     - Whether well-tempered metadynamics was used
   * - ``MTD CV COLUMN``
     - integer
     - ``0``
     - Column index of the metadynamics CV (required if MTD enabled)
   * - ``MTD BIAS FACTOR``
     - real
     - ``1500``
     - Well-tempered metadynamics bias factor (K)
   * - ``SYSTEM TEMP``
     - real
     - ``300``
     - Physical system temperature :math:`T_0` (K)
   * - ``CV TEMP``
     - real
     - ``300``
     - Extended CV (TAMD) temperature :math:`T` (K)
   * - ``TMIN``
     - integer
     - ``1``
     - First MD step included in analysis
   * - ``TMAX``
     - integer
     - ``1000``
     - Last MD step included in analysis
   * - ``REWEIGHTING TOOL``
     - string
     - ``pmf``
     - Analysis mode: ``pmf`` (mean force) or ``prob`` (probability + WHAM)
   * - ``PROBABILITY DIMENSION``
     - integer
     - —
     - Dimensionality of output: ``1`` or ``2`` (max 3)
   * - ``PROBABILITY CV INDEX``
     - int list
     - —
     - CV column indices for probability/PMF output
   * - ``STATISTICAL ERRORS BLOCK SIZE``
     - integer
     - off
     - Block size for statistical error estimation
   * - ``CV PRINT FREQUENCY``
     - integer
     - ``0``
     - Print frequency of CV data in trajectory files
   * - ``MTD PRINT FREQUENCY``
     - integer
     - ``0``
     - Hill deposition frequency in metadynamics files
   * - ``GRIDS``
     - real × 3 × NCV
     - —
     - Grid bounds and spacing for each CV
   * - ``B-SPLINE INTERPOLATION``
     - flag
     - off
     - Enable B-spline interpolation of 1D PMF
   * - ``READ CT``
     - flag
     - off
     - Read pre-computed :math:`c(t)` values from file
   * - ``READ VBIAS``
     - flag
     - off
     - Read pre-computed bias potential from file

Detailed keyword reference
--------------------------

NUMBER OF CV
~~~~~~~~~~~~

:Type: integer
:Required: yes
:Example:

.. code-block:: text

   NUMBER OF CV
   5

Total number of collective variables tracked during the TASS simulation.
Must match the number of CV columns in trajectory files and the number of
``GRIDS`` entries.

CODE NAME
~~~~~~~~~

:Type: string (``CPMD`` | ``PLUMED``)
:Required: yes
:Example:

.. code-block:: text

   CODE NAME
   CPMD

Selects the MD package format for reading trajectory and metadynamics files.
Unit conversions for umbrella force constants differ between packages (see
:doc:`usage`).

NUMBER OF UMBRELLA
~~~~~~~~~~~~~~~~~~

:Type: integer
:Required: yes
:Example:

.. code-block:: text

   NUMBER OF UMBRELLA
   22

Number of umbrella sampling replica windows.  Must match the number of entries
in ``input.inp``.

UCV COLUMN
~~~~~~~~~~

:Type: integer (1-based index)
:Required: yes
:Example:

.. code-block:: text

   UCV COLUMN
   1

Column index of the umbrella sampling CV in the trajectory file.  For
probability reweighting, at least one entry in ``PROBABILITY CV INDEX`` must
equal this value.

MTD ENABLED
~~~~~~~~~~~

:Type: ``y`` or ``n``
:Default: ``n``
:Example:

.. code-block:: text

   MTD ENABLED
   y

Enable metadynamics bias removal.  When ``y``, metadynamics trajectory files
must be listed in ``input.inp``.

MTD CV COLUMN
~~~~~~~~~~~~~

:Type: integer
:Default: ``0``
:Example:

.. code-block:: text

   MTD CV COLUMN
   2

Column index of the metadynamics CV.  Required (non-zero) when ``MTD ENABLED``
is ``y``.

MTD BIAS FACTOR
~~~~~~~~~~~~~~~

:Type: real (K)
:Default: ``1500``
:Example:

.. code-block:: text

   MTD BIAS FACTOR
   1200

Well-tempered metadynamics bias factor :math:`\gamma` used to compute the
time-dependent factor :math:`c(t)` :cite:`barducci2008well`.

SYSTEM TEMP
~~~~~~~~~~~

:Type: real (K)
:Default: ``300``
:Example:

.. code-block:: text

   SYSTEM TEMP
   300

Physical temperature :math:`T_0` of the system (not the extended CV
temperature).

CV TEMP
~~~~~~~

:Type: real (K)
:Default: ``300``
:Example:

.. code-block:: text

   CV TEMP
   1000

Temperature :math:`T` at which the extended CVs are propagated in the TAMD
/d-AFED scheme.  Used in the TASS temperature reweighting factor.

TMIN
~~~~

:Type: integer (MD step)
:Default: ``1``
:Example:

.. code-block:: text

   TMIN
   1000

First molecular dynamics step included in the analysis (after equilibration).

TMAX
~~~~

:Type: integer (MD step)
:Default: ``1000``
:Example:

.. code-block:: text

   TMAX
   14000

Last molecular dynamics step included in the analysis.  If omitted, all
available steps are used.

REWEIGHTING TOOL
~~~~~~~~~~~~~~~~

:Type: string (``pmf`` | ``prob``)
:Default: ``pmf``
:Example:

.. code-block:: text

   REWEIGHTING TOOL
   pmf

Selects the analysis pathway:

* ``pmf`` — compute free energy directly via the mean-force method
  (:doc:`methods/meanforce`)
* ``prob`` — generate unbiased probability histograms and a ``whaminput`` file
  for subsequent WHAM analysis (:doc:`methods/wham`)

PROBABILITY DIMENSION
~~~~~~~~~~~~~~~~~~~~~

:Type: integer (``1`` or ``2``)
:Required: when using ``prob`` or for 2D PMF
:Example:

.. code-block:: text

   PROBABILITY DIMENSION
   2

Number of dimensions for the output probability or free energy surface.

PROBABILITY CV INDEX
~~~~~~~~~~~~~~~~~~~~

:Type: space-separated integers
:Required: yes
:Example:

.. code-block:: text

   PROBABILITY CV INDEX
   1 3

CV column indices along which the probability or PMF is computed.  Provide one
index for 1D output, two for 2D output.

STATISTICAL ERRORS BLOCK SIZE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:Type: integer
:Default: disabled
:Example:

.. code-block:: text

   STATISTICAL ERRORS BLOCK SIZE
   4

When set, performs block averaging to estimate statistical errors in the mean
force and writes ``variance.dat`` and ``delta_G.dat``.  A block size of 4–5 is
typically optimal :cite:`flyvbjerg1989block`.

CV PRINT FREQUENCY
~~~~~~~~~~~~~~~~~~

:Type: integer
:Default: ``0``
:Example:

.. code-block:: text

   CV PRINT FREQUENCY
   1

Frequency (in MD steps) at which CV values are written to the trajectory file
(``cvmdck_mtd`` for CPMD, ``COLVAR`` for PLUMED).  Must match the simulation
output frequency.

MTD PRINT FREQUENCY
~~~~~~~~~~~~~~~~~~~

:Type: integer
:Default: ``0``
:Example:

.. code-block:: text

   MTD PRINT FREQUENCY
   10

Frequency (in MD steps) at which Gaussian hills are deposited in the
metadynamics output files.

GRIDS
~~~~~

:Type: three reals per CV line
:Required: yes
:Example:

.. code-block:: text

   GRIDS
   1.8 6.0  0.01
   0.5 5.0  0.05
   1.0 8.0  0.01
   1.0 10.0 0.01
   1.0 8.0  0.01

For each of the ``NUMBER OF CV`` collective variables, specify:

* ``gridmin`` — lower bound of the histogram grid
* ``gridmax`` — upper bound of the histogram grid
* ``griddif`` — bin width

The number of bins is computed as
:math:`N_\mathrm{bins} = \mathrm{NINT}\bigl((\mathrm{gridmax} - \mathrm{gridmin}) / \mathrm{griddif}\bigr) + 1`.

B-SPLINE INTERPOLATION
~~~~~~~~~~~~~~~~~~~~~~

:Type: flag (no value)
:Default: disabled
:Example:

.. code-block:: text

   B-SPLINE INTERPOLATION

When present, the 1D mean-force PMF is interpolated onto a finer grid using
B-splines.  Output is written to ``interp_free_energy.dat``.

READ CT
~~~~~~~

:Type: flag (no value)
:Default: disabled

Read pre-computed well-tempered metadynamics :math:`c(t)` values from file
instead of recalculating them.

READ VBIAS
~~~~~~~~~~

:Type: flag (no value)
:Default: disabled

Read pre-computed metadynamics bias potential from file instead of
recalculating it.

Complete example
----------------

.. code-block:: text

   NUMBER OF CV
   5
   CODE NAME
   CPMD
   NUMBER OF UMBRELLA
   22
   UCV COLUMN
   1
   MTD ENABLED
   n
   MTD CV COLUMN
   2
   MTD BIAS FACTOR
   1200
   SYSTEM TEMP
   300
   CV TEMP
   1000
   TMIN
   1000
   TMAX
   14000
   REWEIGHTING TOOL
   pmf
   PROBABILITY DIMENSION
   2
   PROBABILITY CV INDEX
   1 3
   STATISTICAL ERRORS BLOCK SIZE
   4
   CV PRINT FREQUENCY
   1
   MTD PRINT FREQUENCY
   10
   GRIDS
   1.8 6.0  0.01
   0.5 5.0  0.05
   1.0 8.0  0.01
   1.0 10.0 0.01
   1.0 8.0  0.01

See also ``input.tass`` in the repository root and ``example/input.tass`` for
working examples paired with ``input.inp``.
