Installation
============

Prerequisites
-------------

* **Fortran compiler**: GNU ``gfortran`` or Intel ``ifort`` / ``ifx``
* **GNU Make**
* **MPI Fortran compiler** (``mpif90`` or ``mpiifort``) — required only for
  building the WHAM executable
* **CMake** ≥ 3.10 (optional, for CMake-based builds)

Supported MD packages for input data: CPMD and PLUMED.

Download
--------

.. code-block:: bash

   git clone https://github.com/rahulumrao/Reweighing-TASS-1.2.git
   cd Reweighing-TASS-1.2

Makefile build (recommended)
----------------------------

The ``configure`` script generates a ``Makefile`` tailored to your Fortran
compiler.

.. code-block:: bash

   ./configure
   # When prompted, enter: gnu   (for gfortran)
   #                    or: intel (for ifort)

Build targets
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Target
     - Description
   * - ``make install``
     - Build ``tass_analysis.x`` (main analysis, serial)
   * - ``make bspline``
     - Build ``1d_bspline.x`` and ``2d_bspline.x`` interpolation tools
   * - ``make wham``
     - Build ``wham.x`` (MPI-parallel WHAM reweighting)
   * - ``make clean``
     - Remove object files and module files
   * - ``make distclean``
     - Remove ``bin/`` and ``lib/`` directories

Executables are installed to ``bin/``; object and module files go to ``lib/``.

Manual compiler selection
~~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, edit ``F90`` and ``FC`` directly in the ``Makefile``:

.. code-block:: makefile

   F90=ifort          # serial Fortran compiler
   FC=mpif90          # MPI wrapper for WHAM
   FCFLAGS = -fallow-argument-mismatch   # gfortran flag for legacy interfaces

GNU compiler notes
~~~~~~~~~~~~~~~~~~

For ``gfortran``, the flag ``-fallow-argument-mismatch`` is set by default to
handle legacy Fortran interfaces in the B-spline library.

Intel compiler notes
~~~~~~~~~~~~~~~~~~~~

Before building with Intel Fortran, source the oneAPI environment:

.. code-block:: bash

   source /opt/intel/oneapi/setvars.sh
   which ifort    # verify compiler is on PATH

MPI setup
~~~~~~~~~

WHAM requires an MPI Fortran compiler.  After installing OpenMPI_ or MPICH_,
add the library to your environment:

.. _OpenMPI: https://www.open-mpi.org/
.. _MPICH: https://www.mpich.org/

.. code-block:: bash

   export PATH=/path/to/mpi/bin:$PATH
   export LD_LIBRARY_PATH=/path/to/mpi/lib:$LD_LIBRARY_PATH

Then build WHAM:

.. code-block:: bash

   make wham

Serial vs. MPI executables
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Executable
     - Parallelism
     - Purpose
   * - ``tass_analysis.x``
     - Serial
     - Main TASS reweighting and mean-force analysis
   * - ``1d_bspline.x``, ``2d_bspline.x``
     - Serial
     - Standalone B-spline interpolation
   * - ``wham.x``
     - MPI
     - WHAM free energy reconstruction from ``whaminput``

CMake build (alternative)
-------------------------

A CMake build system can be used when ``CMakeLists.txt`` is present in the
repository root.  The workflow below follows the same pattern used in the
Reweighing-TASS development branch.

Serial build (GNU)
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   mkdir build && cd build
   cmake -DUSE_GNU_COMPILER=ON ..
   make
   make install

Serial build (Intel)
~~~~~~~~~~~~~~~~~~~~

Ensure ``ifort`` is sourced, then:

.. code-block:: bash

   mkdir build && cd build
   cmake -DUSE_INTEL_COMPILER=ON ..
   make
   make install

Custom installation prefix
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   cmake -DUSE_GNU_COMPILER=ON -DCMAKE_INSTALL_PREFIX=/opt/reweighing-tass ..
   make install

MPI / WHAM with CMake
~~~~~~~~~~~~~~~~~~~~~

When building with CMake, the WHAM target automatically searches for
``mpif90`` (GNU) or ``mpiifort`` (Intel) and links MPI libraries:

.. code-block:: bash

   cmake -DUSE_GNU_COMPILER=ON ..
   make wham.x
   make install

Verify the build
----------------

.. code-block:: bash

   ls bin/
   # tass_analysis.x  [1d_bspline.x  2d_bspline.x  wham.x]

Run ``./configure --help`` or ``./configure --info`` for additional information
about the program.

B-spline library
----------------

The package bundles `bspline-fortran
<https://github.com/jacobwilliams/bspline-fortran>`_ for spline interpolation.
Its object files are compiled automatically as dependencies of the main targets;
no separate installation is required.
