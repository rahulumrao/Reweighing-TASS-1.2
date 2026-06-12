Reweighing-TASS 1.2
====================

**Reweighing-TASS** is a modular Fortran program for post-processing output from
Temperature Accelerated Sliced Sampling (TASS) simulations performed with
`CPMD <https://www.cpmd.org/>`_ or `PLUMED <https://www.plumed.org/>`_.
It removes biases introduced by umbrella sampling, metadynamics, and temperature
acceleration, and reconstructs multidimensional free energy surfaces.

Overview
--------

Temperature Accelerated Sliced Sampling (TASS) combines temperature-accelerated
molecular dynamics (TAMD/d-AFED), umbrella sampling, and metadynamics to explore
high-dimensional collective variable (CV) spaces efficiently
:cite:`awasthi2017tass,awasthi2019wire`.  This reweighting package processes
trajectory data from each umbrella replica and computes:

* **1D and 2D unbiased probability distributions** via histogram reweighting
* **1D and 2D free energy surfaces** via the mean-force (PMF) method
  :cite:`pal2021mftass`
* **WHAM-based free energy reconstruction** from unbiased probabilities (MPI)
* **B-spline interpolation** of free energy profiles for smooth visualization

Features
--------

* Supports **CPMD** and **PLUMED** output formats
* Optional **well-tempered metadynamics** bias removal
* **Mean-force reweighting** for direct PMF calculation along user-defined CVs
* **Probability unbiasing** with automatic ``whaminput`` generation for WHAM
* **Statistical error estimation** via block averaging
  :cite:`flyvbjerg1989block`
* Standalone **1D/2D B-spline** interpolation utilities
* Serial execution for analysis; **MPI parallelization** for WHAM only

.. important::

   Input keyword labels in ``input.tass`` are **case-sensitive**.  Optional
   flags such as ``B-SPLINE INTERPOLATION`` are activated by the presence of
   the keyword line (no value required).

Quick start
-----------

.. code-block:: bash

   ./configure          # select gnu or intel compiler
   make install         # build tass_analysis.x
   cd example && ../bin/tass_analysis.x < input.tass

See :doc:`example` for a full walkthrough with the bundled data set.

Documentation contents
----------------------

.. toctree::
   :maxdepth: 2
   :caption: User guide

   installation
   usage
   example
   input_reference

.. toctree::
   :maxdepth: 2
   :caption: Methods

   methods/index

.. toctree::
   :maxdepth: 1
   :caption: References

   references

Author
------

**Rahul Verma**
Department of Chemistry, IIT Kanpur, India
Email: `vrahul@iitk.ac.in <mailto:vrahul@iitk.ac.in>`_

Original TASS reweighting concepts by Shalini Awasthi (``ashalini@iitk.ac.in``).

.. raw:: html

   <div class="footer-authorship">
     Rahul Verma &middot; Department of Chemistry &middot; IIT Kanpur
   </div>
