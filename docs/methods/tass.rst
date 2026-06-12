TASS
====

Temperature Accelerated Sliced Sampling
---------------------------------------

Temperature Accelerated Sliced Sampling (TASS) is an enhanced sampling method
that combines three techniques to efficiently map high-dimensional free energy
landscapes :cite:`awasthi2017tass,awasthi2019wire`:

1. **Umbrella sampling** — partitions the umbrella CV into overlapping windows
2. **Metadynamics** — fills free energy wells along a (possibly different) CV
3. **TAMD/d-AFED** — accelerates exploration of additional extended CVs at
   elevated temperature

The "sliced" nature of TASS refers to the simultaneous exploration of multiple
CV dimensions: the umbrella CV defines the slice center, while TAMD drives
rapid diffusion within each slice.

Method schematic
----------------

.. code-block:: text

   Biasing (simulation)              Reweighing (post-processing)
   ─────────────────────             ────────────────────────────
   Umbrella:  V_US = ½κ(s−s₀)²  ──►  Remove US bias
   MetaD:     V_MTD(ξ, t)       ──►  Remove MTD bias
   TAMD:      coupling at T>T₀  ──►  Temperature reweighting
                                           │
                                           ▼
                                    PMF or WHAM free energy

Combined bias
-------------

During a TASS simulation, the total biasing potential acting on the system is:

.. math::

   V_\mathrm{total} =
   V_\mathrm{US}(s) +
   V_\mathrm{MTD}(\xi, t) +
   V_\mathrm{coupling}(\mathbf{s})

where :math:`s` is the umbrella CV, :math:`\xi` the metadynamics CV, and
:math:`\mathbf{s}` the vector of extended TAMD CVs.

Reweighing-TASS removes these biases in post-processing to recover the
unbiased probability or free energy along user-selected CV combinations.

Typical workflow
----------------

1. **Simulation** — run TASS in CPMD or PLUMED with multiple umbrella replicas
2. **Prepare inputs** — create ``input.tass`` and ``input.inp`` with replica
   paths and analysis parameters
3. **Reweight** — run ``tass_analysis.x`` with ``REWEIGHTING TOOL = pmf`` or
   ``prob``
4. **Post-process** — optionally run WHAM (``wham.x``) or B-spline interpolation

Supported configurations
------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Configuration
     - Description
   * - US only
     - ``MTD ENABLED = n``; mean-force PMF or WHAM
   * - US + MTD
     - Full TASS with metadynamics; 2D US-vs-MTD probabilities
   * - US + TAMD
     - ``MTD ENABLED = n``; 2D US-vs-TASS-CV probabilities
   * - US + MTD + TAMD
     - Full TASS; mean-force or probability reweighting

See also
--------

* :doc:`umbrella`, :doc:`metadynamics`, :doc:`tamd` — component methods
* :doc:`meanforce` — direct PMF reconstruction
* :doc:`wham` — WHAM-based free energy estimation
