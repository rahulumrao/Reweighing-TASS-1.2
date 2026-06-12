Methods
=======

Reweighing-TASS implements the post-processing steps required to reconstruct
unbiased free energy surfaces from biased TASS simulations.  The enhanced
sampling methods used during the MD simulation are described below; the
reweighting code removes their combined effect.

.. toctree::
   :maxdepth: 1

   umbrella
   metadynamics
   tamd
   tass
   meanforce
   wham

Analysis pipeline
-----------------

During a TASS simulation, three biasing mechanisms may act simultaneously:

1. **Umbrella sampling** — harmonic restraint along a selected CV
2. **Metadynamics** — history-dependent bias along a (possibly different) CV
3. **TAMD/d-AFED** — elevated temperature for extended CVs

Reweighing-TASS removes these biases in sequence and produces either:

* Unbiased **probability distributions** (``REWEIGHTING TOOL = prob``), or
* Direct **mean-force free energies** (``REWEIGHTING TOOL = pmf``)

For probability-based analysis, the WHAM algorithm further combines histograms
from all umbrella windows into a single free energy surface.
