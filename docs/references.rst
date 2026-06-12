References
==========

Primary literature
------------------

The TASS method and its reweighting approaches are described in the following
publications:

.. bibliography::
   :filter: cited

Key references
--------------

Awasthi & Nair (2017)
~~~~~~~~~~~~~~~~~~~~~

Shalini Awasthi and Narayanan N. Nair,
"Exploring high dimensional free energy landscapes: Temperature accelerated sliced sampling,"
*J. Chem. Phys.* **146**, 094108 (2017).
DOI: `10.1063/1.4977704 <https://doi.org/10.1063/1.4977704>`_

Introduces the TASS method combining TAMD, umbrella sampling, and metadynamics
for efficient exploration of high-dimensional CV spaces.

Awasthi & Nair (2019)
~~~~~~~~~~~~~~~~~~~~~

Shalini Awasthi and Narayanan N. Nair,
"Exploring high-dimensional free energy landscapes of chemical reactions,"
*WIREs Comput. Mol. Sci.* **9**, e1398 (2019).
DOI: `10.1002/wcms.1398 <https://doi.org/10.1002/wcms.1398>`_

Review of enhanced sampling methods including TASS and applications to chemical
reaction free energy landscapes.

Pal et al. (2021)
~~~~~~~~~~~~~~~~~

Ankit Pal, Subrata Pal, Shalini Verma, Motoyuki Shiga, and Narayanan N. Nair,
"Mean force based temperature accelerated sliced sampling: Efficient reconstruction of high dimensional free energy landscapes,"
*J. Comput. Chem.* **42**, 2010–2020 (2021).
DOI: `10.1002/jcc.26727 <https://doi.org/10.1002/jcc.26727>`_

Describes the mean-force reweighting algorithm implemented in Reweighing-TASS
for direct PMF reconstruction without WHAM.

Additional references cited in the methods documentation:

.. bibliography::
   :filter: citekey in {"flyvbjerg1989block", "barducci2008well"}

Software dependencies
---------------------

* `bspline-fortran <https://github.com/jacobwilliams/bspline-fortran>`_ —
  B-spline interpolation library by Jacob Williams

Citation
--------

If you use Reweighing-TASS 1.2 in your research, please cite the relevant
method papers above and acknowledge the software:

   Rahul Verma, *Reweighing-TASS* 1.2,
   https://github.com/rahulumrao/Reweighing-TASS-1.2
