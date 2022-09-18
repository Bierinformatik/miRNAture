========
Appendix
========

Pre-processing miRBase data:
----------------------------

The ``MIRfix`` pipeline :cite:p:`Yazbeck:19a` provides the general core of
functions to curate *bona fide* metazoan microRNA annotations. To make use of
this curation process, it is fundamental to organize the input data in a
specific format, as referenced in more detail in :cite:`Yazbeck:19a`. In
summary, it is required:

 - A set of precursor sequences with their associated mature sequences.
 - Genome sequences from which miRNAs were annotated.
 - A relation file that describes the relation between precursors and their
   annotated matures.

Additional parameters are required, but they did not depend from external
information/databases.

On ``miRNAture`` the source of the curation data has been obtained from a re-evaluation
of the annotations deposited on ``miRBase`` v.22.1 :cite:`Kozomara:19`. In this
version, ``miRBase`` accounted for a set of X *canonical* and *non-canonical* miRNA families,
from which Y are constituted by metazoan sequences. Internally, ``miRNAture``
performs an evaluation of the *canonical* model, that relies on the correct
positioning cleavages performed by Drosha and Dicer at the precursor maturation steps. 
Computationally, this is translated on a correct position of the mature sequence, 
accurate delimitation of precursor, and a phylogenetic support, addressed by the 
construction of family *mature-anchored* structural alignments. As previously
reported for the ``Rfam`` miRNA families :cite:`VelandiaDiss:2022`, an iterative
assessment involves a selection of sequences, consistent criteria to evaluate
the miRNAs and their mature products, and generation of probabilistic models
derived from anchored-alignments to search additional candidates that would
incorporate defined curation criteria. This criteria was inherited to perform an 
evaluation of the ``miRBase`` metazoan families and generate the corrected
dataset that ``miRNAture`` uses to evaluate new candidates and their maturation
entities.

As a toy example, the family miR-17 (MIPF0000001) was selected to demonstrate
the assessment steps performed over all pre-calculated dataset used by `miRNAture`.
As reported in ``miRBase`` this family is composed by 239 miRNA precursors derived from 
39 vertebrate species. Through the filtering approach the following subsetting steps are
considered:

- Remove non-metazoan sequences.
- Filter duplicates (which share 100% identity) and select one representant
  sequence.

In this family, 117 duplicated sequences where recognized. For instance the
sequence bta-mir-18a (MI0004740) from *Bos taurus* has shown 21 orthologs, as follows:

.. _table_dup:

.. figure:: ./table_MI0004740.png
   :width: 600
   :align: center
   :alt: Orthologs in other species
   :figclass: align-center

   Identified orthologs from bta-mir-18a on vertebrates.

And some corresponding alignments:

.. _align_dup:

.. figure:: ./alignments_MI0004740.png
   :width: 600
   :align: center
   :alt: Orthologs in other species
   :figclass: align-center
   
   Alignments as evidence of 100% identity.

Remaining 122 families were subject of a structural assessment by ``MIRfix``,
which filtered 4 sequences based on the incorrect miRNA folding in regard their
annotated mature sequences, and one sequence contained a bad positioned mature
sequence in the reported precursor, a successful extension of the precursor based
on the miR and miR* prediction, rescued the candidate.

================================ ==============================================
  Category                           Accession numbers
================================ ==============================================
Bad position mature sequences    MI0004822
Filtered sequences               MI0012797, MI0012947, MI0019542, and MI0013837
================================ ==============================================

At the end of the assessment 118 sequences passed all filters to be considered
into the curation dataset used on ``miRNAture``. 

The same approach curated all metazoan miRNA families from ``miRBase`` (1415),
validating about 79% (1111) of the families and setting the curation dataset used on
``miRNAture``.


Construction of Hidden Markov and Covariance Models:
---------------------------------------------------

As described in :cite:`Velandia:2021`, a set of quality-filtering steps could be
used to construct family structural alignments and their corresponding
covariance models (CMs). In this case, to build new structural alignments from
``miRBase`` sequences, we selected all sequences from metazoan
species and removed all of those from studied organisms. Given that curated
subset, a genetic algorithm was used to maximize the quality the final structural
alignment. To do so, filtering miRNA sequences was done in function of: Identity
percentage (:math:`I`), phylogenetic distribution of sequences (:math:`C`) and quality
(:math:`Q`) [#f3]_, where: :math:`I = (60, 70, 80, 90, 100)`, :math:`C =` (Metazoa, Vertebrata,
Mammalia, Primates) and :math:`Q =` (normal, high). An individual :math:`A_{n}` was
defined as a vector :math:`\overrightarrow{A_{n}}= \begin{pmatrix} I \\ C \\ Q
\end{pmatrix}`, which return a structural alignment using ``MIRfix``, using selected
sequences. The *fitness* function (:math:`F`) to be maximized was defined through
empirical observation over features inferred from generated structural alignment, 
as follows:

.. math::

   F = (N_{seq} + (N_{spe} * (-F_{energy})) + (N_{parts} * 10))

Where :math:`N_{seq}` is the final number of sequences, :math:`N_{spe}` is the number of species,
:math:`F_{energy}` corresponds to folding energy calculated using
``RNAalifold`` :cite:`Lorenz2011` and :math:`N_{parts}` accounts the number of
additional (:math:`> 1`) stem-loops on the reported consensus structure. The initial
population was :math:`A_{p}=40`, used operators were: *Selection* = Tournament,
:math:`n=39`; *Crossover* = Single point, probability=0.7; *Mutation* =
Displacement mutation, probability=0.1. The implementation were performed in
``Python`` v3.7.9 using ``deap`` package :cite:`Fortin:2012`. 

Finally, hidden Markov (HMMs) and covariance (CMs) models were build as
described in :cite:`Velandia:2021` using ``RNAalifold`` :cite:`Lorenz2011`
and ``Infernal`` package v.1.1.2 :cite:``.


.. rubric:: Footnotes

.. [#f3] Confidence of the annotation assigned by ``miRBase``, see `https://www.mirbase.org/blog/2014/03/high-confidence-micrornas/`
