========
Strategy
========

Current approaches to computational miRNA detection relies on homology relationships 
or detection of hairpin-loop candidates with lower folding energy. A complete set of 
tools to automatize this task have been assembled on ``miRNAture``. This current 
approach combines two different sequence-homology modes, using ``blast`` or ``HMMer``, and 
a secondary structure validation step, performed by the ``INFERNAL`` package. Merging and 
consolidating task from multiple search strategies is done automatically by ``miRNAture``,
throwing at the the end of the *Homology searches* stage, a list of regions that
reported highest scores based on selected homology searches. Those candidates
passed designed filters (see in more detail in :numref:`filters`), to be considered as homologs by the
applied computational searches. 

Further structural and microRNA-specific evaluations are covered on the
*Mature evaluation* step, which makes use of an updated version of the original ``MIRfix`` 
pipeline :cite:`Yazbeck:19a`. At this step, ``miRNAture`` evaluates the identity at
family level of the homology candidates found in the previous step. Based on
that, makes use of the reported precursor, mature and genomic information contained on the 
``miRBase`` database. Specially for this step, this curation step is prepared
and reported with each release, allowing the user to perform the best mature positioning
and assignment on their predicted precursor sequences. Please refer to the
:ref:`appendix:Appendix` section to know more details about the curation process of 
the ``miRBase`` database and the generation of associated data. At this point, 
for each of those precursors, ``MIRfix`` will try to:

 - Assign the best-fitting mature sequence from those reported for the
   discovered miRNA family.
 - Predict the position of the miR* and correct the precursor sequence, based on
   the assigned mature.
 - Evaluate on a multiple structural alignment, the fit of the new annotated precursor in
   regard existing annotated miRNAs classified in the same family.

Those results are feeded onto additional evaluation steps, that would assign a
confidence to each precursor with its associated mature(s), namely: High, Medium
or No confidence.

.. _filters:

.. figure:: ./methodmirnature.png
   :width: 600
   :align: center
   :alt: miRNAture filters
   :figclass: align-center

   Designed homology/strucrure filters in ``miRNAture``. Specific programs used for each mode in parenthesis. Ann.: Annotation, SS: Secondary structure. CSS: Consensus secondary structure. ge: gathering cutoff from Rfam family. nBit = Bitscore/ge. ted: tree edit distance between default miRNA and modified multiple stockholm alignments. MFE: Minimum free energy. HSPs: high scoring pairs.

In the *complete* mode, ``miRNAture`` will report the following output files:

 - ``GFF3`` and ``BED`` files of the precursors with their mature sequences.
 - ``Fasta`` sequences from miRNA precurspors.
 - A summary table describing features of found miRNAs, such as: their *loci* number, family 
   classification and their confidence.

Current workflow is depicted in :numref:`workflow`:

.. _workflow:

.. figure:: ../../mirnature3.png
   :width: 600
   :alt: ``miRNAture`` workflow
   :align: center

   General ``miRNAture`` workflow.
