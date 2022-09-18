============
Introduction
============

MicroRNAs (miRNAs) have been characterized as important regulators 
in almost all animals and plants, as well as in unicellular eukaryotes :cite:`Moran:17`. 
Since their discovery in 1993 :cite:`Lee1993` and subsequently their recognition as 
a biological entity later from 2000 :cite:`Reinhart2000`, microRNAs were identified 
as key regulators on the temporal control of heterochronic genes on the nematode 
*Caenorhabditis elegans* as well as in another metazoan species, too :cite:`Pasquinelli:00`. 
This recognition was complemented by a complete characterization of their 
typical features as: a stem-loop structure, well-conservation over multiple metazoan clades
and typical expression patterns as isolated locus or co-expressed as polycistronic miRNA
transcripts :cite:`Lagos-Quintana:2001,Lau:2001,Lee:2001,Reinhart2000`.

Currently, metazoan miRNAs have been recognized as a conserved group of short ncRNAs, 
typically ~ 22 nt, with important roles as post-transcriptional regulators of the gene expression 
affecting a sizeable number of mRNAs and expressed in all developmental
process and diseases :cite:`Bartel:18`. Their canonical biogenesis starts in the form of primary
precursors (pri-miRNAs) transcribed from long non-coding RNAs or protein-coding transcripts,
mostly from introns :cite:`Zeidler:2020`. Later, derived from hairpin-like precursors excised in 
the nucleus (pre-miRNAs), their acting form is subsequently further processed as 
miRNA/miRNA* duplexes on the cytoplasm and incorporated into the RISC complex.
Target specificity is achieved by complementarity between the miR and mRNA sequence.
(see more details in :cite:`Bartel:18`).

Current classification of annotated miRNAs into families are available in miRBase [#f1]_ and
mirGeneDB [#f2]_ databases. As an example, the human genome reported 1984 miRNA precursors in 
the miRBase v.22.1 :cite:`Kozomara:19` and the corresponding mature products were estimated ~2300
:cite:`Alles:19`. Focusing on the *confidently canonical* miRNAs reported in :cite:`Fromm:15`,
the number of miRNAs is 519. Those differences are explained on the basis of multiple 
*miRNA* detection methods as well as the intrinsic definitions to define a *canonical miRNA*.

Despite the small size of the precursors (80-100 nt), sequence comparison
methods are able to detected them, due a high level of sequence conservation
:cite:`Price:11`. In one hand, it is important to point that the use of blast-based homology 
searches alone tend to produce false positives that require extensive curation, which relies
on properties obtained from miRNAs :cite:`Tarver:18`. On the other hand, the
inclusion of the consensus structure complements the homology search methods,
for example using covariance models (CMs) :cite:`Eddy:94, Gardner:2009`. The accuracy and 
sensitivity of CMs depends critically on the sequence alignment and the annotated consensus 
used to build the model. Those observations call for an integrated workflow to perform 
homology search and to evaluate their results consistently.

In this computational approach, focused on the *canonical* miRNAs processed by Drosha and Dicer,
we improve on ideas from MIRfix :cite:`Yazbeck:19a` and integrate it with homology search.
``miRNAture`` is used to identify and annotate homologs of metazoan microRNAs in a homology-based
setting.

.. rubric:: Footnotes

.. [#f1] https://www.mirbase.org/
.. [#f2] https://www.mirgenedb.org/
