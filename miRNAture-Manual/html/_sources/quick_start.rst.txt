===========
Quick Start
===========

Input files
-----------

The most important input file is a DNA sequence. This could be a multifasta
sequence that belongs from a common specie (i.e. complete genome or group of
particular sequences). At the same time, previous to execute **miRNAture** a
you have to download a pre-calculated dataset (as indicated on :ref:`precalculated data` Section) 
that contains default data as CMs, HMMs, and required files to perform mature prediction. 
Once located in your computer, the path might be indicated with the flag ``-dataF``. 

Activate the ``mirnature`` environment
======================================
::

    conda activate mirnature

Run ``miRNAture``
=================
A complete mode should be run as follows::

    ./miRNAture -stage complete -dataF <Precalculated_folder> \
                -speG <Target Genome> -speN <Specie_name> \
                -speT <Tag_specie> -w <Output_dir> \
                -m <Mode> (-str <Blast_strategy>) \
                -blastq <Blast_queries_folder>

But it is always recommended to look up specific parameters. Do not use the
default parameters for all experiments.

Output files
============

Final predicted miRNAs will be written on the ``<Output_dir>`` indicated with the ``-w`` flag.
The final candidates are described on the folder ``Final_miRNA_evaluation/`` as
follows:

::

    Final_miRNA_evaluation/
    ├── Fasta/
    ├── MFE/
    ├── miRNA_annotation_<Tag_specie>_accepted_conf.bed
    ├── miRNA_annotation_<Tag_specie>_accepted_conf.gff3
    ├── miRNAture_summary_<Tag_specie>.txt
    └── Tables/

Inside this folder, **miRNAture** will create 3 folders containing their
correspondent results: sequences in ``fasta`` format (``Fasta/``), minimum free energy and
lengths from described sequences (``MFE/``) and the supporting information ordered in tables
for each annotated candidate (``Tables/``). Additionally, associated genomic positions
for the miRNA candidates are reported in ``BED`` and ``GFF3`` formats and a summary file,
``miRNAture_summary_<Tag_specie>.txt``, that describes overall descriptive statistics from found miRNA
families.

.. _precalculated data:

Pre-calculated datasets
=======================

Pre-calculated data composed by miRNA CMs, HMMs and required input files to perform mature annotation has
to be downloaded before run the full ``miRNAture`` pipeline. Available datasets are listed below:

 - **NEW**: Curated metazoan families from miRBase v.22.1, available to the structural validation stage.
 - Required data to re-annotate human miRNAs: include CMs and HMMs build from miRBase without human sequences. 
   Stored in Zenodo `here <https://zenodo.org/record/4531376#.YCQS8EMo_ys>`_.
