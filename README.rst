.. miRNAture documentation master file, created by
   Cristian Velandia on Tue Jul 7 23:21:30 2020.

miRNAture
=========

Computational detection of microRNA candidates.

Introduction
^^^^^^^^^^^^

``miRNAture`` is used to annotate homologs of metazoan microRNAs, discovered by means of
sequence, structural validation and correct position of their *mature* sequences,
detected on a harboring-hairpin-loop structure. This computational strategy enables 
an automated annotation workflow to search *bona fide* miRNA candidates, based
on pairwise comparisons and hidden Markov profiles, at sequence level. The
secondary structure is evaluated using covariance models and evaluating the
correct folding of the hairpin sequences, based on the position of the *mature*
sequences. 

Strategy
^^^^^^^^
Current approaches to computational miRNA detection relies on homology relationships 
or detection of hairpin-loop candidates with lower folding energy. A complete set of 
tools to automatize this task have been assembled on ``miRNAture``. This current 
approach combines two different sequence-homology modes, using ``blast`` or ``HMMer``, and 
a secondary structure validation step, performed by the ``INFERNAL`` package. Combination 
of different strategies and the modification of default covariance models, let ``miRNAture``
report not only the homology relationships, but also define positions of *mature* sequences 
and the correct miRNA annotation, supported by multiple family specific alignments.
Current workflow is depicted as follows:

.. image:: ../../mirnature2.png
   :width: 600

Installation
^^^^^^^^^^^^
Download the source code or clone the miRNAture project located on Github_.
On the ``Code/`` folder execute the file:: 

    ./activate_environment_conda.sh

Which will create the ``miRNAture`` conda repository (based on the
``miRNAture.yml`` file) and will install all the required programs to run
``miRNAture``.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   tutorial
   license
   help


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _Github: https://github.com/cavelandiah/miRNAture
