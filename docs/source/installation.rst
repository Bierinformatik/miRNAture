============
Installation
============

The easiest way to install ``miRNAture`` is through ``conda``. To do so, please first install
Conda [#conda]_.

To speed up installation of dependencies and packages we suggest to use mamba [#mamba]_, for 
this just run::
    
    conda install mamba -c conda-forge

You can use ``mamba`` as drop-in replacement for ``conda`` by simply replacing the call to ``conda`` 
with a call to ``mamba``.

*Install via ``conda``*

To install ``miRNAture`` from ``conda`` in a specific ``mirnature`` environment simply run::

    mamba create -n mirnature mirnature

if ``mamba`` is available, else run::

    conda create -n mirnature mirnature

*Manual install, resolve dependencies via ``conda``*

Create a `mirnature` ``conda`` environment with the file ``miRNAture.yml``::

    mamba env create -n mirnature -f miRNAture.yml

Activate the environment containing all dependencies::

    conda activate mirnature

followed by the manual steps::

    perl Build.PL
    ./Build
    ./Build test
    ./Build install

which will install ``miRNAture`` in the `mirnature` ``conda`` environment.

.. rubric:: Footnotes

.. [#conda] https://docs.conda.io/projects/conda/en/latest/user-guide/install/
.. [#mamba] https://github.com/mamba-org/mamba
.. _Github: https://github.com/Bierinformatik/miRNAture
