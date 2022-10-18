![miRNAture](https://github.com/Bierinformatik/miRNAture/blob/main/mirnature_logo.png "miRNAture") 
=========
# Computational detection of microRNA candidates
[![License](https://img.shields.io/github/license/cavelandiah/miRNAture_v1)](https://github.com/cavelandiah/miRNAture_v1)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mirnature/README.html)
![Conda](https://img.shields.io/conda/v/bioconda/mirnature)
![Conda](https://img.shields.io/conda/dn/bioconda/mirnature)

## Description

Detection of miRNAs is a difficult problem. Due their small size limits the
available information and current sensitive methods, such as: `blast`, `nhmmer`,
or `cmsearch` are designed to increase sensitivity, but lead to an inevitable
large number of false positives only detected by detailed analysis of specific
features of typical miRNAs and/or conservation patterns in a structure-annotated
multiple sequence alignments.

The **miRNAture** pipeline implements a workflow specific to animal microRNAs
that automatizes homology search and validation steps.
On the homology search it combines two modes: sequence-homology by `blast` or/and 
`nhmmer` using query sequences or hidden markov models (HMMs), and structural 
validation performed by the `INFERNAL` package, using covariance models (CMs).
A merging step produces a final list of homology candidates. Over those
candidates a _Mature annotation_ stage performs a correction of the position of
mature sequences on the detected precursor and a structural evaluation 
in terms of minimum free energy (MFE), precursor length, folding and the
evaluation of anchored family specific-multiple secondary alignment 
(using `MIRfix`). Final sanity checks are performed on the _Evaluation_ stage, 
that reviews all the last mature annotation process, filtering the invalid candidates 
at structure level and reporting valid candidates on GFF3/BED and fasta files 
together with a summarize file that provides overall information about detected
miRNA candidates and families.

## Installation

The easiest way to install **miRNAture** is through `conda`. To do so, please first install
[conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

To speed up installation of dependencies and packages we suggest to use
[mamba](https://github.com/mamba-org/mamba), for this just run:

```
conda install mamba -c conda-forge
```

You can use `mamba` as drop-in replacement for `conda` by simply replacing the
call to `conda` with a call to `mamba`.


### Install via Conda

To install **miRNAture** from `conda` in a specific `mirnature` environment
simply run:

```
mamba create -n mirnature mirnature
```

if `mamba` is available, else run:

```
conda create -n mirnature mirnature
```

### Manual install, resolve dependencies via Conda

Create a `mirnature` `conda` environment with the file `miRNAture.yml`:

```
mamba env create -n mirnature -f miRNAture.yml
```

Activate the environment containing all dependencies:

```
conda activate mirnature
```

followed by the manual steps:

```
perl Build.PL
./Build
./Build test
./Build install
```

which will install **miRNAture** in the `mirnature` `conda` environment.


## Input files

The most important input file is a DNA sequence. This could be a multi-fasta 
sequence that belongs from a common specie (i.e. complete genome or group of 
particular sequences). At the same time, previous to execute **miRNAture** a
_pre-calculated_ dataset (that contains default data as CMs, HMMs, and required 
files to perform mature prediction) must be downloaded and correctly indicated
in the command line options with the flag `-dataF`. 

**New in version 1.1**
A new dataset containing all miRBase HMMs/CMs and validated mature sequences is
recommended to use as first approach to identify miRNAs over target species.
This dataset can be downloaded from **here**.

To run **miRNAture** in its _complete_ mode with default options, just run as:

```
# Activate the mirnature environment
conda activate mirnature

# Run miRNAture
./miRNAture -stage complete -dataF <Precalculated_folder> -speG <Target Genome> -speN <Specie_name> -speT <Tag_specie> -w <Output_dir> -m <Mode> (-str <Blast_strategy>) -blastq <Blast_queries_folder> 
```

## Output files

Final predicted miRNAs will be written on the `<Output_dir>` indicated with the `-w` flag.
The final candidates are described on the folder `Final_miRNA_evaluation/` as
follows:
```
Final_miRNA_evaluation/
├── Fasta/
├── MFE/
├── miRNA_annotation_Lach_accepted_conf.bed
├── miRNA_annotation_Lach_accepted_conf.gff3
├── miRNAture_summary_Lach.txt
└── Tables/
```

Inside this folder, **miRNAture** will create 3 folders containing their
correspondent results: sequences in `fasta` format (`Fasta/`), minimum free
energy and lengths from described sequences (`MFE/`) and the supporting
information ordered in tables for each annotated candidate (`Tables/`).
Additionally, associated genomic positions for the miRNA candidates are reported
in `BED` and `GFF3` formats and a summary file, `miRNAture_summary_*.txt`, that
describes overall descriptive statistics from found miRNA families. 

For detailed instructions how to use **miRNAture** please refer to the Manual pages:

* Through your favourite explorer, open [the manual pages here](http://www.bioinf.uni-leipzig.de/~cristian/miRNAture-Manual/).

## Pre-calculated datasets

Pre-calculated data composed by miRNA CMs, HMMs and required input files to
perform mature annotation has to be downloaded before run the full **miRNAture**
pipeline. Available datasets are listed below:

- **New dataset** containing metazoan curated miRBase v.22.1 families.
  Recommended for use with `miRNAture` v.1.1. Download [here](https://doi.org/10.5281/zenodo.7180160). 
- Required data to re-annotate human miRNAs: include CMs and HMMs build from
  miRBase without human sequences. Stored in Zenodo
  [here](https://zenodo.org/record/4531376#.YCQS8EMo_ys).
