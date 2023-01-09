![miRNAture](https://github.com/Bierinformatik/miRNAture/blob/main/mirnature_logo.png "miRNAture") 
=========
# Computational detection of microRNA candidates
[![License](https://img.shields.io/github/license/cavelandiah/miRNAture_v1)](https://github.com/cavelandiah/miRNAture_v1)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mirnature/README.html)
![Conda](https://img.shields.io/conda/v/bioconda/mirnature)
![Conda](https://img.shields.io/conda/dn/bioconda/mirnature)

## Description

Detection of microRNAs is difficult, although miRNA precursors (about 80-100 nt) are often 
identified by sequence comparison due to their high conservation levels. Small mature products 
(about 22 nt) challenge sensitive homology-search methods, such as: `blast`, `nhmmer`, or `cmsearch`, 
gathering an inevitably large number of false positives. Recognizing true miRNAs with high 
confidence through fast and efficient computational approaches relies on detailed analysis of 
specific features from typical miRNAs and/or their conservation patterns in a structure-annotated 
multiple sequence alignment.

The **miRNAture** pipeline implements a workflow specific to animal miRNAs that automatizes homology 
search and miRNA-specific validation steps. The homology search combines two modes: sequence-homology by 
`blast` and/or `nhmmer` using query sequences or Hidden Markov Models (HMMs), and structural validation 
performed by the `INFERNAL` package, using pre-defined or user covariance models (CMs). 
In addition to standard homology-searches, **miRNAture** implements a merging step that produces a 
non-redundant list of homology candidates. Over those candidates a _Mature annotation_ stage (using `MIRfix`)
performs a correction of the position of mature sequences on the detected precursor and an exhaustive 
structural evaluation in terms of: minimum free energy (MFE), precursor length and folding, and the creation
of a *mature*-anchored family-specific multiple secondary alignment. Final sanity checks are performed 
during the _Evaluation_ stage, which reviews all the latest precursor and mature annotation steps, 
filtering out the invalid candidates at a structural level and reporting valid candidates on GFF3/BED and 
fasta files with summarized information about detected miRNA candidates and families.

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
sequence that belongs from a common species (i.e. complete genome or group of 
particular sequences). At the same time, previous to execute **miRNAture** a
_pre-calculated_ dataset (that contains default data as CMs, HMMs, and required 
files to perform mature prediction) must be downloaded and correctly indicated
in the command line options with the flag `-dataF`. 

**New in version 1.1**
A new dataset containing all miRBase v.22.1 HMMs/CMs and validated mature sequences is
recommended to use as first approach to identify miRNAs over target species.
This dataset can be downloaded from [here](https://doi.org/10.5281/zenodo.7180160).

To run **miRNAture** in its _complete_ mode with default options, just execute:

```
# Activate the mirnature environment
conda activate mirnature

# Run miRNAture
./miRNAture -stage complete -dataF <Precalculated_folder> -speG <Target Genome> -speN <Species_name> -speT <Species_tag> -w <Output_dir> -m <Mode> (-str <Blast_strategy>) -blastq <Blast_queries_folder> 
```

## Output files

Final predicted miRNAs will be written on the `<Output_dir>` indicated with the `-w` flag.
The final candidates are described on the folder `Final_miRNA_evaluation/` as
follows:
```
Final_miRNA_evaluation/
├── Fasta/
├── MFE/
├── miRNA_annotation_<Species_tag>_accepted_conf.bed
├── miRNA_annotation_<Species_tag>_accepted_conf.gff3
├── miRNAture_summary_<Species_tag>.txt
└── Tables/
```

Inside this folder, **miRNAture** will create 3 subfolders containing the
following results: sequences in `fasta` format (`Fasta/`), minimum free
energy and lengths from described sequences (`MFE/`) and the supporting information 
summarized in tabular format (`Tables/`).
Additionally, associated genomic positions for the miRNA candidates are reported
in `BED` and `GFF3` formats and a summary file, `miRNAture_summary_*.txt`, that
reports statistics from annotated miRNA families. 

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
