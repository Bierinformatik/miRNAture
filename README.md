![miRNAture](https://github.com/Bierinformatik/miRNAture/blob/main/mirnature_logo.png "miRNAture") 
=========
## Computational detection of microRNA candidates
[![License](https://img.shields.io/github/license/cavelandiah/miRNAture_v1)](https://github.com/cavelandiah/miRNAture_v1) 
=========

### Description

Detection of miRNAs is a difficult problem due their small size limits the
available information. Current sensitive methods as: parameter optimized
`blast`, `nhmmer`, or `cmsearch` runs designed to increase sensitivity, but
leading to an inevitable large number of false positives only detected by
detailed analysis of specific features of typical miRNAs and/or analysis of
conservation patterns in a structure-annotated multiple sequence alignments.

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

### Installation

The easiest way to install **miRNAture** is through `conda`. To do so, please first install
[conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

To speed up installation of dependencies and packages we suggest to use
[mamba](https://github.com/mamba-org/mamba), for this just run:

```
conda install mamba -c conda-forge
```

You can use `mamba` as drop-in replacement for `conda` by simply replacing the call to `conda` with a call to
`mamba`.


#### Install via Conda

To install **miRNAture** from `conda` in a specific `mirnature` environment simply run:

```
mamba create -n mirnature mirnature
```

if `mamba` is available, else run:

```
conda create -n mirnature mirnature
```

#### Manual install, resolve dependencies via Conda

Create a `mirnature` `conda` environment with the file `miRNAture.yml`:

```
mamba env create -n mirnature -f Code/miRNAture.yml`
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


### Input files

The most important input file is a DNA sequence. This could be a multifasta
sequence (i.e. complete genome or group of particular sequences) that belongs
from a common specie. Here I describe the general command line options to run
**miRNAture** in its _complete_ mode:

```
# Activate the mirnature environment
conda activate mirnature

# Run miRNAture
./miRNAture -stage complete -speG <Target Genome> -speN <Specie_name> -speT <Tag_specie> -w <Output_dir> -m <Mode> (-str <Blast_strategy>) -blastq <Blast_queries_folder> -rep default,150,100
```

### Output files

Final predicted miRNAs will be written on the `<Output_dir>` indicated with the `-w` flag.
The final candidates are described on the folder `Final_miRNA_evaliation/` as
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

Inside the folder `Final_miRNA_evaluation/` will be created 3 folders containing their
correspondent results: sequences in `fasta` format (`Fasta/`), minimum free energy and 
lengths from described sequences (`MFE/`) and the supporting information ordered in tables
for each annotated candidate (`Tables/`). Additionally, associated genomic positions 
for the miRNA candidates are reported in `BED` and `GFF3` formats and a summary file 
describes overall descriptive statistics from found miRNA families. 

For detailed instructions how to use **miRNAture** please refer to the Manual pages:

* In this ![PDF](miRNAture-Manual/latex/mirnature.pdf) version.
* Through your favorite explorer, open the `index.html` file located at `miRNAture-Manual/html/index.html`.

### Precalculated datasets

Pre-calculated data composed by miRNA CMs, HMMs and required input files to perform mature annotation has
to be downloaded for the user to run the full pipeline. Available datasets are listed as follow:

- Required data to re-annotate human miRNAs: include CMs and HMMs build from miRBase without human 
sequences. Stored in Zenodo [here](https://zenodo.org/record/4531376#.YCQS8EMo_ys).
	
