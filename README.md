![miRNAture](https://github.com/cavelandiah/miRNAture_v1/blob/main/mirnature_logo.png "miRNAture") 
=========
## Computational detection of microRNA candidates
[![License](https://img.shields.io/github/license/cavelandiah/miRNAture_v1)](https://github.com/cavelandiah/miRNAture_v1) 
=========

### Description
Current approaches to computational miRNA detection relies on homology
relationships or detection of hairpin-loop candidates with lower folding energy.
A complete set of tools to automatizate this task have been assembled on
**miRNAture**. This current approach combines two different sequence-homology
modes, using `blast` or `HMMer`, and a secondary structure validation step,
performed by the `INFERNAL` package.  Combination of different
strategies and the modification of default covariance models, let **miRNAture**
report not only the homology relationships, but also define positions of mature
sequences and the correct miRNA annotation, supported by multiple family
specific alignments.

### Input files
The most important input file is a dna sequence. This could be a multifasta
sequence (i.e. complete genome or group of particular sequences) that belongs
from a common specie. Here I describe the general command line options to run
**miRNAture** in its _complete_ mode:

`./miRNAture -stage complete -speG <Target Genome> -speN <Specie_name> -speT <Tag_specie> -w <Output_dir> -mfx <MIRFix_path> -m <Mode> (-str <Blast_strategy>) -pe <SLURM_activation> -blastq <Blast_queries_folder>`

### Output files
A standard miRNA search will generate a detailed table of the final miRNA
annotation on the subject sequence(s). This result will be complemented by the
correspondent fasta sequences, not only from the reported hairpin, but from the
predicted mature, miR and miR\*, sequences. 


### General workflow
![workflow](https://github.com/cavelandiah/miRNAture_v1/blob/main/miRNAture2.png "miRNAture") 
