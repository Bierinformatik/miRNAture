package Main::miRNAture;

use v5.26.2;
use strict;
use warnings;

our $VERSION = 1.00;

=head1 NAME

C<miRNAture> - Computational detection of microRNA candidates

=head1 SYNOPSIS

./miRNAture [-options]

=head1 OPTIONS

=over 12

=item -help           

Print this documentation

=item -man 

Prints the manual page and exits.

=item -blstq <PATH>

Path of blast query sequences in fasta format to be searched on the subject sequence.

=item -mfx <PATH>

Alternative path of the MIRfix.py program.

=item -m <Blast, HMM, Other_CM, Infernal, Final>

Homology search modes: Blast, HMM, Other_CM, Infernal and/or Final. It is possible to perform individual analysis, but it is always recommended to include the I<Final> option.

=item -pe <0|1>

Activate SLURM resource manager to submit parallel jobs (1) or not (0). 

=item -rep <default,200,100|relax,Number_Loci,Candidates_to_evaluate>

Setup number of maximum loci number that will be evaluated by the mature's annotation stage. It will detect the families that reported high number of loci, > 200 by default or a number X defined by the user. Then, it will select the top 100 candidates in terms of alignment scores to be evaluated in the next stage. The options are: 'default' and 'relaxed':
    -'default,200,100': Will label miRNA families with > 200 loci as repetititive. From them, the best 100 will be evaluated by the annotation stage, meanwhile the remaining will be reported as 'potential' candidates on the homology part.
    -'relax,Number_Loci,Candidates_to_evaluate': Allow the user to select the threshold to consider the families as repetitive, defined in the second number (Number_Loci). The number of candidates prone to be evaluated are defined by the third number (Candidates_to_evaluate).

=item -str <1,2,3,4,5,6,7,8,9,10>

This flag is blast mode specific. It corresponds to blast strategies that would be used to search miRNAs. It must be indicated along with -m I<Blast> flag. 


=item -stage <'homology','validation','evaluation','summarise','complete'>

Selects the stage to be run on miRNAture. The options are: 'homology', 'validation', 'evaluation', 'summarise' or 'complete'.

=item -speG <PATH>

Path of target sequences to be analyzed in fasta format.

=item -speN <Genera_specie>

Specie or sequence source's scientific name. The format must be: I<Genera_specie>, separated by '_'.

=item -speT <TAG_NAME>

Experiment tag. Will help to identify the generated files along miRNA output files.

=item -sublist <FILE_WITH_CM_NAMES>

Target list of CMs to be searched on subject genome/sequences. If not indicated, miRNAture will run all RFAM v14.4 metazoan miRNA models.

=item -w <OUT_PATH>

Working directory path to write all miRNAture results.

=back

=head1 DESCRIPTION

B<miRNAture> detects I<bona fide> miRNA candidates through sequence homology 
searches and validation steps using structural alignments with pre-defined or/and 
modified miRNA-specific covariance models. The miRNAture pipeline is composed of
three modules: (1) Homology search operating on miRNA precursors, (2) prediction 
of the positioning of mature miRNAs within the precursor mature annotation, and 
(3) an Evaluation scheme designed to identify false positive miRNA annotations.
This multi-stage approach generates annotation files in BED/GFF3 from precursors 
and detected mature regions and corresponding FASTA files. At the same time, a 
summary file with the MFE, precursor length and number of loci of each annotated 
miRNA family.

=head1 AUTHORS

I<Cristian A. Velandia Huerto>
I<Joerg Fallmann>
I<Peter F. Stadler>

=head1 LICENSE

GPL-3.0

=head1 BUGS, CAVEATS, COMPLAINS or DONATIONS 

Write directly to cristian at bioinf.uni-leipzig.de

=cut 

1;
