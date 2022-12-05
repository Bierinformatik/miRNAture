package Bio::miRNAture;

use v5.31.1;
use strict;
use warnings;

our $VERSION = v1.1;

=head1 NAME

C<miRNAture> - Computational detection of microRNA candidates

=head1 VERSION

miRNAture v.1.1 (Sept 21th, 2022)

=head1 SYNOPSIS

./miRNAture [-options]

=head1 AUTHORS:

I<Cristian A. Velandia Huerto>
I<Joerg Fallmann>
I<Peter F. Stadler>

=head1 OPTIONS

=over 12

=item -h/-help           

Print this documentation.

=item -man 

Prints an extended help page and exits.

=item -blstq/-blastQueriesFolder <PATH>

Path of blast query sequences in FASTA format to be searched on the subject sequence.

=item -dataF/-datadir <PATH>

Path to pre-calculated data directory containing RFAM and miRBase covariance, hidden markov models, 
and necessary files to run MIRfix. 

=item -debug_mode/dbug <0|1>

Activate the Perl DEBUGGER to run miRNAture.

=item -mfx/-mirfix_path <PATH>

Alternative path of the MIRfix.py program.

=item -m/-mode <blast,hmm,rfam,mirbase,user,final>

Homology search modes: blast, hmm, rfam, mirbase, infernal, user and/or final. It is possible to perform individual analysis, but it is always recommended to include the I<final> option to merge multiple results.

=item -ps/-parallel_slurm <0|1>

Activate SLURM resource manager to submit parallel jobs (1) or not (0). 

=item -pe/-parallel <0|1>

Activate parallel jobs processing (1) or not (0). 

=item -rep/-repetition_cutoff <relax,Number_Loci,Candidates_to_evaluate>

Setup number of maximum loci number that will be evaluated by the mature's annotation stage. By default, 
miRNAture will detect miRNA families that report high number of loci (> 200 loci). Then, it will select
the top 100 candidates in terms of alignment scores, as candidates for the validation stage (default,200,100). 
The designed values could be modified by the following flag: 'relax,Number_Loci,Candidates_to_evaluate'.
This option allows to the user to select the threshold values to detect repetitive families. The first 
parameter is <relax>, which tells miRNAture to change the default configuration. The next one, <Number_Loci> 
is the threshold of loci number to classify a family as repetitive. The last one, <Candidates_to_evaluate>, 
is the number of candidates prone to be evaluated in the next evaluation section. The rest candidates are
included as homology 'potential' candidates.

=item -str/-strategy <1,2,3,4,5,6,7,8,9,10>

This flag is blast mode specific. It corresponds to blast strategies that would be used to search miRNAs. It must be indicated along with -m I<Blast> flag. 

=item -stg/-stage <'homology','no_homology','validation','evaluation','summarise','complete'>

Selects the stage to be run on miRNAture. The options are: 'homology', 'no_homology', 'validation', 'evaluation', 'summarise' or 'complete'.

=item -speG/-species_genome <PATH>

Path of target sequences to be analyzed in FASTA format.

=item -speN/-species_name <Genera_species>

Species or sequence source's scientific name. The format must be: I<Genera_species>, separated by '_'.

=item -speT/-species_tag <TAG_NAME>

Experiment tag. Will help to identify the generated files along miRNA output files.

=item -sublist/-subset_models <FILE_WITH_CM_NAMES>

Target list of CMs to be searched on subject genome/sequences. If not indicated, miRNAture will run all RFAM v14.4 metazoan miRNA models.

=item -usrM/-user_models <PATH>

Directory with additional hidden Markov (HMMs) or covariance models (CMs) provided by the user. This must
be contain a corresponding HMMs/ and CMs/ folders, which the user models will be identified.

=item -w/-workdir <OUT_PATH>

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

=head1 BUGS, CAVEATS, COMPLAINS or DONATIONS 

Write directly to cristian at bioinf.uni-leipzig.de

=head1 LICENSE

GPL-3.0

=cut 

1;
