#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Long qw(HelpMessage GetOptionsFromArray);
use Pod::Usage;
use Cwd qw(getcwd realpath);
use YAML::Tiny;
use FindBin qw($Bin);
use lib "$Bin/../lib";

use MiRNAture::Evaluate; #Evaluation CMsearch
use MiRNAture::ConfigFile; #Analyse configuration file class
use MiRNAture::Blast;
use MiRNAture::_BlastSearch;
use MiRNAture::BlastPrepareQueries;
use MiRNAture::HMMsearch; #HMM subroutine library
use MiRNAture::ToolBox;
use MiRNAture::HMM;
use MiRNAture::CM;
use MiRNAture::Others;
use MiRNAture::FinalCandidates;
use MiRNAture::LogFile;

### Input data
my $nameC = ""; #cmlist file
my @path_cm = ""; #cm_path
my @path_hmm = ""; #hmm_path
my $mode = ""; #run mode
my $specie = ""; #target specie genomes
my $strategy = ""; #Blast strategy
my $blastQueriesFolder = ""; # Path blast queries
my $work_folder = ""; #output folder
my $name_specie = ""; #Scientific name of specie
my $data_folder = "";
my $parallel_run = "";
my $rep_cutoff = ""; # Homology cutoff to speed-up the searches 

my $help = 0; 
my $man = 0;
my @cm_model_others = ""; #Possible path to new CMs != from RFAM ones
my @original_ARGV = @ARGV;
my @strategy;
GetOptions (
    'cmlist|l=s' => \$nameC,
    'cmpath|cmp=s{3}' => \@path_cm,
    'hmmpath|hmmp=s{3}' => \@path_hmm,
    'mode|m=s' => \$mode,
    'specie|spe=s' => \$specie,
    'specie_name|n_spe=s' => \$name_specie,
    'workdir|w=s' => \$work_folder,
    'repetition_cutoff|rep=s' => \$rep_cutoff,
    'help|h' => \$help,
    man => \$man,
    'strategy|str=s' => \@strategy,
    'blast_queries|blstq=s' => \$blastQueriesFolder,
    'data_folder|data=s' => \$data_folder,
    'parallel|pe=i' => \$parallel_run,
    'new_models|nmodels=s{2}' => \@cm_model_others,
) or pod2usage(2);
@strategy = split (/,/, join(',',@strategy));
#pod2usage(1) if $help;
#pod2usage(-exitstatus => 0, -verbose => 2) if $man;
#pod2usage("$0: No options given.") if ((scalar @ARGV) == 0);

my $startComplete = time();
my $current_dir = getcwd;
#### Flags Evaluation
evaluate_input_flags($nameC, $mode, $work_folder, $specie, $parallel_run,$rep_cutoff);
my $configuration_mirnature = read_config_file("$work_folder/../miRNAture_configuration_$specie.yaml");
## Working Paths
get_basic_files($configuration_mirnature->[3]->{Default_folders}->{Data_folder});

my ($start_hmm, $start_other, $start_infernal, $start_blast);
my %genomes;
my $name = $nameC;
$name =~ s/(\/Data\/|\.\/Data\/|\/.*\/|\.\/)(.*)/$2/g;
my $tag = ((strftime "%H%M%S%d%m%Y", localtime) + (int(rand(10)))); #Today date + random number 0..10.
print_process("Processing: $name\t$name_specie");

my $input_line = join " ", @original_ARGV; 

##Create Configuration and Log File
my $configuration_file = MiRNAture::ConfigFile->new(
    tag => $tag, 
    list_file => $nameC,
    data_folder => $data_folder,
);
$configuration_file->include_running_mode($mode);
%genomes = $configuration_file->read_genomes_paths();

my $log_file = MiRNAture::LogFile->new(
    log_file_input => "$work_folder/LOGs/miRNAture_homology_log_${specie}_${mode}_$tag.log",
    command_line => $input_line,		
);
$log_file->create_log_folder($work_folder);
$log_file->create_file;

#Key Output locations: Hard-Coded
my @species;
my $outHMM = "$work_folder/HMMs";
my $outBlast = "$work_folder/Blast";
my $outInfernal = "$work_folder/Infernal";
my $outOther = "$work_folder/Other";
my $out_final_path = "$work_folder";

##My Result files:
my ($blast_output, $hmm_output, $infernal_output, $other_output);

##Load basic databases
my $basicFiles = "$data_folder/Basic_files";
my ($bitscores, $len_r, $names_r, $names_r_inverse, $families_names);
my $user_data_path = $configuration_mirnature->[3]->{Default_folders}{"User_folder"};

# Infer user scores
if (length $user_data_path > 0){ #If path defined check that exists scores file
       	if (!-e "$user_data_path/all_user_scores.txt" || -z "$user_data_path/all_user_scores.txt"){ #Don't exist, infer
	       	print_process("Inferring scores from your miRNA CMs");
	       	infer_data_from_cm("$user_data_path/CMs", $user_data_path, "all_user_scores.txt"); #Here miRNAture generates the score's file: all_Other_scores.txt in the new CM folder.
      	}
}

if ($mode =~ m/OTHER_CM/){ #Load specific scores for user-provided CMs and mirbase.
    if (length $user_data_path > 0){ #If path defined check that exists scores file
	    ($bitscores, $len_r, $names_r, $names_r_inverse, $families_names) = load_all_databases("Additional", $basicFiles, $user_data_path);
    } else { #Only with miRNAture models
	    ($bitscores, $len_r, $names_r, $names_r_inverse, $families_names) = load_all_databases("mirbase", $basicFiles, $user_data_path); 
    }
} elsif ($mode =~ m/Infernal/) { #Here get all the scores from RFAM database
    ($bitscores, $len_r, $names_r, $names_r_inverse, $families_names) = load_all_databases("Basic", $basicFiles, $user_data_path); 
} else { #Here get all the scores from RFAM
	if (length $user_data_path > 0){ #If path defined check that exists scores file
		($bitscores, $len_r, $names_r, $names_r_inverse, $families_names) = load_all_databases("Joined", $basicFiles, $user_data_path); 
	} else {
		($bitscores, $len_r, $names_r, $names_r_inverse, $families_names) = load_all_databases("JoinedN", $basicFiles, $user_data_path); 
	}
}

## Start all
open my $LIST, "< $nameC" or die;
print_process("Processing $specie");

if (exists $genomes{$specie}){
    my $genome_path_complete = $genomes{$specie};
    my ($Zvalue, $genome_size_bp) = calculate_Z_value($genome_path_complete, "Genome");
    my $minBitscore = calculate_minimum_bitscore($genome_size_bp);
    #print "$Zvalue and $genome_size_bp and $minBitscore\n";
    if ($configuration_file->mode eq "BLAST"){
        print_process("Running Mode: ".$configuration_file->mode." searches");	
        write_line_log($log_file, "# Running Mode: ".$configuration_file->mode." at ".localtime."\n");
        $start_blast = time;
        detect_blast_queries($blastQueriesFolder);
        index_query_genome($genomes{$specie}, $configuration_mirnature->[2]->{Program_locations}->{makeblastdb});
        create_folders("$work_folder", "Blast");
        for (my $i = 0; $i<=$#strategy; $i++){
            print_process("Running on strategy $strategy[$i]");
            my $start_blast_str = time;
            my $blast_experiment = MiRNAture::_BlastSearch->new(
                blast_str => $strategy[$i],
                output_folder => "$work_folder/Blast",
                query_folder => $blastQueriesFolder,
                genome_subject => $configuration_mirnature->[3]->{Specie_data}->{Genome},
                subject_specie => $specie,
                data_folder => $data_folder,
                current_directory => $current_dir,
                path_covariance => \@path_cm,
                bitscores_CM => $bitscores,
                length_CM => $len_r,
                names_CM => $names_r,
                parallel_running => $parallel_run,
                blast_program_path => $configuration_mirnature->[2]->{Program_locations}->{blastn},
                cmsearch_program_path => $configuration_mirnature->[2]->{Program_locations}->{cmsearch},
                makeblast_program_path => $configuration_mirnature->[2]->{Program_locations}->{makeblastdb},
                models_list => $configuration_mirnature->[3]->{Default_folders}->{List_cm_miRNAs},
		user_data => $user_data_path,
            );	
            if ($blast_experiment->blast_str =~ /^\d+$/){
                my ($id_process_running, $molecules, $query_species, $families, $files_relation) = $blast_experiment->searchHomologySequenceBlast; #Run all blastn jobs, returns array by Str
                $blast_experiment->wait_processes($id_process_running, $parallel_run); #Wait until complete all processes from Str
                $blast_experiment->searchHomologyBlast($molecules, $query_species, $families, $files_relation, $Zvalue, $minBitscore);	
            } elsif ($blast_experiment->blast_str =~ /^ALL$/){
                $blast_experiment->join_all_blast_str;
                $blast_experiment->clean_empty;
            }
            my $diff = time - $start_blast_str;
            write_line_log($log_file, "# Running time for strategy ".$blast_experiment->blast_str." : ".$diff." s\n");
        }
        my $diff = time - $start_blast;
        write_line_log($log_file, "# Total time for ".$configuration_file->mode." : ".$diff." s\n");
    } elsif ($configuration_file->mode eq "HMM"){
        print_process("Running Mode: HMM searches");	
        write_line_log($log_file, "# Running Mode: ".$configuration_file->mode." at ".localtime."\n");
        $start_hmm = time;
        while (<$LIST>){ #In this case, list is the CM names
            chomp;
            #print_process("Running on $specie\t$_");
            my $hmm_experiment = MiRNAture::HMM->new(
                hmm_model => $_,	
                genome_subject => $genomes{$specie},
                subject_specie => $specie,
                output_folder => $outHMM,
                path_hmm_models => \@path_hmm,
                path_covariance => \@path_cm,
                bitscores_CM => $bitscores,
                length_CM => $len_r,
                names_CM => $names_r,
                families_names_CM => $families_names,
                nhmmer_program_path => $configuration_mirnature->[2]->{Program_locations}->{nhmmer},
                cmsearch_program_path => $configuration_mirnature->[2]->{Program_locations}->{cmsearch},
            );
            $hmm_experiment->create_folders_hmm($work_folder);
            $hmm_experiment->search_homology_HMM($Zvalue, $minBitscore);						
            $hmm_experiment->clean_empty;
        }
        #my $diff = $start - time;
        #LogFile::write_line_log("# Running homology search time: ".$diff." s\n");
    } elsif ($configuration_file->mode eq "INFERNAL"){
        print_process("Running Mode: ".$configuration_file->mode." searches");	
        write_line_log($log_file, "# Running Mode: ".$configuration_file->mode." at ".localtime."\n");
        $start_infernal = time;
        while (<$LIST>){ #In this case, list is the CM names
            chomp;
            next if $_ !~ /^RF/;
            #print_process("Running on $specie\t$_");
            my $cm_experiment = MiRNAture::CM->new(
                cm_model => $_,
                genome_subject => $genomes{$specie},
                subject_specie => $specie,
                output_folder => $outInfernal,
                path_covariance => \@path_cm,
                bitscores_CM => $bitscores,
                length_CM => $len_r,
                names_CM => $names_r,
                families_names_CM => $families_names,
                cmsearch_program_path => $configuration_mirnature->[2]->{Program_locations}->{cmsearch},
            );
            $cm_experiment->create_folders_cm;
            $cm_experiment->search_homology_CM($Zvalue,$minBitscore);
            $cm_experiment->clean_empty;
        }		
        #my $diff = $start - time;
        #LogFile::write_line_log("# Running homology search time: ".$diff." s\n");
    } elsif ($configuration_file->mode eq "OTHER_CM"){
        print_process("Running Mode: ".$configuration_file->mode." searches");	
        write_line_log($log_file, "# Running Mode: ".$configuration_file->mode." at ".localtime."\n");
        $start_other = time;
        while (<$LIST>){ #In this case, list is the CM names
            chomp;
            next if $_ =~ /^RF/;
            my $other_experiment = MiRNAture::Others->new(
                cm_model => $_,
                genome_subject => $genomes{$specie},
                subject_specie => $specie,
                output_folder => $outOther,
                path_covariance => \@cm_model_others, 
                bitscores_CM => $bitscores, 
                length_CM => $len_r,
                names_CM => $names_r,
                families_names_CM => $families_names,
                cmsearch_program_path => $configuration_mirnature->[2]->{Program_locations}->{cmsearch},
            );					
            $other_experiment->create_folders_other;
            $other_experiment->search_homology_other($Zvalue,$minBitscore);
            $other_experiment->clean_empty;
        }
        #my $diff = $start - time;
        #write_line_log("# Running search time: ".$diff." s\n");
    }
} else {
    print_error("Your sequence tag is not correct respect to genomes file");
}

print_process("Merging candidates on ".$configuration_file->mode);	
## Merging candidates by Strategy
if ($configuration_file->mode eq "HMM"){
    print_process("\t".$configuration_file->mode);	
    define_final_CMs("$outHMM/$specie/Infernal/Final", $specie, "HMM", $len_r);
    my $diff = time - $start_hmm;
    write_line_log($log_file, "# Total running time: ".$diff." s\n");
} elsif ($configuration_file->mode eq "INFERNAL"){
    print_process("\t".$configuration_file->mode);	
    define_final_CMs("$outInfernal/$specie/Final", $specie, "INFERNAL", $len_r);
    my $diff = time - $start_infernal;
    write_line_log($log_file, "# Total running time: ".$diff." s\n");
} elsif ($configuration_file->mode eq "OTHER_CM"){
    print_process("\t".$configuration_file->mode);	
    define_final_CMs("$outOther/$specie/Final", $specie, "INFERNAL", $len_r);
    my $diff = time - $start_other;
    write_line_log($log_file, "# Total running time: ".$diff." s\n");
}

if ($configuration_file->mode eq "Final"){
    print_process("Refining final candidates on ".$specie);	
    ####Preliminar Results ###
    $blast_output = "$outBlast/$specie/Infernal/Final/all_RFAM_${specie}_ALL.truetable.joined.table"; #Blast
    $hmm_output = "$outHMM/$specie/Infernal/Final/all_RFAM_${specie}.truetable.clean.joined.table"; #HMMs
    $infernal_output = "$outInfernal/$specie/Final/all_RFAM_${specie}.truetable.joined.table"; #Direct Infernal
    $other_output = "$outOther/$specie/Final/all_RFAM_${specie}.truetable.joined.table"; #Other direct CMs
    #$infernalDefault = "$outInfernalDefault/$specie/Final/all_RFAM_${specie}.truetable.joined.table"; #Direct Infernal
    ###
    if ($other_output){
        # Refill hashes with both score files
	if (length $user_data_path > 0){ #If path defined check that exists scores file
		($bitscores, $len_r, $names_r, $names_r_inverse, $families_names) = load_all_databases("Joined", $basicFiles, $user_data_path); 
	} else {
		($bitscores, $len_r, $names_r, $names_r_inverse, $families_names) = load_all_databases("JoinedN", $basicFiles, $user_data_path); 
	}
}
    my $final_candidates = MiRNAture::FinalCandidates->new(
        blast_results => $blast_output,
        hmm_results => $hmm_output,
        infernal_results => $infernal_output,
        other_results => $other_output,
        output_folder => $work_folder,	
        subject_specie => $specie,
        genome_subject => $configuration_mirnature->[3]->{Specie_data}{"Old_Genome"}, #$genomes{$specie}, # Here the original one
        names_CM => $names_r_inverse,
        length_CM => $len_r,
        families_names_CM => $families_names,
        specie_name => $name_specie,
        specie_genome_new_database => $configuration_mirnature->[3]->{Specie_data}{Database_names_genome},
        repetition_rules => $rep_cutoff,
    );
    $final_candidates->create_folders_final;
    $final_candidates->generate_final_output;
    $final_candidates->get_fasta_sequences;
    $final_candidates->generate_gff_homology;
    my $end = time(); 
    my $diff = $end - $startComplete;
    write_line_log($log_file, "# Total running time: ".$diff." s\n");
}

__END__

=head1 NAME
C<miRNAture> - Computational homology searches on a genomic DNA

=head1 SYNOPSIS

./run_homology.pl [-options]

Options:

-help	 brief help message

-man 	 full documentation

Author:

I<Cristian A. Velandia Huerto>

=head1 OPTIONS

=over 12

=item -help           

print this documentation

=item -man 

Prints the manual page and exits.

=item -cmlist            

list of CM to evaluate over target genome

=item -mode 

running mode: BLAST, HMM, INFERNAL, ALL, Final

=item -specie

target specie, write tag referenced on genomes file

=item -workdir 

current working directory

=item -strategy 

if selected BLAST as -mode write the strategy number  

=back 

=head1 DESCRIPTION

B<miRNAture> scan sequence fasta file for homologus metazoan microRNAs

=head1 AUTHOR

I<Cristian A. Velandia Huerto>

=head1 BUGS, CAVEATS, COMPLAINS or DONATIONS 
Write directly to cristian at bioinf.uni-leipzig.de

=cut 
