#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long 'HelpMessage';
use Pod::Usage;
use Cwd qw(getcwd);
use FindBin qw($Bin);
use lib "$Bin/../lib";

use MiRNAnchor::Main;
use MiRNAnchor::Tools;
use MiRNAnchor::Validate;
use MiRNAnchor::Check;
use MiRNAnchor::Classify;

### Input data
my $miRNA_fasta_path = ""; #Path of fasta sequences from candidates grouped by RFAM family
my $MIRFIX_path = ""; #Path of the exe file of MIRfix
my $working_path = ""; #Out folder
my $RFAM_mirbase_source = ""; #Path with all pre-calculated data from RFAM-miRBase
my $parallel_running = ""; #Decide if run with SLURM or not
my $origin_genome = ""; #Path of the genome from the subject specie
my $original_genome = ""; #Path of the original genome to perform blast mapping
my $final_miRNA_table = ""; #miRNA table database generate by miRNAture
my $tag_specie = "";
my $user_data = "";
my $help = 0; 
my $man = 0;

my $current_dir = getcwd;
my $address = "miRNA_validation";
my @candidates_collection;

GetOptions (
    'candidates|c=s' => \$miRNA_fasta_path,
    'mirfix|m=s' => \$MIRFIX_path,
    'outdir|o=s' => \$working_path,
    'extdata|e=s' => \$RFAM_mirbase_source,
    'parallel|p=s' => \$parallel_running,
    'genome|g=s' => \$origin_genome,
    'original_genome|og=s' => \$original_genome,
    'table_database|db=s' => \$final_miRNA_table,
    'tag_specie|tag=s' => \$tag_specie,
    'user_models|usrM=s' => \$user_data,
    'help|h' => \$help,
    man => \$man,
) or pod2usage(2);

my $start_config = MiRNAnchor::Main->new (
    current_dir => $current_dir,
    fasta_sequences => $miRNA_fasta_path,
    mirfix_path => $MIRFIX_path,
    output_folder => $working_path,
    precalculated_path => $RFAM_mirbase_source,
    parallel_running_confirmation => $parallel_running,
    subject_genome_path => $origin_genome,
    subject_original_genome_path => $original_genome,
    tag_spe_query => $tag_specie
);

my $db_models_relation = load_correspondence_models_rfam($current_dir, $user_data); #Load databases correspondence to do validation.
$start_config->check_existence_folder_output;
my $genomes_file = $start_config->get_genome_validation_list;
$start_config->recognize_families_homology($db_models_relation); #Include in final files the annotation family
my $RFAM_families = $start_config->identify_RFAM_families;  #Identify families already with the corresponding validation family
my $outTempAddress = ".mirfixTempIndividual";
my $genome_location = load_genomes_location_file($genomes_file); #Load database genomes

my %all = ();
my @all_results;

# Iterate Rfam Model
foreach my $rfam_defined_models (sort keys %{ $RFAM_families }){ #restrict only to the targeted ones
    my $eval = detect_multifamily_rfam($rfam_defined_models);
    my $rfam_precalculated_specific = $start_config->precalculated_path."/".$rfam_defined_models; #Data from working RFAM family.
    if (!-d $rfam_precalculated_specific){
        # Consider to put here an error, because all data were calculated
        warn "The corrected family $rfam_defined_models doesn't exists, not possible to validate the subject miRNAs\n";
        next;
    }
    #Detect flipped sequences
    #my $flipped = obtain_flip_sequences("$rfam_precalculated_specific/Output/$rfam_defined_models.out/$rfam_defined_models.json");
    my $candidates = process_query_fasta_specific($start_config->fasta_sequences->stringify, $rfam_defined_models, $start_config->output_folder, $start_config->tag_spe_query);
    my @candidates_group;
    # Iterate each sequence and run MIRfix
    foreach my $header_fasta (sort keys %{$candidates}){ #Add new genomes from query fasta
        my $specific_query_path;
        my $code_header = (split /\s+|\t/, $header_fasta)[1];
        if ($eval == 1){
            #Here detect the available groups
            my @groups = (0, 1); #Until now only 2 sets by corrected family.
            foreach my $num (@groups){
                my $subset = "${rfam_defined_models}_${num}"; #Set Name
                create_environment_detailed_subset($start_config->output_folder, $rfam_defined_models, $code_header, $address, $subset, $tag_specie);
                copy_rfam_files_group_subset($start_config->output_folder, $rfam_defined_models, $rfam_precalculated_specific, $code_header, $address, $subset, $tag_specie);
                create_genome_file_subset($start_config->output_folder, $rfam_defined_models, $rfam_precalculated_specific, $code_header, $address, $genomes_file, $subset, $tag_specie);
                process_query_fasta_group_subset($rfam_defined_models, $start_config->output_folder, $code_header, $header_fasta, $candidates, $address, $subset, $tag_specie); #Include subject fasta sequences
                $specific_query_path = $start_config->output_folder."/$address/$tag_specie/$rfam_defined_models/$subset/$code_header";
                include_subject_genome($specific_query_path, "BaseFiles", $subset, $genome_location, $start_config->subject_genome_path);
                clean_redundancy_file_subset($start_config->output_folder, $rfam_defined_models, $code_header, $address, $subset, $tag_specie);
                setup_all_to_run_mirfix_group_subset($specific_query_path, $subset, $start_config->mirfix_path, $current_dir, $start_config->parallel_running_confirmation, 
							$code_header, $outTempAddress, $start_config->tag_spe_query, $subset); #1: block copy again the mirfix input files 
                if ($start_config->parallel_running_confirmation == 1){ #Support of SLURM
                    run_mirfix_individual_subset($rfam_defined_models, $current_dir, $code_header, $outTempAddress, $start_config->tag_spe_query, $subset);
                } # Else is included on setup_all_to_run_mirfix_group
            }
        } else {
            create_environment_detailed($start_config->output_folder, $rfam_defined_models, $code_header, $address, $tag_specie);
            copy_rfam_files_group($start_config->output_folder, $rfam_defined_models, $rfam_precalculated_specific, $code_header, $address, $tag_specie);
            create_genome_file($start_config->output_folder, $rfam_defined_models, $rfam_precalculated_specific, $code_header, $address, $genomes_file, $tag_specie);
            process_query_fasta_group($rfam_defined_models, $start_config->output_folder, $code_header, $header_fasta, $candidates, $address, $tag_specie); #Include subject fasta sequences
            $specific_query_path = $start_config->output_folder."/$address/$tag_specie/$rfam_defined_models/$code_header";
            include_subject_genome($specific_query_path, "BaseFiles",$rfam_defined_models, $genome_location, $start_config->subject_genome_path);
            clean_redundancy_file($start_config->output_folder, $rfam_defined_models, $code_header, $address, $tag_specie);
            setup_all_to_run_mirfix_group($specific_query_path, $rfam_defined_models, $start_config->mirfix_path, $current_dir, $start_config->parallel_running_confirmation, $code_header, $outTempAddress, $start_config->tag_spe_query); #1: block copy again the mirfix input files 
            if ($start_config->parallel_running_confirmation == 1){ #Support of SLURM
                run_mirfix_individual($rfam_defined_models, $current_dir, $code_header, $outTempAddress, $start_config->tag_spe_query);
            } # Else is included on setup_all_to_run_mirfix_group
        }
        #Flipped sequences are now considered
        #remove_flip_sequences_group($start_config->output_folder, $rfam_defined_models, $rfam_precalculated_specific, $flipped, $code_header, $address);
    }
    if ($start_config->parallel_running_confirmation == 1){ #Support of SLURM
        wait_slurm_process($rfam_defined_models, $start_config->tag_spe_query);
    }
    
#Iterate again to collect all validated candidates after processing
    foreach my $header_fasta (sort keys %{$candidates}){ #Add new genomes from query fasta
        my $code_header = (split /\s+|\t/, $header_fasta)[1];
        if ($eval == 1){ # Multifamily 
            #Here detect the available groups
            my @groups = (0, 1); #Until now only 2 sets by corrected family.
            foreach my $num (@groups){
                my $subset = "${rfam_defined_models}_${num}"; #Set Name
                my $validated = validate_mirnas_group($start_config->output_folder."/".$address."/".$start_config->tag_spe_query, $rfam_defined_models, $candidates, $code_header, "Subset", $subset);
                if ($validated eq "NA"){
                    #open (my $DISCARDED, ">>", $start_config->output_folder."/".$address."/".$start_config->tag_spe_query."/".$rfam_defined_models."/".$rfam_defined_models.".discarded");
                    open (my $DISCARDED, ">>", $start_config->output_folder."/".$address."/".$start_config->tag_spe_query."/".$rfam_defined_models."/".$rfam_defined_models.".discarded");
                    my $id = (split /\s+|\t/, $header_fasta)[1];
                    print $DISCARDED "$id\t$header_fasta\tALIGNMENT_BROKEN\n";
                    close $DISCARDED;
                    next;
                } else {
                    push @all_results, $validated;
                }
            }
        } else {
            my $validated = validate_mirnas_group($start_config->output_folder."/".$address."/".$start_config->tag_spe_query, $rfam_defined_models, $candidates, $code_header, "Normal", "NA");
            if ($validated eq "NA"){
                open (my $DISCARDED, ">>", $start_config->output_folder."/".$address."/".$start_config->tag_spe_query."/".$rfam_defined_models."/".$rfam_defined_models.".discarded");
                my $id = (split /\s+|\t/, $header_fasta)[1];
                print $DISCARDED "$id\t$header_fasta\tALIGNMENT_BROKEN\n";
                close $DISCARDED;
                next;
            } else {
                push @all_results, $validated;
            }
        }
    }
    my %all;
    #Here build index of validated candidates inside evaluated family
    foreach my $final_registers (@all_results){
        my @split = split /\s+|\t/, $final_registers;
        my $info_final = join " ", @split[2..$#split];
        $all{$split[0]}{$split[-2]}{$split[1]}{$split[-1]} = $info_final;
    }
    # Print description candidates RF.out file
    print_results(\%all, $start_config->output_folder."/".$address."/".$start_config->tag_spe_query, $rfam_defined_models);
    @all_results = ();
    if ($eval == 1){
        my $data = $all{"Accepted"}{$rfam_defined_models};
        # Here select a valid 0 or 1 candidate from split if any accepted
        copy_best_element_subset($data, $start_config->output_folder."/".$address."/".$start_config->tag_spe_query, $rfam_defined_models, $current_dir);
    }
    %all = ();
    ##### Checkup 
    my $family_check = MiRNAnchor::Check->new(
        final_miRNA_database => $final_miRNA_table,
        source_miRNAs_fasta => $miRNA_fasta_path,
        output_folder => $working_path."/".$address."/".$start_config->tag_spe_query."/".$rfam_defined_models,
	    current_dir => $current_dir,
        accession_number_RFAM => $rfam_defined_models,
        subject_genome_path => $origin_genome,		
        subject_genome_path_original => $original_genome,
        subject_tag => $start_config->tag_spe_query,
        discarded_file => $working_path."/".$address."/".$start_config->tag_spe_query."/".$rfam_defined_models."/".$rfam_defined_models.".discarded",
        accepted_file => $working_path."/".$address."/".$start_config->tag_spe_query."/".$rfam_defined_models."/".$rfam_defined_models.".accepted",
        mature_description => $working_path."/".$address."/".$start_config->tag_spe_query."/".$rfam_defined_models."/".$rfam_defined_models.".out",
    );
    $family_check->process_fasta_sequences; #Start all required files and variables
    $family_check->check_candidate; #Based on produced files, run the steps of evaluation for each ACC number from CM models
    $family_check->concatenate_results("$working_path/$address", $start_config->tag_spe_query); #all_accepted/all_discarded
}

my $classification = MiRNAnchor::Classify->new(
    database_mirnas => $final_miRNA_table,
    output_folder => $working_path."/".$address."/".$start_config->tag_spe_query,
    current_dir => $current_dir,
    accepted_mirnas => $working_path."/".$address."/all_accepted_".$tag_specie."_miRNAs.txt", # Concatenated pre-accepted
    discarted_mirnas => $working_path."/".$address."/all_discarded_".$tag_specie."_miRNAs.txt", # Concatenated pre-discarded
    mature_description => $working_path."/".$address."/all_mature_".$tag_specie."_miRNAs_description.txt", # Concatenated mature info
    precalculated_path => $RFAM_mirbase_source,
    tag_spe_query => $tag_specie,
    accepted_file => $working_path."/".$address."/accepted_".$tag_specie.".miRNAs.txt", #Out accepted file
    accepted_noStr_file => $working_path."/".$address."/accepted_".$tag_specie."_noalign.txt", #Out accepted file with cand that destroyed alignment
    discarded_file => $working_path."/".$address."/discarded_".$tag_specie.".miRNAs.txt", #Out discarded file
    gff_high_file => $working_path."/".$address."/miRNA_annotation_".$tag_specie."_high_conf.gff3",
    bed_high_file => $working_path."/".$address."/miRNA_annotation_".$tag_specie."_high_conf.bed",
    gff_med_file => $working_path."/".$address."/miRNA_annotation_".$tag_specie."_medium_conf.gff3",
    bed_med_file => $working_path."/".$address."/miRNA_annotation_".$tag_specie."_medium_conf.bed",
    gff_NO_file => $working_path."/".$address."/miRNA_annotation_".$tag_specie."_NO_conf.gff3",
    bed_NO_file => $working_path."/".$address."/miRNA_annotation_".$tag_specie."_NO_conf.bed",
    gff_ACCEPTED_file => $working_path."/".$address."/miRNA_annotation_".$tag_specie."_accepted_conf.gff3",
    bed_ACCEPTED_file => $working_path."/".$address."/miRNA_annotation_".$tag_specie."_accepted_conf.bed",
);

$classification->process_all_candidates; #Generate evaluation at STO align level
$classification->generate_output_files; ## Finally create GFF3/BED files
$classification->organise_mess("$working_path/$address");

exit;
