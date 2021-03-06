package Mir::miRNAture;

use Moose;
use MooseX::Types::Path::Class;
use YAML::Tiny;
use Data::Dumper;
use Term::ANSIColor;
use lib "lib/MiRNAture";

with 'MiRNAture::ToolBox';
our $VERSION = v1.0.0;

has 'all_parameters' => (
	is => 'ro',
	isa => 'Object',
	required => 1
);

sub create_folder_environment {
	my $shift = shift;
	my $config_file = shift;
	my $variables = shift;
	create_folders($shift->all_parameters->[3]->{"Default_folders"}->{"Output_folder"}, "miRNA_prediction");
	$variables->[4]->{"User_results"}{"Output_miRNAture_folder"} = $shift->all_parameters->[3]->{"Default_folders"}->{"Output_folder"}."/miRNA_prediction";
	$variables->write($config_file);
	return;
}

sub test_name {
    my $name = shift;
    my $count = 1;
NAME:
    my $test_name = "$name-$count";
    if (-e "$test_name.sh"){
        $count++;
        goto NAME;        
    } else {
        return $test_name;
    }    
}

sub generate_file_mirnature {
	my $shift = shift;
	my $variable = shift;
	my $config_file = shift;
	my $program = $variable->[2]->{Program_locations}->{miRNAture};
    my $outfile = $variable->[4]->{"User_results"}->{"Output_miRNAture_folder"}."/mirnature_run_".$variable->[3]->{Specie_data}->{Tag}."_".$variable->[3]->{Homology_options}->{Mode};
    if (-e "$outfile.sh" && !-z "$outfile.sh" ){
        $outfile = test_name($outfile);
    }
    open (my $OUT, ">", "$outfile.sh");
	open (my $GENOME_OUT, ">>", $variable->[3]->{"Default_folders"}->{"Data_folder"}."/genomes.txt");
	my $header = "# miRNAture\n# ".$variable->[0]->{"miRNAture"}->{"Author"}."\n# ".$variable->[0]->{"miRNAture"}->{"Date"}."\n# ".$variable->[0]->{"miRNAture"}->{"Version"};
	my $headerUser = "# Session data:\n# Hostname: ".$variable->[1]->{"Data_user"}->{"Hostname"}."\n# Date:".$variable->[1]->{"Data_user"}->{"Running_date"}."\n# User:".$variable->[1]->{"Data_user"}->{"User"}."\n# Program: miRNAture";
	print $GENOME_OUT $variable->[3]->{"Specie_data"}->{"Tag"}."="."\"".$variable->[3]->{"Specie_data"}->{"Genome"}."\"\n";
	print $OUT "#!/bin/bash\n\n####\n$header\n####\n$headerUser\n####\n";
	my $parameters = build_parameters_mirnature($variable);
	my $mode = $variable->[3]->{"Homology_options"}->{"Mode"};
	if ($mode =~ m/Blast/i){
		print $OUT "#BLAST searches\n";
		print $OUT $program." ".$$parameters{blast}."\n";
	} 
	if ($mode =~ m/HMM/i){
		print $OUT "#HMM searches\n";
		print $OUT $program." ".$$parameters{hmm}."\n";
	}
	if ($mode =~ m/Infernal/i){
		print $OUT "#Infernal searches\n";
		print $OUT $program." ".$$parameters{infernal}."\n";
	}
	if ($mode =~ m/Other/i){
		print $OUT "#Other searches\n";
		print $OUT $program." ".$$parameters{others}."\n";
	}
	if ($mode =~ m/Final/i){
		print $OUT "#Final searches\n";
		print $OUT $program." ".$$parameters{final}."\n";
	}
	$variable->[4]->{"User_results"}{"miRNAture_program"} = $outfile.".sh";
	$variable->write($config_file);
	close $OUT;
	close $GENOME_OUT;
	return;
}

sub build_parameters_mirnature {
	my $variable = shift;
	my %parameters_homology;
	my $blast_strategies = join ",", @{$variable->[3]->{"Homology_options"}->{"Blast_strategies"}};
	my $parameters_blast = " -l ".$variable->[3]->{"Default_folders"}->{"List_cm_miRNAs"}." -m BLAST -str ".$blast_strategies." -blstq ".$variable->[3]->{"Default_folders"}->{"Blast_queries"}." -pe ".$variable->[3]->{"Homology_options"}->{"Parallel"}." -spe ".$variable->[3]->{"Specie_data"}->{"Tag"}." -n_spe ".$variable->[3]->{"Specie_data"}->{"Name"}." -w ".$variable->[4]->{"User_results"}->{"Output_miRNAture_folder"}." -data ".$variable->[3]->{"Default_folders"}->{"Data_folder"}. " -cmp ".$variable->[3]->{"Default_folders"}->{"CM_folder"}." ".$variable->[3]->{"Default_folders"}->{"Other_CM_folder"}." ".$variable->[3]->{"Default_folders"}->{"User_CM_folder"}." -hmmp ".$variable->[3]->{"Default_folders"}->{"HMM_folder"}." ".$variable->[3]->{Default_folders}{"Other_HMM_folder"}." ".$variable->[3]->{"Default_folders"}->{"User_HMM_folder"}." -rep ".$variable->[3]->{Homology_options}->{Repetition_threshold};
	my $parameters_hmm = "-l ".$variable->[3]->{"Default_folders"}->{"List_cm_miRNAs"}." -m HMM -pe 0 -spe ".$variable->[3]->{"Specie_data"}->{"Tag"}." -n_spe ".$variable->[3]->{"Specie_data"}->{"Name"}." -w ".$variable->[4]->{"User_results"}->{"Output_miRNAture_folder"}." -data ".$variable->[3]->{"Default_folders"}->{"Data_folder"}. " -cmp ".$variable->[3]->{"Default_folders"}->{"CM_folder"}." ".$variable->[3]->{"Default_folders"}->{"Other_CM_folder"}." ".$variable->[3]->{"Default_folders"}->{"User_CM_folder"}." -hmmp ".$variable->[3]->{"Default_folders"}->{"HMM_folder"}." ".$variable->[3]->{Default_folders}{"Other_HMM_folder"}." ".$variable->[3]->{"Default_folders"}->{"User_HMM_folder"}." -rep ".$variable->[3]->{Homology_options}->{Repetition_threshold};
	my $parameters_infernal = "-l ".$variable->[3]->{"Default_folders"}->{"List_cm_miRNAs"}." -m INFERNAL -pe 0 -spe ".$variable->[3]->{"Specie_data"}->{"Tag"}." -n_spe ".$variable->[3]->{"Specie_data"}->{"Name"}." -w ".$variable->[4]->{"User_results"}->{"Output_miRNAture_folder"}." -data ".$variable->[3]->{"Default_folders"}->{"Data_folder"}. " -cmp ".$variable->[3]->{"Default_folders"}->{"CM_folder"}." -hmmp ".$variable->[3]->{"Default_folders"}->{"HMM_folder"}." -rep ".$variable->[3]->{Homology_options}->{Repetition_threshold};
	my $parameters_others = "-l ".$variable->[3]->{"Default_folders"}->{"List_cm_miRNAs"}." -m OTHER_CM -pe 0 -nmodels ".$variable->[3]->{"Default_folders"}->{"Other_CM_folder"}." ".$variable->[3]->{"Default_folders"}->{"User_CM_folder"}." -spe ".$variable->[3]->{"Specie_data"}->{"Tag"}." -n_spe ".$variable->[3]->{"Specie_data"}->{"Name"}." -w ".$variable->[4]->{"User_results"}->{"Output_miRNAture_folder"}." -data ".$variable->[3]->{"Default_folders"}->{"Data_folder"}. " -cmp ".$variable->[3]->{"Default_folders"}->{"Other_CM_folder"}." ".$variable->[3]->{"Default_folders"}->{"User_CM_folder"}." -hmmp ".$variable->[3]->{Default_folders}{"Other_HMM_folder"}." ".$variable->[3]->{"Default_folders"}->{"User_HMM_folder"}." -rep ".$variable->[3]->{Homology_options}->{Repetition_threshold};
	my $parameters_final_homology = "-l ".$variable->[3]->{"Default_folders"}->{"List_cm_miRNAs"}." -m Final -pe 0 -spe ".$variable->[3]->{"Specie_data"}->{"Tag"}." -n_spe ".$variable->[3]->{"Specie_data"}->{"Name"}." -w ".$variable->[4]->{"User_results"}->{"Output_miRNAture_folder"}." -data ".$variable->[3]->{"Default_folders"}->{"Data_folder"}." -rep ".$variable->[3]->{Homology_options}->{Repetition_threshold};
	$parameters_homology{"blast"} = $parameters_blast;	
	$parameters_homology{"hmm"} = $parameters_hmm;
	$parameters_homology{"infernal"} = $parameters_infernal;
	$parameters_homology{"others"} = $parameters_others;
	$parameters_homology{"final"} = $parameters_final_homology;
	return \%parameters_homology;
}

sub run_miRNAture {
	my $shift = shift;
	my $variables = shift;
	my $run_file = $variables->[4]->{"User_results"}->{"miRNAture_program"};
	if (-e $run_file && !-z $run_file){
		print_result("I am ready to run miRNAture");
		system "chmod 755 $run_file";
		system "$run_file";
	} else {
		print_error("The file $run_file does not exists, fatal error");
	}
	return;
}

sub clean_cache {
    my $shift = shift;
    my $variables = shift;
    my $mode = shift;
    my $folder = $variables->[3]->{"Default_folders"}->{"Current_dir"}; 
    my $data_folder = $variables->[3]->{"Default_folders"}->{"Data_folder"}; 

    print_process("Cleaning temporal files...");
    if ($mode =~ /^All$/){
	    if (-e "$folder/.used_ids.txt"){
		    system("rm $folder/.used_ids.txt");
	    }
	    system("rm -r $folder/.blastTemp/");
	    system("rm -r $folder/.infernalTemp/");
	    if (-e "$folder/.mirfixTempIndividual/"){
		    system("rm -r $folder/.mirfixTempIndividual/");
	    }
	    system("rm $data_folder/genomes.txt");
	    system("rm $data_folder/Validation_mature_data/all_genomes_list.txt");
    } elsif ($mode =~ /^Homology/){
	    if (-e "$folder/.used_ids.txt"){
		    system("rm $folder/.used_ids.txt");
	    }
	    system("rm -r $folder/.blastTemp/");
	    system("rm -r $folder/.infernalTemp/");
    } elsif ($mode =~ /^Validation/){
	    if (-e "$folder/.mirfixTempIndividual/"){
		    system("rm -r $folder/.mirfixTempIndividual/");
	    }
	    system("rm $data_folder/genomes.txt");
	    system("rm $data_folder/Validation_mature_data/all_genomes_list.txt");
    } else {
	    ;
    }
    return;
}

sub move_results {
	my $shift = shift;
	my $variables = shift;
	my $folder = $variables->[3]->{"Default_folders"}->{"Current_dir"};
	my $output = $variables->[3]->{Default_folders}{"Output_folder"};
	_move_to_results($folder."/miRNAture_configuration_".$variables->[3]->{Specie_data}{"Tag"}.".yaml", $output);
}

sub _move_to_results {
	my ($file, $output) = @_;
	system("mv $file $output");
	return;
}

no Moose;
__PACKAGE__->meta->make_immutable;
