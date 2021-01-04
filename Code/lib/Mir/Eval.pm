package Mir::Eval;

use Moose;
use MooseX::Types::Path::Class;
use YAML::Tiny;
use Data::Dumper;
use Term::ANSIColor;
use lib "lib/miRNAture";

with 'miRNAture::ToolBox';

has 'all_parameters' => (
	is => 'ro',
	isa => 'Object',
	required => 1
);

sub create_folder_environment {
	my $shift = shift;
	my $config_file = shift;
	my $variables = shift;
	create_folders($shift->all_parameters->[3]->{"Default_folders"}->{"Output_folder"}, "Final_miRNA_evaluation");
	$variables->[4]->{"User_results"}{"Evaluation_results_folder"} = $shift->all_parameters->[3]->{"Default_folders"}->{"Output_folder"}."/Final_miRNA_evaluation";
	$variables->write($config_file);
	return;
}

sub test_input_files {
	my $shift = shift;
	my $config_file = shift;
	my $variables = shift;
    my $final_accepted_str = $shift->all_parameters->[3]->{"Default_folders"}->{"Output_folder"}."/miRNA_validation/accepted_".$shift->all_parameters->[3]->{Specie_data}->{Tag}.".miRNAs.txt";
    my $final_accepted_nostr = $shift->all_parameters->[3]->{"Default_folders"}->{"Output_folder"}."/miRNA_validation/accepted_".$shift->all_parameters->[3]->{Specie_data}->{Tag}."_noalign.txt";
    my $final_discarded = $shift->all_parameters->[3]->{"Default_folders"}->{"Output_folder"}."/miRNA_validation/discarded_".$shift->all_parameters->[3]->{Specie_data}->{Tag}.".miRNAs.txt";
    if (-e $final_accepted_str){
        #$variables->[4]->{"User_results"}{"Validated_miRNAs_with_structure"} = $final_accepted_str;
        $variables->[4]->{"User_results"}{"High_confidence_miRNAs"} = $final_accepted_str;
    } else {
        print_error("The accepted miRNA with structure file was not created");
    }
    if (-e $final_accepted_nostr){
        #$variables->[4]->{"User_results"}{"Validated_miRNAs_without_structure"} = $final_accepted_nostr;
        $variables->[4]->{"User_results"}{"Medium_confidence_miRNAs"} = $final_accepted_nostr;
    } else {
        print_error("The accepted miRNA without structure file was not created");
    }
    if (-e $final_discarded){
        #$variables->[4]->{"User_results"}{"Discarded_miRNAs"} = $final_discarded;
        $variables->[4]->{"User_results"}{"Low_confidence_miRNAs"} = $final_discarded;
    } else {
        print_error("The discarded structure file was not created");
    }
    $variables->write($config_file);
    return;
}

sub perform_evaluation {
    my $shift = shift;
    my $variables = shift;
    my $genome = $shift->all_parameters->[3]->{Specie_data}->{Old_Genome};
    my $specie_name = $shift->all_parameters->[3]->{Specie_data}->{Name};
    my $specie_tag = $shift->all_parameters->[3]->{Specie_data}->{Tag};
    #my @modes = ('validatedStr', 'validatedNoStr', 'discarded');
    my @modes = ('high_confidence', 'medium_confidence', 'low_confidence');
    my $input_table;
    my $out_fasta = $variables->[4]->{"User_results"}{"Evaluation_results_folder"};
    foreach my $md (@modes){
        if ($md eq 'low_confidence'){
            #$input_table = $variables->[4]->{"User_results"}{"Discarded_miRNAs"};
            $input_table = $variables->[4]->{"User_results"}{"Low_confidence_miRNAs"};
        } elsif ($md eq 'high_confidence'){
            #$input_table = $variables->[4]->{"User_results"}{"Validated_miRNAs_with_structure"};
            $input_table = $variables->[4]->{"User_results"}{"High_confidence_miRNAs"};
        } else {
            #$input_table = $variables->[4]->{"User_results"}{"Validated_miRNAs_without_structure"};
            $input_table = $variables->[4]->{"User_results"}{"Medium_confidence_miRNAs"};
        }
        my $fasta_file_out = getSequencesFasta_final($specie_name, $genome, $out_fasta, $input_table, $md);
        my $file_measure = perform_measure_MFE($fasta_file_out);
        if ($file_measure eq "NA"){
            print_result("The generation of MFE scores for $md mode is not available, there are no candidates.");
        } else {
            cross_table_final($file_measure, $input_table, "$out_fasta/${md}_${specie_tag}_final.table");
        }
    }
    return;
}

no Moose;
__PACKAGE__->meta->make_immutable;