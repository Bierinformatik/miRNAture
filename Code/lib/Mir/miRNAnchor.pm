package Mir::miRNAnchor;

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

sub include_folder_environment {
    my $shift = shift;
    my $config_file = shift;
    my $variables = shift;
    create_folders($shift->all_parameters->[3]->{"Default_folders"}->{"Output_folder"}, "miRNA_validation");
    # all_RFAM_Latch_Final.ncRNAs_homology.txt.db
    my $database_file = $variables->[3]->{Default_folders}->{Output_folder}."/miRNA_prediction/Final_Candidates/all_RFAM_".$variables->[3]->{Specie_data}->{Tag}."_Final.ncRNAs_homology.txt.db";
    if (-e $database_file && !-z $database_file){
        $variables->[4]->{User_results}{Database_miRNA_file_result} = $database_file;
    } else {
        print_error("Seems that miRNAture did not detected candidates, makes no sense to run miRNAnchor");
    }
    $variables->[4]->{"User_results"}{"Output_miRNAnchor_folder"} = $shift->all_parameters->[3]->{"Default_folders"}->{"Output_folder"}."/miRNA_validation";
    $variables->write($config_file);
    return;
}

sub generate_file_mirnanchor {
    my $shift = shift;
    my $variable = shift;
    my $config_file = shift;
    my $program = $variable->[2]->{Program_locations}->{miRNAnchor};
    open (my $OUT, ">", $variable->[4]->{"User_results"}->{"Output_miRNAnchor_folder"}."/mirnanchor_run_".$variable->[3]->{Specie_data}->{Tag}.".sh");
    my $header = "# miRNAture\n# ".$variable->[0]->{"miRNAture"}->{"Author"}."\n# ".$variable->[0]->{"miRNAture"}->{"Date"}."\n# ".$variable->[0]->{"miRNAture"}->{"Version"};
    my $headerUser = "# Session data:\n# Hostname: ".$variable->[1]->{"Data_user"}->{"Hostname"}."\n# Date:".$variable->[1]->{"Data_user"}->{"Running_date"}."\n# User:".$variable->[1]->{"Data_user"}->{"User"}."\n# Program: miRNAnchor";
    print $OUT "#!/bin/bash\n\n####\n$header\n####\n$headerUser\n####\n";
    my $parameters = build_parameters_mirnanchor($variable);
    print $OUT "#Validation miRNAs\n";
    print $OUT $program." ".$parameters."\n";
    $variable->[4]->{"User_results"}{"miRNAnchor_program"} = $variable->[4]->{"User_results"}{"Output_miRNAnchor_folder"}."/mirnanchor_run_".$variable->[3]->{Specie_data}->{Tag}.".sh"; 
    $variable->write($config_file);
    close $OUT;
    return;
}

sub build_parameters_mirnanchor {
    my $variable = shift;
    # my $dbfile = "all_RFAM_".$variable->[3]->{Specie_data}->{Tag}."_Final.ncRNAs_homology.txt.db";
    my $parameters = "-c ".$variable->[3]->{Default_folders}->{Output_folder}."/miRNA_prediction/Final_Candidates/Fasta -m ".$variable->[2]->{"Program_locations"}{"MIRfix"}." -o ".$variable->[3]->{"Default_folders"}->{"Output_folder"}." -e ".$variable->[3]->{"Default_folders"}->{"Pre_calculated_validation_data"}." -og ".$variable->[3]->{"Specie_data"}->{"Old_Genome"}." -g ".$variable->[3]->{"Specie_data"}->{"Genome"}." -db ".$variable->[4]->{User_results}->{Database_miRNA_file_result}." -p ".$variable->[3]->{Homology_options}->{Parallel}." -tag ".$variable->[3]->{"Specie_data"}->{"Tag"};
    return $parameters;
}

sub run_miRNAnchor {
    my $shift = shift;
    my $variables = shift;
    my $run_file = $variables->[4]->{"User_results"}->{"miRNAnchor_program"};
    if (-e $run_file && !-z $run_file){
        print_result("I am ready to run miRNAnchor");
        system "chmod 755 $run_file";
        system "$run_file";
    } else {
        print_error("The file $run_file does not exists, fatal error");
    }
    return;
}

no Moose;
__PACKAGE__->meta->make_immutable;
