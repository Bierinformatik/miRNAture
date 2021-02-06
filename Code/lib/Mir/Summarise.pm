package Mir::Summarise;

use Moose;
use MooseX::Types::Path::Class;
use YAML::Tiny;
use Data::Dumper;
use Term::ANSIColor;
use File::Copy;
use lib "lib/MiRNAture";
use lib 'src/Statistics-R-0.34/lib';
use Statistics::R;

with 'MiRNAture::ToolBox';

has 'all_parameters' => (
	is => 'ro',
	isa => 'Object',
	required => 1
);

sub generate_summary_file {
	my $shift = shift;
	my $config_file = shift;
	my $variables = shift;
	my $validatedStr = $variables->[4]->{"User_results"}{"Evaluation_results_folder"}."/high_confidence_".$variables->[3]->{Specie_data}{Tag}."_final.table";
	my $validatedNoStr = $variables->[4]->{"User_results"}{"Evaluation_results_folder"}."/medium_confidence_".$variables->[3]->{Specie_data}{Tag}."_final.table";
	my $discarded = $variables->[4]->{"User_results"}{"Evaluation_results_folder"}."/NO_confidence_".$variables->[3]->{Specie_data}{Tag}."_final.table";
	my $out_summary_file = $variables->[4]->{"User_results"}{"Evaluation_results_folder"}."/miRNAture_summary_".$shift->all_parameters->[3]->{Specie_data}->{Tag}.".txt";
	calculate_summary_file($validatedStr, $validatedNoStr, $discarded, $out_summary_file);
	get_final_output($variables);
	return;
}

sub calculate_summary_file {
	my ($input_fileStr, $input_fileNoStr, $input_discarded, $outfile) = @_;
	open (my $OUT, ">", $outfile) or die "Not possible to generate summary file\n";
	print $OUT "###############\n# miRNAture v.1.0 Feb 1, 2021 \n# ".localtime()."\n###############\n#\n";
	my $R = Statistics::R->new();
	#Packages
	$R->run(q`suppressPackageStartupMessages(library(dplyr))`);
	#Input
    if (-e $input_fileStr && !-z $input_fileStr){
        $R->run(qq`tabStr = read.table(file='$input_fileStr', h=F)`);
    } else {
        $R->run(qq`tabStr = data.frame(matrix(,ncol=27, nrow=1))`);
        $R->run(qq`colnames(tabStr) = sub("X", "V", colnames(tabStr))`)
    }
    if (-e $input_fileNoStr && !-z $input_fileNoStr){
    	$R->run(qq`tabNoStr = read.table(file='$input_fileNoStr', h=F)`);
    } else {
        $R->run(qq`tabNoStr = data.frame(matrix(,ncol=27, nrow=1))`);
        $R->run(qq`colnames(tabNoStr) = sub("X", "V", colnames(tabNoStr))`)
    }
    if (-e $input_discarded && !-z $input_discarded){
        $R->run(qq`tabDiscarded = read.table(file='$input_discarded', h=F)`);
    } else {
        $R->run(qq`tabDiscarded = data.frame(matrix(,ncol=27, nrow=1))`);
        $R->run(qq`colnames(tabDiscarded) = sub("X", "V", colnames(tabDiscarded))`)
    }
	#Parsing and cleaning data
	#Chr, Start, End, Strand, Family, MFE, Len,Homology_Method
	$R->run(q`subsetDataStr = tabStr %>% select(c("V2","V7","V8","V3","V12","V4","V26","V27","V15"))`);
	$R->run(q`subsetDataNoStr = tabNoStr %>% select(c("V2","V7","V8","V3","V12","V4","V26","V27","V15"))`);
	$R->run(q`subsetDataDiscarded = tabDiscarded %>% select(c("V2","V7","V8","V3","V12","V4","V26","V27","V15"))`);
	#Column names
	$R->run(q`colnames(subsetDataStr) = c("Chr", "Start", "End", "Strand", "Family", "ACC", "MFE", "LEN", "Homology_Method")`);
	$R->run(q`colnames(subsetDataNoStr) = c("Chr", "Start", "End", "Strand", "Family", "ACC", "MFE", "LEN", "Homology_Method")`);
	$R->run(q`colnames(subsetDataDiscarded) = c("Chr", "Start", "End", "Strand", "Family", "ACC", "MFE", "LEN", "Homology_Method")`);
	#Categories
	$R->run(q`subsetDataStr$Class = "High"`);
	$R->run(q`subsetDataNoStr$Class = "Medium"`);
	$R->run(q`subsetDataDiscarded$Class = "NO"`);
	$R->run(q`subsetDataStr$ConservationStr = "SS_Conserved"`);
	$R->run(q`subsetDataNoStr$ConservationStr = "SS_modified"`);
	$R->run(q`subsetDataDiscarded$ConservationStr = "NA"`);
	#Merging ALL
	$R->run(q`subsetData = rbind(subsetDataStr,subsetDataNoStr,subsetDataDiscarded)`);
	$R->run(q`subsetData = subsetData %>% filter(!is.na(Chr))`);
    $R->run(q`familycountHigh = subsetData %>% subset(Class == "High") %>% group_by(Family, ACC, Class) %>% summarise(Total = n(), AverageMFE = mean(MFE), AverageLength = mean(LEN)) %>% arrange(Total)`); #Number loci by family
    $R->run(q`familycountMed = subsetData %>% subset(Class == "Medium") %>% group_by(Family, ACC, Class) %>% summarise(Total = n(), AverageMFE = mean(MFE), AverageLength = mean(LEN)) %>% arrange(Total)`); #Number loci by family

	$R->run(q`familycountHigh2 = as.data.frame(familycountHigh)`);
    $R->run(q`familycountMed2 = as.data.frame(familycountMed)`);

	$R->run(q`countfamiliesHigh = nrow(familycountHigh)`); #number of validated miRNA families detected
    $R->run(q`countfamiliesMed = nrow(familycountMed)`); #number of validated miRNA families detected
    
    #Count Accepted
	$R->run(q`total = subsetData %>% group_by(Class) %>% summarise(Total_miRNAs = n())`); #Count of total miRNAs
	$R->run(q`accepted = subsetData %>% subset(Class == "High" | Class == "Medium") %>% group_by(ConservationStr) %>% summarise(Total_miRNAs = n())`); #Count of total accepted miRNAs
    $R->run(q`acceptedHigh = subsetData %>% subset(Class == "High") %>% group_by(ConservationStr) %>% summarise(Total_miRNAs = n())`); #Count of total accepted miRNAs HIGH
    $R->run(q`acceptedMed = subsetData %>% subset(Class == "Medium") %>% group_by(ConservationStr) %>% summarise(Total_miRNAs = n())`); #Count of total accepted miRNAs MED

	$R->run(q`strand = subsetData %>% subset(Class == "High" | Class == "Medium") %>% group_by(Strand) %>% summarise(Total = n())`); #Count Strand
	$R->run(q`numberChr = subsetData %>% subset(Class == "High" | Class == "Medium") %>% group_by(Chr) %>% summarise(Total = n()) %>% summarise(Total_Genomic_Unit = sum(Total))`); 
	$R->run(q`methods = subsetData %>% subset(Class == "High" | Class == "Medium") %>% group_by(Homology_Method) %>% summarise(Total = n())`);
	my $total_miRNAs = $R->get('total');
	my $total_accepted = $R->get('accepted');
    my $total_accepted_High = $R->get('acceptedHigh');
    my $total_accepted_Med = $R->get('acceptedMed');


	my $number_familiesHigh = $R->get('countfamiliesHigh');
	my $table_detail_loci_High = $R->get('familycountHigh2');
    my $number_familiesMed = $R->get('countfamiliesMed');
	my $table_detail_loci_Med = $R->get('familycountMed2');

	my $table_strand = $R->get('strand');
	my $number_chr = $R->get('numberChr');
	my $methodCount = $R->get('methods');
	$R->stop();

	print $OUT "# Number of total detected miRNA loci\n"; 
	process_table_small($total_miRNAs, $OUT, "Num_loci");
	print $OUT "##\n";
	print $OUT "# Accepted miRNAs ##\n";
    process_table_small($total_accepted, $OUT, "Accepted_miRNAs");
    print $OUT "## High confidence miRNAs ##\n";
    print $OUT "## miRNA loci:\n";
    process_table_small($total_accepted_High, $OUT, "miRNA_loci_high");
	print $OUT "## miRNA families:\nNumber_families_high\t$number_familiesHigh\n";
	print $OUT "## Detail of miRNA loci:\n";
	print $OUT "## Section\tmiRNA_Family\tACC\tLoci\tAvg_MFE\tAvg_LEN\n";
	process_table($table_detail_loci_High, $OUT, "Detail_loci_high");
    print $OUT "## Medium confidence miRNAs ##\n";
    print $OUT "## miRNA loci:\n";
    process_table_small($total_accepted_Med, $OUT, "miRNA_loci_medium");
	print $OUT "## miRNA families:\nNumber_families_medium\t$number_familiesMed\n";
	print $OUT "## Detail of miRNA loci:\n";
	print $OUT "## Section\tmiRNA_Family\tACC\tLoci\tAvg_MFE\tAvg_LEN\n";
	process_table($table_detail_loci_Med, $OUT, "Detail_loci_medium");
	print $OUT "##\n";
	print $OUT "# Distribution of miRNA over the target genomic sequence:\n";
	print $OUT "# Genomic units with at least one miRNA: $$number_chr[-1]\n";
	print $OUT "# Section\tStrand\tLoci\n";
	process_table_small($table_strand, $OUT, "Strand");
	print $OUT "##\n";
	print $OUT "# miRNA detection methodology\n";
	print $OUT "# Section\tMethod\tLoci\n";
	process_table_small($methodCount, $OUT, "Detection_method");
	close $OUT;
	return;
}

### Subs
sub process_table {
	my $table = shift;
	my $OUT = shift;
    my $label = shift;
	my $number = scalar @$table;
    if ($number < 12){
        for (my $i=6; $i <= $number-1; $i += 6){
		    print $OUT "$label\t$$table[$i-5]\t$$table[$i-4]\t$$table[$i-2]\t$$table[$i-1]\t$$table[$i]\n";	
	    }
    } else {
        for (my $i=12; $i <= $number-1; $i += 7){
		    print $OUT "$label\t$$table[$i-5]\t$$table[$i-4]\t$$table[$i-2]\t$$table[$i-1]\t$$table[$i]\n";	
	    }
    }
	return;
}

sub process_table_small {
	my $table = shift;
	my $OUT = shift;
    my $label = shift;
	my $number = scalar @$table;
	for (my $i=12; $i <= $number-1; $i += 3){
		print $OUT "$label\t$$table[$i-1]\t$$table[$i-0]\n";	
	}
	return;
}

sub get_final_output {
	my $shift = shift;
	my $working_path = $shift->[4]->{"User_results"}{"Evaluation_results_folder"};
	create_folders($working_path, "Tables");
	create_folders($working_path, "Fasta");
	create_folders($working_path, "MFE");
	opendir DH, $working_path;
	while (my $file = readdir DH) {
		if ($file =~ /\.fasta$/){
			move("$working_path/$file", "$working_path/Fasta");
		} elsif ($file =~ /\.table/){
			move("$working_path/$file", "$working_path/Tables");
		} elsif ($file =~ /\.mfe/){
			move("$working_path/$file", "$working_path/MFE");
		} else {
			;
		}
	}
	my $tag = $shift->[3]->{Specie_data}{Tag};
	my $gffAccepted = $shift->[4]->{"User_results"}{"Output_miRNAnchor_folder"}."/GFF3/miRNA_annotation_".$tag."_accepted_conf.gff3";
	my $bedAccepted = $shift->[4]->{"User_results"}{"Output_miRNAnchor_folder"}."/BED/miRNA_annotation_".$tag."_accepted_conf.bed";
	copy($gffAccepted, $working_path);
	copy($bedAccepted, $working_path);
	return;
}

no Moose;
__PACKAGE__->meta->make_immutable;
