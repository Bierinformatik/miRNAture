package MiRNAture::FinalCandidates;

use Moose;
use MooseX::Types::Path::Class;
use Data::Dumper;

has 'blast_results' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
);
has 'hmm_results' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,

);
has 'infernal_results' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
);
has 'other_results' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
);
has 'user_results' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
);

has 'output_folder' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1,
);

has 'specie_name' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'subject_specie' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'genome_subject' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'specie_genome_new_database' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'length_CM' => (
	is => 'rw',
	isa => 'HashRef[Str]',
);

has 'names_CM' => (
	is => 'rw',
	isa => 'HashRef[Str]',
);

has 'families_names_CM' => (
	is => 'rw',
	isa => 'HashRef[Str]',
);

has 'final_out_table' => (
	is => 'rw',
	isa => 'Path::Class::File',
	coerce => 1,
	handles => { '_outtable' => 'table_file' }
);

has 'repetition_rules' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
    default => 'default,200,100',
);

with 'MiRNAture::ToolBox'; #Import set of subroutines
with 'MiRNAture::Merging';
with 'MiRNAture::Evaluate';
with 'MiRNAture::Cleaner';

sub create_folders_final {
	my $shift = shift;
	create_folders($shift->output_folder,"Final_Candidates");#Create Final Folder	
	create_folders($shift->output_folder."/Final_Candidates", "Fasta"); #Create Final fasta folder	
    # Genomes folder to save genome intervals to be searched on mirfix:  <30-08-21, cavelandiah> #
	create_folders($shift->output_folder."/Final_Candidates/Fasta", "Genomes"); #Create Final fasta folder	
	return;
}

sub generate_final_output {
	my $shift = shift;
	my @all_files = ($shift->blast_results->stringify, $shift->hmm_results->stringify, $shift->infernal_results->stringify, $shift->other_results->stringify, $shift->user_results->stringify);
	my @final_files;
	foreach my $possible_out_files (@all_files){
		if (-e $possible_out_files && !-z $possible_out_files){
			push @final_files, $possible_out_files;
		}
	}
	my $number_final_files = scalar @final_files;
	if ($number_final_files == 0){
		print_result("Were not detected blast, hmm, mirbase, rfam or user candidates. Not possible to merge anything");
		exit(0);
	} else {
		generate_final_ncRNAs($shift->blast_results->stringify, $shift->hmm_results->stringify, $shift->infernal_results->stringify, $shift->other_results->stringify, $shift->user_results->stringify, $shift->output_folder."/Final_Candidates", $shift->subject_specie);		
		format_final_table($shift->output_folder."/Final_Candidates", $shift->subject_specie, $shift->families_names_CM, $shift->names_CM, $shift->specie_genome_new_database);
        #my $final_file = $shift->output_folder->stringify."/Final_Candidates/all_RFAM_".$shift->subject_specie."_final.truetable.joined.table";
        my $temporal_file = $shift->output_folder->stringify."/Final_Candidates/all_RFAM_".$shift->subject_specie."_final.ncRNAs_homology.txt.temp";
        my $mode = (split /\,/, $shift->repetition_rules)[0]; 
        my $repeat_threshold = (split /\,/, $shift->repetition_rules)[1];
        my $number_best = (split /\,/, $shift->repetition_rules)[2];
        my $final_file = perform_detection_repeated_loci($temporal_file, $mode, $repeat_threshold, $number_best);
		$shift->final_out_table($final_file);
	}
	return;
}

sub create_final_report {
	my $shift = shift;
    return;
}

sub get_fasta_sequences {
	my $shift = shift;
	getSequencesFasta($shift->subject_specie, $shift->genome_subject, "NA", $shift->output_folder->stringify."/Final_Candidates/Fasta/", "5", $shift->length_CM, $shift->names_CM, $shift->final_out_table->stringify, $shift->specie_name); #Header mode == 5 Final table
	return;
}

sub get_small_genomes {
    # Generate genome anchors from predicted miRNAs to be validated with MIRfix.
    # This will extend reported coordinates +-300 nt.
	my $shift = shift;
    getSequencesFastaSubGenome($shift->specie_name, $shift->genome_subject, $shift->output_folder->stringify."/Final_Candidates/Fasta/Genomes", $shift->final_out_table);
	return;
}

sub generate_gff_homology {
    my $shift = shift;
    my $target_file = $shift->output_folder->stringify."/Final_Candidates/all_RFAM_".$shift->subject_specie."_final.ncRNAs_homology.txt.db"; #Complete_target file
    open my $OUT, "> ${target_file}.gff3" or die;
    open my $IN, "< $target_file";
    print $OUT "##gff-version 3\n";
    while (<$IN>){
        chomp;
        my @all = split /\s+|\t/, $_;
        print $OUT "$all[1]\tmiRNAture\tpre_miRNA\t$all[6]\t$all[7]\t$all[8]\t$all[2]\t\.\tID=$all[0];NAME=$all[-1];Alias=$all[3];Support=Homology:$all[14]\n"; 
    }
    close $OUT;
    close $IN;
    return;
}

no Moose;
__PACKAGE__->meta->make_immutable;
