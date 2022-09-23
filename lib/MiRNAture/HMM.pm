package MiRNAture::HMM;

use Moose;
use MooseX::Types::Path::Class;
use Data::Dumper;

with 'MiRNAture::HMMsearch';

#has 'hmm_model' => (
#	is => 'ro',
#	isa => 'Str',
#	required => 1,
#);

has 'genome_subject' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'subject_specie' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'output_folder' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1,
);

has 'list_models' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1,
);

has 'parallel_linux' => (
	is => 'ro',
	isa => 'Int',
	required => 1,
);

has 'path_hmm_models' => (
	is => 'ro',
	#isa => 'Path::Class::Dir',
	isa => 'ArrayRef[Str]',
	required => 1,
);

has 'path_covariance' => (
	is => 'ro',
	#isa => 'Path::Class::Dir',
	isa => 'ArrayRef[Str]',
	required => 1,
);

has 'bitscores_CM' => (
	is => 'rw',
	isa => 'HashRef[Str]',
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

has 'nhmmer_program_path' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1,
);

has 'cmsearch_program_path' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1,
);

with 'MiRNAture::ToolBox'; 
with 'MiRNAture::Cleaner';

sub create_folders_hmm {
	my ($shift, $work_folder) = @_;
	create_folders($work_folder, "HMMs")
}

sub search_homology_HMM {
	my ($shift, $zvalue, $minBitscore, $maxthreshold) = @_;
	searchHomologyHMM($shift->genome_subject, $shift->subject_specie, $shift->output_folder->stringify, $shift->path_hmm_models, $shift->path_covariance, $shift->bitscores_CM, $shift->length_CM, $shift->names_CM, $shift->families_names_CM, $shift->nhmmer_program_path->stringify, $shift->cmsearch_program_path->stringify, $zvalue, $minBitscore, $maxthreshold, $shift->list_models);
	##searchHomologyHMM($shift->hmm_model, $shift->genome_subject, $shift->subject_specie, $shift->output_folder->stringify, $shift->path_hmm_models, $shift->path_covariance, $shift->bitscores_CM, $shift->length_CM, $shift->names_CM, $shift->families_names_CM, $shift->nhmmer_program_path->stringify, $shift->cmsearch_program_path->stringify, $zvalue, $minBitscore, $maxthreshold);
	my @result_files = check_folder_files($shift->output_folder->stringify."/".$shift->subject_specie, "\.tab\.true\.table");
	for (my $i = 0; $i <= $#result_files; $i++) {
		my $hmm = $result_files[$i];
		$hmm =~ s/([A-Za-z]+\.)(.*)(\.tab\.true\.table)/$2/g;
		getSequencesFasta($shift->subject_specie, $shift->genome_subject, $hmm, $shift->output_folder->stringify."/".$shift->subject_specie, "2", $shift->length_CM, $shift->names_CM); #Header mode == 2 HMM
		searchStructureHMM($hmm, $shift->genome_subject, $shift->subject_specie, $shift->output_folder->stringify, $shift->path_hmm_models, $shift->path_covariance, $shift->bitscores_CM, $shift->length_CM, $shift->names_CM, $shift->families_names_CM, $shift->nhmmer_program_path->stringify, $shift->cmsearch_program_path->stringify, $zvalue, $minBitscore, $maxthreshold, $shift->list_models);
	}
	#create_folders($shift->output_folder."/".$shift->subject_specie."/Infernal","Final");#Create folder specific to specie
	#my @result_files_str = check_folder_files($shift->output_folder->stringify."/".$shift->subject_specie."/Infernal", "\.tab");
	#for (my $i = 0; $i <= $#result_files_str -1 ; $i++) {
	#	my $molecule = get_family_name($hmm, $families_names);
	#	classify_2rd_align_results($shift->subject_specie, $result_files_str[$i], $shift->output_folder->stringify."/".$specie."/Infernal", "$infernal_out_path/$specie.$hmm.tab" ,"hmm", $molecule, $bitscores, $len_r, $names_r, $minBitscore, $maxthreshold); #Obtain true candidates
	#}
	return;
}	

sub clean_empty {
	my $shift = shift;
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_specie); #Blast/Specie
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_specie."/Infernal"); #Blast/Specie/Infernal
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_specie."/Infernal/Final"); #Blast/Specie/Infernal/Final
	return;
}

no Moose;
__PACKAGE__->meta->make_immutable;
