package MiRNAture::HMM;

use Moose;
use MooseX::Types::Path::Class;
use Data::Dumper;

with 'MiRNAture::HMMsearch';

has 'genome_subject' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'subject_species' => (
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
	isa => 'ArrayRef[Str]',
	required => 1,
);

has 'path_covariance' => (
	is => 'ro',
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
	searchHomologyHMM($shift->genome_subject, $shift->subject_species, $shift->output_folder->stringify, $shift->path_hmm_models, $shift->path_covariance, $shift->bitscores_CM, $shift->length_CM, $shift->names_CM, $shift->families_names_CM, $shift->nhmmer_program_path->stringify, $shift->cmsearch_program_path->stringify, $zvalue, $minBitscore, $maxthreshold, $shift->list_models);
	my @result_files = check_folder_files($shift->output_folder->stringify."/".$shift->subject_species, "\.tab\.true\.table");
	for (my $i = 0; $i <= $#result_files; $i++) {
		my $hmm = $result_files[$i];
		$hmm =~ s/([A-Za-z]+\.)(.*)(\.tab\.true\.table)/$2/g;
		getSequencesFasta($shift->subject_species, $shift->genome_subject, $hmm, $shift->output_folder->stringify."/".$shift->subject_species, "2", $shift->length_CM, $shift->names_CM); #Header mode == 2 HMM
		searchStructureHMM($hmm, $shift->genome_subject, $shift->subject_species, $shift->output_folder->stringify, $shift->path_hmm_models, $shift->path_covariance, $shift->bitscores_CM, $shift->length_CM, $shift->names_CM, $shift->families_names_CM, $shift->nhmmer_program_path->stringify, $shift->cmsearch_program_path->stringify, $zvalue, $minBitscore, $maxthreshold, $shift->list_models);
	}
	return;
}	

sub clean_empty {
	my $shift = shift;
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_species); #Blast/Species
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_species."/Infernal"); #Blast/Species/Infernal
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_species."/Infernal/Final"); #Blast/Species/Infernal/Final
	return;
}

no Moose;
__PACKAGE__->meta->make_immutable;
