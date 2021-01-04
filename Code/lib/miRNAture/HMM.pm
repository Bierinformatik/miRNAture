package miRNAture::HMM;

use Moose;
use MooseX::Types::Path::Class;
use Data::Dumper;

with 'miRNAture::HMMsearch';

has 'hmm_model' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

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

has 'path_hmm_models' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	required => 1,
	coerce => 1,
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

with 'miRNAture::ToolBox'; #Use General tools role
with 'miRNAture::Cleaner';

sub create_folders_hmm {
	my ($shift, $work_folder) = @_;
	create_folders($work_folder, "HMMs")
}

sub search_homology_HMM {
	my ($shift, $zvalue, $minBitscore) = @_;
	searchHomologyHMM($shift->hmm_model, $shift->genome_subject, $shift->subject_specie, $shift->output_folder->stringify, $shift->path_hmm_models, $shift->path_covariance, $shift->bitscores_CM, $shift->length_CM, $shift->names_CM, $shift->families_names_CM, $shift->nhmmer_program_path->stringify, $shift->cmsearch_program_path->stringify, $zvalue, $minBitscore);
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
