package miRNAture::Others;

use Moose;
use MooseX::Types::Path::Class;
use Data::Dumper;

with 'miRNAture::Othersearch';
with 'miRNAture::ToolBox';
with 'miRNAture::Evaluate';
with 'miRNAture::Cleaner';

has 'cm_model' => (
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
	required => 1,
	coerce => 1,
);

has 'path_covariance' => (
	is => 'ro',
	isa => 'Path::Class::Dir',	
	required => 1,
        coerce => 1,
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

has 'cmsearch_program_path' => (
	is => 'ro',
	isa => 'Path::Class::File',	
	required => 1,
        coerce => 1,
);

sub create_folders_other {
	my $shift = shift;
	create_folders($shift->output_folder->stringify, "");
	create_folders($shift->output_folder->stringify, $shift->subject_specie); #Create Final Folder
	create_folders($shift->output_folder->stringify."/".$shift->subject_specie, "Final"); #Create Final Folder
	return;
}

sub search_homology_other {
	my $shift = shift;
    my $zscore = shift;
    my $minBitscore = shift;
	searchOthershomology($shift, $zscore, $minBitscore);
	return;
}	

sub clean_empty {
	my $shift = shift;
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_specie); #Blast/Specie
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_specie."/Final"); #Blast/Specie/Infernal
	return;
}

no Moose;
__PACKAGE__->meta->make_immutable;