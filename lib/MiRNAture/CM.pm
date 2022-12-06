package MiRNAture::CM;

use Moose;
use MooseX::Types::Path::Class;
use Data::Dumper;

with 'MiRNAture::CMsearch';
with 'MiRNAture::ToolBox';
with 'MiRNAture::Evaluate';
with 'MiRNAture::Cleaner';

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
	required => 1,
	coerce => 1,
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

has 'cmsearch_program_path' => (
	is => 'ro',
	isa => 'Path::Class::File',	
	required => 1,
	coerce => 1,
);

sub create_folders_cm {
	my $shift = shift;
	create_folders($shift->output_folder->stringify, "");
	create_folders($shift->output_folder->stringify, $shift->subject_species); #Create Final Folder
	create_folders($shift->output_folder->stringify."/".$shift->subject_species, "Final"); #Create Final Folder
	return;
}

sub search_homology_CM {
	my ($shift, $zscore, $minBitscore,$maxthreshold) = @_;
	searchCMhomology($shift, $zscore,$minBitscore,$maxthreshold);
	my @result_files = check_folder_files($shift->output_folder->stringify."/".$shift->subject_species, $shift->subject_species."\.tab");
	for (my $i = 0; $i <= $#result_files; $i++) {
		my $cm_model = $result_files[$i];
		$cm_model =~ s/(RF[0-9]+)(.*)(\.tab)/$1/g;
		my $molecule = get_family_name($cm_model, $shift->families_names_CM);
		classify_2rd_align_results($shift->subject_species, $cm_model, $shift->output_folder."/".$shift->subject_species, $shift->output_folder."/".$shift->subject_species."/".$cm_model."\_".$shift->subject_species.".tab" ,"rfam", $molecule, $shift->bitscores_CM, $shift->length_CM, $shift->names_CM, $minBitscore, $maxthreshold);
	}
	return;
}	

sub clean_empty {
	my $shift = shift;
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_species); #Blast/Species
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_species."/Final"); #Blast/Species/Infernal
	return;
}

no Moose;
__PACKAGE__->meta->make_immutable;
