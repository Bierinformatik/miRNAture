package MiRNAture::ConfigFile;

use Moose;
use MooseX::Types::Path::Class;
use MiRNAture::ToolBox;
use File::Path;
#use File::Copy;
use Data::Dumper;

has 'tag' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'list_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'data_folder' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1
);

has 'genomes' => (
	is => 'rw',
	isa =>  'Str',
	handles => { '_species' => 'species'}
);

has 'mode' => (
	is => 'rw',
	isa => 'Str',
	handles => { '_mode' => 'mode' }
);


sub include_running_mode {
	my $shift = shift;
	my $mode = shift;
	$shift->mode($mode);
	return;
}	

=head1 read_genomes_paths 
    Title: read_genomes_paths 
    Usage: read_genomes_paths(format_file, name_config_file);
    Function: Fills up the path information table along all
	  genomes. All of this information is retrieved from
	  Data/Basic_files/genomes.txt file.
    Returns: Hash with keys as Species tags, defined by the user, and values as paths. 
=cut 

sub read_genomes_paths {
    my ($shift, $tag, $path) = @_;
    my %genomes; 
    $genomes{$tag} = $path;
    return \%genomes;
}

=head1 subset_search_models 
    Title: subset_search_models
    Usage: subset_search_models(path_cms, path_hmms);
    Function: If defined a short list, select CMs to be searched on the 
	target genome.
    Returns: 0
=cut 

sub subset_search_models {
	my ($shift, $path_cm, $path_hmm, $path_user) = @_;
	my $new_cm = $shift->data_folder->stringify."/Selected_CM_models";
	my $new_hmm = $shift->data_folder->stringify."/Selected_HMM_models";
	# Create required folders to put temporal models
	if (!-e $new_cm){
		create_folders($shift->data_folder->stringify, "Selected_CM_models");
	} else {
		rmtree $new_cm;
		create_folders($shift->data_folder->stringify, "Selected_CM_models");
	}
	if (!-e $new_hmm){
		create_folders($shift->data_folder->stringify, "Selected_HMM_models");
	} else {
		rmtree $new_hmm;
		create_folders($shift->data_folder->stringify, "Selected_HMM_models");
	}
	my $numC = scalar @$path_cm;
	my $numH = scalar @$path_hmm;
	my $numU = scalar @$path_user;
	#Copy CMs
	open my $IN, "<", $shift->list_file->stringify;
	while(<$IN>){
		chomp;
		for (my $j = 0; $j <= $numC -1; $j++) {
			my $test_file = "$$path_cm[$j]/$_.cm";
			if (-e $test_file && !-l "$new_cm/$_.cm"){
				#copy($test_file, $new_cm);
				symlink($test_file, "$new_cm/$_.cm") || print_error("Cannot create CM reference");
			}
		}
		#Copy HMMs
		for (my $j = 0; $j <= $numH -1; $j++) {
			my $test_file = "$$path_hmm[$j]/$_.hmm";
			if (-e $test_file && !-l "$new_hmm/$_.hmm"){
				#copy($test_file, $new_hmm);
				symlink($test_file, "$new_hmm/$_.hmm") || print_error("Cannot create CM reference");
			}
		}
		#Copy CM and HMM user
		for (my $j = 0; $j <= $numU -1; $j++) {
			next if $_ =~ /^RF/;
			next if $_ =~ /^MIPF/;
			my $test_file = "$$path_user[$j]/HMMs/$_.hmm";
			if (-e $test_file && !-l "$new_hmm/$_.hmm"){
				#copy($test_file, $new_hmm);
				symlink($test_file, "$new_hmm/$_.hmm") || print_error("Cannot create CM reference");
			}
			my $test_file2 = "$$path_user[$j]/CMs/$_.cm";
			if (-e $test_file2 && !-l "$new_cm/$_.cm"){
				#copy($test_file2, $new_cm);
				symlink($test_file2, "$new_cm/$_.cm") || print_error("Cannot create CM reference");
			}
		}

	}
	my @new_cm;
	my @new_hmm;
	push @new_cm, $new_cm;
	push @new_hmm, $new_hmm;
	return (\@new_cm, \@new_hmm);
}

no Moose;
__PACKAGE__->meta->make_immutable;
