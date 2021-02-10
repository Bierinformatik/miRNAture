package MiRNAture::ConfigFile;

use Moose;
use MooseX::Types::Path::Class;
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

sub read_genomes_paths { #File
	my $shift = shift;
	open (my $genomes_paths, "<", $shift->data_folder->stringify."/genomes.txt") or die "The file Data/genomes.txt> does not exists\n"; 
	my $targetGenomes;
	my %genomes;
	while (<$genomes_paths>){
		chomp;
		next if ($_ =~ /^#|^$/); 
		my @splitline = split /\=/, $_;
		$splitline[1] =~ s/"//g;
		$genomes{$splitline[0]} = $splitline[1];
	}
	foreach my $keys (sort keys %genomes){
		$targetGenomes .= "$keys ";
	}
	$targetGenomes =~ s/(.*)(\s+)/$1/g;
	$shift->genomes($targetGenomes);
	return %genomes;
}

no Moose;
__PACKAGE__->meta->make_immutable;
