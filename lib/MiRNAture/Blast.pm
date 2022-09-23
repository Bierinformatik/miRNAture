package MiRNAture::Blast;
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(index_query_genome detect_blast_queries);

use Moose;
use MooseX::Types::Path::Class;
use Data::Dumper;
use MiRNAture::ToolBox;

has 'blast_str' => (
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

has 'query_folder' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	required => 1,
	coerce => 1,
);

has 'genome_subject' => (
	is => 'ro',
	isa => "Str",
	required => 1,
);

has 'subject_specie' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'makeblast_program_path' => (
	is => 'ro',
	isa => 'Path::Class::File',
	required => 1,
	coerce => 1,
);

with 'MiRNAture::ToolBox'; 
with 'MiRNAture::BlastPrepareQueries'; 

sub create_folders_blast {
	my $shift = shift;
	create_folders($shift->output_folder, "");
	return;
}

sub index_query_genome {
	#my $shift = shift;
	my ($genome, $makeblastpath) = @_;
	#Only create in case it doesn't exists
	if (!-e $genome.".nhr"){
		print_process("Generating BLAST DB for the target genome");
		make_blast_database($genome,$makeblastpath);
	} else {
		;
	}
	return;
}

sub resolve_mergings_blast {
	my $shift = shift;
	resolveBlastMergings();
	return;
}

no Moose;
__PACKAGE__->meta->make_immutable;
