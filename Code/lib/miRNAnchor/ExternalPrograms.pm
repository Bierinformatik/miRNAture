package miRNAnchor::ExternalPrograms;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(build_parameters run_mirfix clean_mirfix_logs time_measure total_time);

use Moose;
use MooseX::Types::Path::Class;
use Data::Dumper;

has 'cores' => (
	is => 'ro',
	isa => 'Int',
	required => 1,
);

has 'output_path' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1,
);

has 'families_path' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	required => 1,
	coerce => 1,
);

has 'list_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	required => 1,
	coerce => 1,
);

has 'genomes_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	required => 1,
	coerce => 1,
);

has 'mapping_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	required => 1,
	coerce => 1,
);

has 'mature_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	required => 1,
	coerce => 1,
);

has 'extension' => (
	is => 'ro',
	isa => 'Int',
);

has 'log_level' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

sub build_parameters {
	my $shift = shift;
	my $parameters = "-j ".$shift->cores." -o ".$shift->output_path->stringify."/ -i ".$shift->families_path->stringify."/ -f ".$shift->list_file->stringify." -g ".$shift->genomes_file->stringify." -m ".$shift->mapping_file->stringify." -a ".$shift->mature_file->stringify." -e ".$shift->extension." --loglevel ".$shift->log_level;
	return $parameters;
}

sub run_mirfix {
	my ($parameters, $path) = @_;
	my $starttime = time_measure();
	my $endCode = system "python3 $path $parameters";
	if ($endCode != 0){
        return 0;
    }
    my $final_time = total_time($starttime);
	return $final_time;
}

sub clean_mirfix_logs {
	my ($path_c, $out_path, $test) = @_;
	#$test = RF00053_MI0012744
	my $new_test = (split /\_/, $test)[0];
	my $log = "$path_c/log_$new_test";
	my $error = "$path_c/error_$new_test";
	my $log_move = "$out_path/log_$test";
	my $error_move = "$out_path/error_$test";
	if (-e $log){
		move($log, $log_move);
	}
	if (-e $error){
		move($error, $error_move);
	}
	return;
}

sub time_measure {
	my $start_program = time;
	return $start_program;
}

sub total_time {
	my $startt = shift;
	my $duration = time - $startt;
	return $duration;
}

no Moose;
__PACKAGE__->meta->make_immutable;
