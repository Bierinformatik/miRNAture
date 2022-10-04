package MiRNAture::LogFile;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(write_line_log create_file _create_new_file _header_log create_header create_log_folder);
use Moose;
use MooseX::Types::Path::Class;
use Data::Dumper;
use POSIX 'strftime';

has 'log_file_input' => (
	is => 'ro',
	isa => 'Path::Class::File',
	required => 1,
	coerce => 1
);

has 'command_line' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

sub write_line_log {
	my $shift = shift;
	my $line = shift;
	if (-e $shift->log_file_input->stringify){
		open(my $OUT, ">>", $shift->log_file_input->stringify); 
		print $OUT "$line";
		close $OUT;
	}
	return;
}

sub create_log_folder {
	my $shift = shift;
	my $path = shift;
	my $directory = $path."/LOGs";
	if (!-e $directory) {
		system "mkdir $directory";	
	} else {
		;
	}
	return;
}

sub create_file {
	my $shift = shift;
	_create_new_file($shift);
	_header_log($shift);
	return;
}

sub _create_new_file {
	my $shift = shift;
	open(my $OUT, ">", $shift->log_file_input->stringify) or die;
	close $OUT;
	return;
}
sub _header_log {
	my $shift = shift;
	my $header = create_header();
	write_line_log($shift, $header);
	write_line_log($shift, "# User:".`whoami`);
	write_line_log($shift, "# Host:".`hostname`);
	write_line_log($shift, "# Start Time: ".localtime."\n");
	write_line_log($shift, "# Executed Command line: ./miRNAture.pl ".$shift->command_line."\n");
	write_line_log($shift, "#\n");
	return;
}

sub create_header {
	my $header = "######\n###  miRNAture : search miRNAs against a sequence\n###  Authors: Cristian A. Velandia-Huerto, Jöerg Fallmann and Peter F. Stadler\n###  Version: v.1.10 Sept 21th, 2022\n###  Bioinformatics Department, Universität Leipzig\n###  Freely distributed under a BSD open source license.\n######\n";
	return $header;
}

no Moose;
__PACKAGE__->meta->make_immutable;
