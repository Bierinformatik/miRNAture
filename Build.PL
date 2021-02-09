use 5.26.2;
use strict;
use warnings;
use Module::Build;
Module::Build->VERSION('0.4004');

my $builder = Module::Build->new(
	module_name         => 'Bio::miRNAture',
	license             => 'gpl_3',
	dist_author         => q{Cristian A Velandia Huerto <cristian@bioinf.uni-leipzig.de>},
	dist_version_from   => 'lib/Bio/miRNAture.pm',
	release_status      => 'stable',
	configure_requires => {
		'Module::Build' => '0.4224',
	},
	test_requires => {
		'Test::More' => '0',
	},
	requires => {
		'perl' => 'v5.26.2',
		'Getopt::Long' => 2.5, 
		'Pod::Usage' => 1.69,
		'Cwd' =>  3.75,
		'FindBin' => 1.51,
		'File::Share' => 0.25,
		'Exporter' => 5.72,
		'File::Basename' => 2.85,
		'File::Copy' => 2.32,
		'File::Copy::Recursive' => 0.45,
		'JSON' => 4.02,
		'List::MoreUtils' => 0.428,
		'List::Util' => 1.52,
		'Math::BigFloat' => 1.999806,
		'Moose' => 2.2011,
		'MooseX::Types::Path::Class' => 0.09,
		'POSIX' => 1.76, 
		'Scalar::Util' => 1.52,
		'Statistics::R' => 0.34,
		'Term::ANSIColor' => 4.06,
		'YAML::Tiny' =>  1.73,
	},
	script_files => [qw(script/miRNAture script/miRNAture.pl script/miRNAnchor.pl script/evaluate_conserved_str.py)],
	add_to_cleanup => [ 'Bio-miRNAture-*' ],
	share_dir => 'share',
);

$builder->create_build_script();
