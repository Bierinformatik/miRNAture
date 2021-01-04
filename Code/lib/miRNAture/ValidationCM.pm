package miRNAtureValidationCM;

use Moose;
use Data::Dumper;

has 'name_model' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'length_model' => (
	is => 'ro',
	isa => 'Int',
	required => 1,
);

has 'acc_model' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'score_model' => (
	is => 'ro',
	isa => 'Int',
);

has 'family_model' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

no Moose;
__PACKAGE__->meta->make_immutable;
