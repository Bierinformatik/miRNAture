#!perl -T
use 5.26.2;
use strict;
use warnings;
use Test::More;

plan tests => 33;

sub not_in_file_ok {
    my ($filename, %regex) = @_;
    open( my $fh, '<', $filename )
        or die "couldn't open $filename for reading: $!";

    my %violated;

    while (my $line = <$fh>) {
        while (my ($desc, $regex) = each %regex) {
            if ($line =~ $regex) {
                push @{$violated{$desc}||=[]}, $.;
            }
        }
    }

    if (%violated) {
        fail("$filename contains boilerplate text");
        diag "$_ appears on lines @{$violated{$_}}" for keys %violated;
    } else {
        pass("$filename contains no boilerplate text");
    }
}

sub module_boilerplate_ok {
    my ($module) = @_;
    not_in_file_ok($module =>
        'the great new $MODULENAME'   => qr/ - The great new /,
        'boilerplate description'     => qr/Quick summary of what the module/,
        'stub function definition'    => qr/function[12]/,
    );
}

TODO: {
  local $TODO = "Need to replace the boilerplate text";

  not_in_file_ok(README =>
    "The README is used..."       => qr/The README is used/,
    "'version information here'"  => qr/to provide version information/,
  );

  not_in_file_ok(Changes =>
    "placeholder date/time"       => qr(Date/time)
  );

  module_boilerplate_ok('lib/Bio/miRNAture.pm');
  module_boilerplate_ok('lib/Mir/ConfigurationFile.pm');
  module_boilerplate_ok('lib/Mir/Eval.pm');
  module_boilerplate_ok('lib/Mir/miRNAnchor.pm');
  module_boilerplate_ok('lib/Mir/miRNAture.pm');
  module_boilerplate_ok('lib/MiRNAnchor/Check.pm');
  module_boilerplate_ok('lib/MiRNAnchor/Classify.pm');
  module_boilerplate_ok('lib/MiRNAnchor/ExternalPrograms.pm');
  module_boilerplate_ok('lib/MiRNAnchor/Main.pm');
  module_boilerplate_ok('lib/MiRNAnchor/Tools.pm');
  module_boilerplate_ok('lib/MiRNAnchor/Validate.pm');
  module_boilerplate_ok('lib/MiRNAture/Blast.pm');
  module_boilerplate_ok('lib/MiRNAture/BlastPrepareQueries.pm');
  module_boilerplate_ok('lib/MiRNAture/_BlastSearch_BlockDetection.pm');
  module_boilerplate_ok('lib/MiRNAture/_BlastSearch.pm');
  module_boilerplate_ok('lib/MiRNAture/_BlastSearch_ResolveBlastMergings.pm');
  module_boilerplate_ok('lib/MiRNAture/Cleaner.pm');
  module_boilerplate_ok('lib/MiRNAture/CM.pm');
  module_boilerplate_ok('lib/MiRNAture/CMsearch.pm');
  module_boilerplate_ok('lib/MiRNAture/ConfigFile.pm');
  module_boilerplate_ok('lib/MiRNAture/Evaluate.pm');
  module_boilerplate_ok('lib/MiRNAture/FinalCandidates.pm');
  module_boilerplate_ok('lib/MiRNAture/HMM.pm');
  module_boilerplate_ok('lib/MiRNAture/HMMsearch.pm');
  module_boilerplate_ok('lib/MiRNAture/LogFile.pm');
  module_boilerplate_ok('lib/MiRNAture/Merging.pm');
  module_boilerplate_ok('lib/MiRNAture/Othersearch.pm');
  module_boilerplate_ok('lib/MiRNAture/Others.pm');
  module_boilerplate_ok('lib/MiRNAture/ToolBox.pm');
  module_boilerplate_ok('lib/MiRNAture/ValidationCM.pm');
  module_boilerplate_ok('lib/Mir/Summarise.pm');


}

