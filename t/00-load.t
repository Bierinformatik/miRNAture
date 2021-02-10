#!perl -T
use 5.26.2;
use strict;
use warnings;
use Test::More;

plan tests => 31;

BEGIN {
    use_ok( 'Bio::miRNAture' ) || print "Bail out!\n";
    use_ok( 'Mir::ConfigurationFile' ) || print "Bail out!\n";
    use_ok( 'Mir::Eval' ) || print "Bail out!\n";
    use_ok( 'Mir::miRNAnchor' ) || print "Bail out!\n";
    use_ok( 'Mir::miRNAture' ) || print "Bail out!\n";
    use_ok( 'MiRNAnchor::Check' ) || print "Bail out!\n";
    use_ok( 'MiRNAnchor::Classify' ) || print "Bail out!\n";
    use_ok( 'MiRNAnchor::ExternalPrograms' ) || print "Bail out!\n";
    use_ok( 'MiRNAnchor::Main' ) || print "Bail out!\n";
    use_ok( 'MiRNAnchor::Tools' ) || print "Bail out!\n";
    use_ok( 'MiRNAnchor::Validate' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::Blast' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::BlastPrepareQueries' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::_BlastSearch_BlockDetection' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::_BlastSearch' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::_BlastSearch_ResolveBlastMergings' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::Cleaner' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::CM' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::CMsearch' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::ConfigFile' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::Evaluate' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::FinalCandidates' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::HMM' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::HMMsearch' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::LogFile' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::Merging' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::Othersearch' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::Others' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::ToolBox' ) || print "Bail out!\n";
    use_ok( 'MiRNAture::ValidationCM' ) || print "Bail out!\n";
    use_ok( 'Mir::Summarise' ) || print "Bail out!\n";
}

diag( "Testing Bio::miRNAture $Bio::miRNAture::VERSION, Perl $], $^X" );
