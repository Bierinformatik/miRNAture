package MiRNAture::Cleaner;

use Moose::Role;
use Data::Dumper;

sub delete_empty_files {
	my $folder = shift;
	if (-d $folder){
		system("find $folder -size 0 -delete");
	}
	return;
}

sub delete_temporal_files {
	my ($folder, $target) = @_;
	system("rm -f $folder/$target");
}


no Moose::Role;
1;
