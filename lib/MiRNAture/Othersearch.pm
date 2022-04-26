package MiRNAture::Othersearch;

use Moose::Role;
use Data::Dumper;

with 'MiRNAture::ToolBox';

sub searchOthershomology {
	my $shift = shift;
	my $zscore = shift;
	my $minBitscore = shift;
	my $maxthreshold = shift;
	my $models = $shift->path_covariance;
	##foreach my $pathCM (@$models){
	cmsearch_mirbase_parallel($shift->genome_subject, $shift->subject_specie, $shift->output_folder."/".$shift->subject_specie, $models, 1, $shift->cmsearch_program_path->stringify, $zscore);
	##cmsearch($shift->cm_model, $shift->genome_subject, $shift->subject_specie, $shift->output_folder."/".$shift->subject_specie, $pathCM, 1, $shift->cmsearch_program_path->stringify, $zscore);
	##}
	return;
}

sub cmsearch_mirbase_parallel {
	my ($genome, $genomeTag, $outFolder, $path_cm, $mode, $cmsearch_path, $zscore) = @_; #modes: 1 X.cm 0 other name
	existenceProgram($cmsearch_path);
	$genome =~ s/"//g;
	foreach my $path_cm_specific (@$path_cm){
		next if $path_cm_specific =~ /^$/;
		next if $path_cm_specific =~ /RFAM/;
		system("parallel $cmsearch_path -E 0.015 --notrunc -Z $zscore --nohmmonly --tblout $outFolder/{/.}_$genomeTag.tab -o $outFolder/{/.}_$genomeTag.out $path_cm_specific/{/.}\.cm $genome 1> /dev/null ::: $path_cm_specific/MIPF*\.cm");
	}
	return;
}


no Moose::Role;
1;
