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
	foreach my $pathCM (@$models){
		cmsearch($shift->cm_model, $shift->genome_subject, $shift->subject_specie, $shift->output_folder."/".$shift->subject_specie, $pathCM, 1, $shift->cmsearch_program_path->stringify, $zscore);
	}
	my $molecule = "NA";
	classify_2rd_align_results($shift->subject_specie, $shift->cm_model, $shift->output_folder."/".$shift->subject_specie, $shift->output_folder."/".$shift->subject_specie."/".$shift->cm_model."\_".$shift->subject_specie."\.tab", "OTHER_CM", $molecule, $shift->bitscores_CM, $shift->length_CM, $shift->names_CM, $minBitscore, $maxthreshold);
	return;
}

no Moose::Role;
1;
