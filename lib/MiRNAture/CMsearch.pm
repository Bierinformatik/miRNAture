package MiRNAture::CMsearch;

use Moose::Role;
use Data::Dumper;
with 'MiRNAture::CMsearch';
with 'MiRNAture::ToolBox';
with 'MiRNAture::Evaluate';

sub searchCMhomology {
	my ($shift, $zscore,$minBitscore, $maxthreshold) = @_;
	cmsearch_rfam_parallel($shift->genome_subject, $shift->subject_specie, $shift->output_folder."/".$shift->subject_specie, $shift->path_covariance, 1, $shift->cmsearch_program_path->stringify, $zscore);
	##cmsearch($shift->cm_model, $shift->genome_subject, $shift->subject_specie, $shift->output_folder."/".$shift->subject_specie, $shift->path_covariance, 1, $shift->cmsearch_program_path->stringify, $zscore);
	return;
}

sub cmsearch_rfam_parallel {
	my ($genome, $genomeTag, $outFolder, $path_cm, $mode, $cmsearch_path, $zscore) = @_; #modes: 1 X.cm 0 other name
	existenceProgram($cmsearch_path);
	$genome =~ s/"//g;
	foreach my $path_cm_specific (@$path_cm){
		next if $path_cm_specific =~ /^$/;
		system("parallel $cmsearch_path -E 0.015 --notrunc -Z $zscore --nohmmonly --tblout $outFolder/{/.}_$genomeTag.tab -o $outFolder/{/.}_$genomeTag.out $path_cm_specific/{/.}\.cm $genome 1> /dev/null ::: $path_cm_specific/RF*\.cm");
	}
	return;
}

sub cmsearch {
	#In: out, CM model, genome
	my ($nameCM, $genome, $genomeTag, $outFolder, $path_cm, $mode, $cmsearch_path, $zscore) = @_; #modes: 1 X.cm 0 other name
	existenceProgram($cmsearch_path);
	my $nameCMFinal;
	if($mode == 1){ #From other source
		$nameCMFinal = "${nameCM}.cm";
	} elsif ($mode == 0){ #From metazoa
		$nameCMFinal = "${nameCM}_metazoa.cm";
	}
	$genome =~ s/"//g;
	foreach my $path_cm_specific (@$path_cm){
		if (-e "$path_cm_specific/${nameCMFinal}" && !-z "$path_cm_specific/${nameCMFinal}"){
			my $param = "--cpu 5 -E 0.015 --notrunc -Z $zscore --nohmmonly --tblout $outFolder/${nameCM}_$genomeTag.tab -o $outFolder/${nameCM}_$genomeTag.out $path_cm_specific/${nameCMFinal} $genome";
			system "$cmsearch_path $param 1> /dev/null";
		} 
	}
	return;
}

no Moose::Role;
1;
