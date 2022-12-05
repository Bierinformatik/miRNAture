package MiRNAture::HMMsearch;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(searchHomologyHMM cmsearch_specific_sequence define_final_CMs runNhmmer obtainTrueCandidates update_coordinates_to_genome check_format_input_header);
use Moose::Role;
use Data::Dumper;
use MiRNAture::ToolBox;
with 'MiRNAture::Merging';
with 'MiRNAture::ToolBox';
with 'MiRNAture::Evaluate';

my %len;

sub searchHomologyHMM {
	my ($genome, $species, $outHMM, $current_HMM_models, $CM_path, $bitscores, $len_r, $names_r, $families_names, $nhmmer_path, $cmsearch_path, $zvalue, $minBitscore, $maxthreshold, $list_fam_file) = @_;
	if (!-d "$outHMM/$species"){ #Check if already exists species-specific directory
		create_folders($outHMM,$species);#Create folder specific to species
	}
	runNhmmer_parallel($current_HMM_models, $genome, $species, $outHMM, $nhmmer_path, $zvalue, $list_fam_file);
	return;
}

sub searchStructureHMM {
	my ($hmm, $genome, $species, $outHMM, $current_HMM_models, $CM_path, $bitscores, $len_r, $names_r, $families_names, $nhmmer_path, $cmsearch_path, $zvalue, $minBitscore, $maxthreshold, $list_fam_file) = @_;
	create_folders("$outHMM/$species","Infernal");#Create folder specific to species
	my $infernal_out_path = "$outHMM/$species/Infernal";
	if (-z "$outHMM/$species/$species.$hmm.tab.true.table.fasta" || !-e "$outHMM/$species/$species.$hmm.tab.true.table.fasta"){
		# Were not geneneated candidates
		return;	
	} else {
	$zvalue = $zvalue/2; # Because it is only one strand
	cmsearch_specific_sequence($hmm, $CM_path, $infernal_out_path, $species, $cmsearch_path, $zvalue);
	my @result_cmseach = check_folder_files($infernal_out_path, "tab");
	my $molecule = get_family_name($hmm, $families_names);
	create_folders("$outHMM/$species/Infernal","Final");#Create folder specific to species
	classify_2rd_align_results($species, $hmm, $infernal_out_path, "$infernal_out_path/$species.$hmm.tab" ,"hmm", $molecule, $bitscores, $len_r, $names_r, $minBitscore, $maxthreshold); #Obtain true candidates
	}
	return;
}

sub cmsearch_specific_sequence_parallel {
	my ($cm_path, $out_path_infernal, $genome_tag, $cmsearch_path, $Zscore) = @_;
	existenceProgram($cmsearch_path);
	# Here have to resolve the path of CMs:  <25-04-22, cavelandiah> #
	system("parallel --will-cite --rpl '.. s/(\.\.\/|\/.*\/)([A-Za-z]+\.)(.*)(\.tab\.true\.table\.fasta)/$3/;' $cmsearch_path -E 0.015 --notrunc -Z $Zscore --nohmmonly --noali --toponly --tblout $out_path_infernal/${genome_tag}.{..}.tab {} 1> /dev/null ::: $out_path_infernal/../$genome_tag.*.tab.true.table.fasta");
	return;
}

sub cmsearch_specific_sequence {
	my ($hmm,$cm_models_path, $out_path_infernal, $genome_tag, $cmsearch_path, $Zscore) = @_;
	existenceProgram($cmsearch_path);
	my $cm = $hmm;
	my $fasta = "$out_path_infernal/../$genome_tag.$hmm.tab.true.table.fasta";
	foreach my $cm_models_path_specific (@$cm_models_path){
		my $model = "$cm_models_path_specific/${cm}.cm";
		if (-e $model && !-z $model){
			my $param = "--cpu 5 -E 0.015 --notrunc -Z $Zscore --nohmmonly --noali --toponly --tblout $out_path_infernal/${genome_tag}.${cm}.tab $model $fasta";
			system "$cmsearch_path $param 1> /dev/null";
		}
	}
	return;
}

sub define_final_CMs {
	my ($outpath, $species, $str, $len_r) = @_;
	my @result_cmsearch;
	unless (-d $outpath){
		print_result("I did not found true candidates for the CM evaluation on HMM searches");
		return;
	}
	if ($str =~ /^hmm$/){
		@result_cmsearch = check_folder_files($outpath, "true\\.table");
	} elsif ($str =~ /^rfam$/){
		@result_cmsearch = check_folder_files($outpath, "true\\.table");
		concatenate_true_cand($species, $outpath,\@result_cmsearch, $str); #Concantenate all true
		resolve_mergings($species, $outpath, "3", $str);
		return;
	} elsif ($str =~ /^mirbase$|^user$/){
		@result_cmsearch = check_folder_files($outpath, "true\\.table");
		concatenate_true_cand($species, $outpath,\@result_cmsearch, $str); #Concantenate all true
		resolve_mergings($species, $outpath, "3", $str);
		return;
	} else {
		@result_cmsearch = check_folder_files($outpath, "$species\_$str\\..*\\.tab\\.true\\.table");
	}
	my $numFiles = scalar @result_cmsearch;
	if ($numFiles == 0){
		print_result("I did not found true candidates for the CM evaluation on HMM searches");
		return;
	} else {
		concatenate_true_cand($species, $outpath, \@result_cmsearch, $str); #Concantenate all true
		update_coordinates_to_genome($outpath, $species, $str, $len_r);
		resolve_mergings($species, $outpath, "1", $str);
	}
	return;
}

sub runNhmmer {
	#In: out, CM model, genome
	my ($genome, $HMM, $species, $out_folder, $path_hmm, $nhmmer_path, $zvalue) = @_;
	existenceProgram($nhmmer_path);
	my $param;
	foreach my $hmms_path_specific (@$path_hmm){
		if (-e "$hmms_path_specific/$HMM.hmm"){
			$param = "--cpu 8 -E 0.015 -Z $zvalue --noali --tblout $out_folder/$species/$species.$HMM.tab $hmms_path_specific/$HMM.hmm $genome";
			system "$nhmmer_path $param 1> /dev/null";
		} else {
			next;
		}
	}
	return;
}

sub runNhmmer_parallel {
	#In: out, CM model, genome
	my ($path_hmm, $genome, $species, $out_folder, $nhmmer_path, $zvalue, $list_cm_file) = @_;
	existenceProgram($nhmmer_path);
	my $param;
	foreach my $hmms_path_specific (@$path_hmm){
		next if $hmms_path_specific =~ /^$|^NA$/;
		system("parallel --will-cite $nhmmer_path -E 0.01 -Z $zvalue --noali --tblout $out_folder/$species/$species\.{/.}\.tab\.true\.table $hmms_path_specific/{/.}\.hmm $genome 1> /dev/null ::: $hmms_path_specific/*\.hmm");
	}
	return;
}

sub obtainTrueCandidates {
	my ($species, $out_path, $hmm_query) = @_;
	my $threshold_hmm = 0.01;
	if (!-e "$out_path/$species.$hmm_query.tab"){
		return;
	}
	open my $IN, "< $out_path/$species.$hmm_query.tab";
	open my $OUTTTRUE, "> $out_path/$species.$hmm_query.tab.true.table" or die;
	open my $OUTFALSE, "> $out_path/$species.$hmm_query.tab.false.table" or die;
	while (<$IN>){
		chomp;
		next if $_ =~ /^\#/;
		my @spl = split /\s+|\t/, $_;
		if ($spl[12] > $threshold_hmm){
			print $OUTFALSE "$_\n";
		} else {
			print $OUTTTRUE "$_\n";
		}
	}
	close $OUTTTRUE; close $OUTFALSE; close $IN;
	return;
}

sub update_coordinates_to_genome { #Only for HMMs and Blast strategies
	my ($input_folder, $species, $strategy, $len_r) = @_;
	my ($IN, $OUT);
	if ($strategy =~ /^hmm$/){
		open $IN, "< $input_folder/all_RFAM_$species.truetable" or die;
		open $OUT, "> $input_folder/all_RFAM_$species.truetable.clean" or die;
	} else {
		open $IN, "< $input_folder/all_RFAM_${species}_${strategy}.truetable" or die;
		open $OUT, "> $input_folder/all_RFAM_${species}_${strategy}.truetable.clean" or die;
	}
	#dvexscf65180.-.7909.8051.RF00001.3.117.119      -       RF00001 -       cm      1       119     14      132     +       no      1       0.58    0.0     130.5   2.5e-30 !       -
	my ($realEnd, $realStart);
	my (@parts, @other);
	while (<$IN>){
		chomp;
		#scaffold1545-size16374.-.2678.2796.RF00001.3.117.119	-	RF00001	-	cm	1	119	1	119	+	no	1	0.58	0.0	130.5	9.4e-32	!	-	119
		#JH126831.1.+.75490.75749.14      -         RF01021              -          cm        1       94       99    190 +    no    1 0.21   8.8   34.9   2.7e-07 !   -
		@parts = split /\s+|\t/, $_;
		@other = split /\./, $parts[0]; #Header
		@other = check_format_input_header(\@other);
		my $distance = abs($other[3] - $other[2]) + 1;
		my $sense = $other[1];
		if($sense eq "+"){
			if($parts[7] >= $parts[8]){
				die "This makes no sense, negative is not correct, hasn't evaluated!\n";
			} else {
				$realStart = ($other[2] + $parts[7]);
				$realEnd = ($other[2] + $parts[8]);
			}		
		} elsif ($sense eq "-"){
			if($parts[7] >= $parts[8]){
				die "This makes no sense, negative is not correct, hasn't evaluated!\n";
			} else {
				$realStart = ($other[2] + $parts[7]);
				$realEnd = ($other[2] + $parts[8]);
			}
		}
		my $len_cm = $$len_r{$parts[3]};
		# Here prints the line with the group reference of the merged query candidates
		# on the blast block file species_str.miRNA.tab.db.location.database: $parts[-1]
		if ($strategy =~ /^hmm$/){
			print $OUT "$other[0]\t$other[4]\t$other[1]\t$parts[5]\t$parts[6]\t$realStart\t$realEnd\t$parts[14]\t$parts[15]\t$parts[10]\t$parts[2]\t$len_cm\t$strategy\n";
		} else {
			print $OUT "$other[0]\t$other[-1]\t$other[1]\t$parts[5]\t$parts[6]\t$realStart\t$realEnd\t$parts[14]\t$parts[15]\t$parts[10]\t$parts[2]\t$len_cm\t$strategy\n";
		}
	}
	close $IN;
	close $OUT;
	return;
}

sub check_format_input_header {
	my $header_array = shift;
	#my @header_array =  split /\./, $header;
	my $len = scalar @$header_array;
	if ($$header_array[1] =~ /^[0-9]/ && $len >= 5){ #Discriminate "NW_003107091.1.+.410955.411080.xtr101.21.89.108" and accept "chr1.-.3351215.3351292.cin280.1.53.53"
		splice @$header_array, 1, 1;
	} else {
		;
	}
	return @$header_array;
}

no Moose::Role;
1;
