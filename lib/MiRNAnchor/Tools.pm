package MiRNAnchor::Tools;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(check_folder_files2 load_genomes_location_file create_folder create_folder_environment copy_files copy_folders change_name create_genome_list_mirfix include_subject_genome setup_all_to_run_mirfix_group write_mirfix_specific_family_individual run_mirfix_individual wait_slurm_process check_if_file_exists print_error print_result print_process filter_wrong_alignments detect_multifamily_rfam setup_all_to_run_mirfix_group_subset run_mirfix_individual_subset load_correspondence_models_rfam read_data_share end_close evaluated_family);

use Moose::Role;
use Data::Dumper;
use File::Copy;
use File::Copy::Recursive qw(dircopy);
use Term::ANSIColor qw(:constants); 
use MiRNAnchor::ExternalPrograms;
use File::Find;

=head1 check_folder_files2
	Title: check_folder_files
	Usage: check_folder_files(<PATH_DIR_TO_SEARCH>, PATTERN);
	Function: Get list files that have specific pattern. 
	Returns: List of files in selected folder that contains the
		menctioned pattern in its name. 
=cut 

sub check_folder_files2 {
	my ($dir, $prefix) = @_;
	opendir(DIR, $dir);
	my @files = grep(/$prefix$/, readdir(DIR));
	closedir(DIR);
	return @files;
}

=head1 load_genomes_location_file 
	Title: load_genomes_location_file
	Usage: load_genomes_location_file(genome_file);
	Function: Get list files that have specific pattern. 
	Returns: List of files in selected folder that contains the
		menctioned pattern in its name. 
=cut 

sub load_genomes_location_file {
	my $genome_list_file = shift;
	open my $GENOMESPATH, "< $genome_list_file" or die;
	my @complete_genomes_path;
	while (<$GENOMESPATH>){
		chomp;
		push @complete_genomes_path, $_;
	}
	return \@complete_genomes_path;
}

=head1 create_folder
	Title: create_folder
	Usage: create_folder(FOLDER_TO_CREATE)
	Function: Create folder from perl, using mkdir.
	Returns: Create desired folder, with chmod 755.

=cut 

sub create_folder {
	my $dir = shift;
	if (!-d $dir){ #If DIR doesn't exists
		mkdir($dir, 0755);
	}
	return;
}

=head1 create_folder_environment
	Title: create_folder_environment
	Usage: create_folder_environment(FOLDER_TO_CREATE)
	Function: Create folder from perl, using mkdir.
	Returns: Create desired folder, with chmod 755.

=cut 

sub create_folder_environment {
	my ($working_path, $mode, $class, $acc_RFAM, $acc_mirbase) = @_;
	my ($dir, $dir2, $dir3);
	my $dir0 = "$working_path/$acc_RFAM";
	create_folder($dir0);
	if ($mode eq "RFAM"){
		my $dir1 = "$working_path/$acc_RFAM/$class";
		create_folder($dir1);
		my $dir2 = "$working_path/$acc_RFAM/$class/$acc_mirbase";
		create_folder($dir2);
		my @folders = ("BaseFiles","Families","Output"); 
		foreach my $folder (@folders){
			$dir3 = "$working_path/$acc_RFAM/$class/$acc_mirbase/$folder";
			create_folder($dir3);
		}
	} elsif ($mode eq "Validation"){
		my @folders = ("BaseFiles","Families","Output"); 
		foreach my $folder (@folders){
			$dir3 = "$working_path/$acc_RFAM/$folder";
			create_folder($dir3);
		}
	} elsif ($mode eq "Subset"){
		my @folders = ("BaseFiles","Families","Output"); 
		foreach my $folder (@folders){
			$dir3 = "$working_path/$acc_RFAM/$folder";
			create_folder($dir3);
		}
	} else {
		die "Not yet implemented!\n";
	}
	return;
}

=head1 copy_files 
	Title: copy_files
	Usage: copy_files(ORIGIN, DESTINITY)
	Function: Copy files from Perl.
	Returns: Copied file, in the indicated path.

=cut 

sub copy_files {
	my ($old, $location) = @_;
	my $new_name = change_name($old);
	my $new_copy = "$location/$new_name";
	copy($old, $new_copy) or die "Copy failed for $old or $location";
	return;
}

sub copy_folders {
	my ($path, $new) = @_;
	dircopy($path, $new); # Copy $dir1 to $dir2 recursively
	return;
}

=head1 change_name
	Title: copy_files
	Usage: copy_files(ORIGIN, DESTINITY)
	Function: Copy files from Perl.
	Returns: Copied file, in the indicated path.

=cut 

sub change_name {
	my $name = shift;
	$name =~ s/^(\/.*\/|\.\.\/|\.\/|.*\/)(.*)$/$2/g;
	my $new_name = $name;
	if ($name =~ m/\_mapping\_/ && $name !~ m/\_mapping\_original/){
		$new_name =~ s/^(.*)(\_mapping\_updated.txt)$/$1.mapping_original.txt/g;# RF00103.mapping_original.txt RF00103_mapping_updated.txt
	} elsif ($name =~ m/\_mature\_/ && $name !~ m/\_mature\_original/){
		$new_name =~ s/^(.*)(_mature_updated.fa)$/$1_mature_original.fa/g;
	} elsif ($name =~ m/\.mapping\_original/){
		$new_name =~ s/^(.*)(\.mapping_original.txt)$/$1.mappingtest.txt/g;
	} elsif ($name =~ m/\_mature\_original/){
		$new_name =~ s/^(.*)(_mature_original.fa)$/$1_maturetest.fa/g;
	} elsif ($name =~ m/Final\.fasta/){
		$new_name =~ s/(.*)(\-Final\.fasta)/$1.fa/g;
	} else {
		;
	}
	return $new_name;
}

=head1 create_genome_list_mirfix
	Title: create_genome_list_mirfix
	Usage: create_genome_list_mirfix()
	Function: Assign the list of genomes required to run MIRFix.
	Returns: At the pointed out folder, a new created genome file is created.

=cut 

sub create_genome_list_mirfix {
	my ($working_path, $destination, $name_tag, $genome_path_list, $specie_pattern) = @_; 
	my @candidates = grep(m/$specie_pattern/, @$genome_path_list); 
	my $outT = "$working_path/$destination/${name_tag}_genomes_list.txt";
	open my $OUT6, ">> $outT" or die;
	foreach my $ln (@candidates){
		print $OUT6 "$ln\n";
	}
	close $OUT6;
	return;
}

=head1 filter_wrong_alignments
	Title: filter_wrong_alignments
	Usage: filter_wrong_alignments
	Function: Detect the CMs identified as misleading and avoid validation 
    Returns: Warning message
=cut 

sub filter_wrong_alignments  {
	my $model = shift;
	my @filtered = ("RF00929","RF00178", "RF00639", "RF00668", "RF00672", "RF00753", "RF00788", "RF00866", "RF00900", "RF00910", "RF00990", "RF01011", "RF01916");
	foreach my $modelAvoid (@filtered){
		if ($modelAvoid =~ /^$model$/){
			return 2;
		} 
	}
	return 0;
}

=head1 detect_multifamily_rfam
	Title: detect_multifamily_rfam
	Usage: detect_multifamily_rfam
	Function: Detect miRNA families from RFAM splitted in two sets 
    Returns: Positive or negative reporter variable
=cut 

sub detect_multifamily_rfam {
	my $model = shift;
	my @filtered = ("RF00178", "RF00639", "RF00668", "RF00711");
	foreach my $modelAvoid (@filtered){
		if ($modelAvoid =~ /^$model$/){
			return 1;
		} 
	}
	return 0;
}

sub include_subject_genome {
	my ($working_path, $destination, $name_tag, $genome_path_list, $new_genome) = @_; 
	my $outT = "$working_path/$destination/${name_tag}_genomes_list.txt";
	open my $OUT6, ">> $outT" or die;
	print $OUT6 "$new_genome\n";
	close $OUT6;
	return;
}

sub setup_all_to_run_mirfix_group {
	my ($location, $name, $MIRFIX_path, $current_dir, $gridengine, $code, $address, $tagSpe) = @_;
	if (-e "$location/BaseFiles/${name}_genomes_list.txt" && !-z "$location/BaseFiles/${name}_genomes_list.txt"){ #Copied mirfix files, skip copy
		print_process("RNA structure evaluation for: $name-$code");
	} else {
		print_error("The required input files to validate the miRNAs are missing, check input files from Rfam");
	}
	##Parameters constructor
	my $param_mirfix = MiRNAnchor::ExternalPrograms->new(
		cores => "20",
		output_path => "$location/Output/",
		families_path => "$location/Families/",
		list_file => "$location/BaseFiles/${name}_list.txt",
		genomes_file => "$location/BaseFiles/${name}_genomes_list.txt",
		mapping_file => "$location/BaseFiles/${name}.mappingtest.txt",
		mature_file => "$location/BaseFiles/${name}_maturetest.fa",
		extension => "10",
		log_level => "DEBUG"
	);
	my $parameters = $param_mirfix->build_parameters;	
	if ($gridengine == 0){
		#Run without SLURM
		open my $OUTTIME, "> $location/time_$name.dat" or die;
		my $time = run_mirfix($parameters, $MIRFIX_path);
		clean_mirfix_logs($current_dir, $location, $name);
		print $OUTTIME "$name: $time s\n";
		close $OUTTIME;
	} elsif ($gridengine == 1){
		write_mirfix_specific_family_individual($current_dir, $name, $parameters, $MIRFIX_path, $code, $address, $tagSpe);
	}
	return;
}

sub write_mirfix_specific_family_individual {
	my ($dir_now, $name, $param, $MIRFIX_path, $code, $outAddress, $tagSpe) = @_;
    ##create_folder($dir_now); # Create base working dir
    ##create_folder("$dir_now/$outAddress"); #Create folder specific to specie
    ##create_folder("$dir_now/$outAddress/LOGs"); #Create folder specific to specie
	my $OUTTEMPC;
	if (-e "$outAddress/${name}_${code}.sh"){
		open $OUTTEMPC, "> $outAddress/${name}_${code}_${tagSpe}.sh";
	} else {
		open $OUTTEMPC, "> $outAddress/${name}_${code}_${tagSpe}.sh";
		print $OUTTEMPC "\#!\/bin\/bash\n";
	}
	print $OUTTEMPC "python3 $MIRFIX_path $param\n";
	chmod 0755, $OUTTEMPC;
	close $OUTTEMPC;
	return;
}

sub setup_all_to_run_mirfix_group_subset {
	my ($location, $name, $MIRFIX_path, $current_dir, $gridengine, $code, $address, $tagSpe, $subset) = @_;
	if (-e "$location/BaseFiles/${subset}_genomes_list.txt" && !-z "$location/BaseFiles/${subset}_genomes_list.txt"){ #Copied mirfix files, skip copy
		print_process("$name-$code input files already copied");
	} else {
		print_error("The required input files to validate the miRNAs are missing, check input files from Rfam");
		#ExternalPrograms::copyStartFiles("$location/BaseFiles", $name);
	}
	##Parameters constructor
	my $param_mirfix = MiRNAnchor::ExternalPrograms->new(
		cores => "20",
		output_path => "$location/Output/",
		families_path => "$location/Families/",
		list_file => "$location/BaseFiles/${subset}_list.txt",
		genomes_file => "$location/BaseFiles/${subset}_genomes_list.txt",
		mapping_file => "$location/BaseFiles/${subset}.mappingtest.txt",
		mature_file => "$location/BaseFiles/${subset}_maturetest.fa",
		extension => "10",
		log_level => "DEBUG"
	);
	my $parameters = $param_mirfix->build_parameters;	
	if ($gridengine == 0){
		#Run without SLURM
		open my $OUTTIME, "> $location/time_$subset.dat" or die;
		my $time = run_mirfix($parameters, $MIRFIX_path);
		clean_mirfix_logs($current_dir, $location, $subset);
		print $OUTTIME "$name: $time s\n";
		close $OUTTIME;
	} elsif ($gridengine == 1){
		write_mirfix_specific_family_individual_subset($current_dir, $name, $parameters, $MIRFIX_path, $code, $address, $tagSpe, $subset);
	}
	return;
}

sub write_mirfix_specific_family_individual_subset {
	my ($dir_now, $name, $param, $MIRFIX_path, $code, $outAddress, $tagSpe, $subset) = @_;
	create_folder($dir_now); # Create base working dir
	create_folder("$dir_now/$outAddress"); #Create folder specific to specie
	create_folder("$dir_now/$outAddress/LOGs"); #Create folder specific to specie
	my $OUTTEMPC;
	if (-e "$outAddress/${name}_${code}_${tagSpe}_$subset.sh"){
		open $OUTTEMPC, "> $outAddress/${name}_${code}_${tagSpe}_$subset.sh";
	} else {
		open $OUTTEMPC, "> $outAddress/${name}_${code}_${tagSpe}_$subset.sh";
		print $OUTTEMPC "\#!\/bin\/bash\n";
	}
	print $OUTTEMPC "python3 $MIRFIX_path $param\n";
	chmod 0755, $OUTTEMPC;
	close $OUTTEMPC;
	return;
}

sub run_mirfix_individual {
	my ($name, $dir_now, $code,$outAddress,$tagSpe) = @_;	
	my $short_spe = substr($tagSpe, 0, 3);
	my $short = substr($name,-5);
	$short = "${short_spe}${short}";
	system "sbatch --job-name=${short} --nodes=1 --ntasks=1 --cpus-per-task=20 --time=24:00:00 --mem=2G --output=$outAddress/LOGs/out_${name}_${code}_${tagSpe}.out --error=$outAddress/LOGs/error_${name}_${code}_${tagSpe}.out $outAddress/${name}_${code}_${tagSpe}.sh &";	
	#system "qsub -S /bin/sh -cwd -l h_vmem=4G -pe smp 5 -m bes -V -q normal.q -o $dir_now/$outAddress/LOGs/out_${name}_${code}.out -e $dir_now/$outAddress/LOGs/error_${name}_${code}.out $dir_now/$outAddress/${name}_${code}.sh &";
	sleep(2); #I detected some delay on the sge system while I summited the job and it started to be recognized as job on job table
	return;
}

sub run_mirfix_individual_subset {
	my ($name, $dir_now, $code, $outAddress, $tagSpe, $subset) = @_;	
	system "sbatch --job-name=${subset}_${code}_${tagSpe}_$subset --nodes=1 --ntasks=1 --cpus-per-task=20 --time=24:00:00 --mem=2G --output=$outAddress/LOGs/out_${subset}_${code}_${tagSpe}_$subset.out --error=$outAddress/LOGs/error_${subset}_${code}_${tagSpe}_$subset.out $outAddress/${subset}_${code}_${tagSpe}_$subset.sh &";	
	#system "qsub -S /bin/sh -cwd -l h_vmem=4G -pe smp 5 -m bes -V -q normal.q -o $dir_now/$outAddress/LOGs/out_${name}_${code}.out -e $dir_now/$outAddress/LOGs/error_${name}_${code}.out $dir_now/$outAddress/${name}_${code}.sh &";
	sleep(2); #I detected some delay on the sge system while I summited the job and it started to be recognized as job on job table
	return;
}

sub wait_slurm_process {
	my $reference = shift;	
	my $specie = shift;
	$reference = substr($reference,-5);
	my $short_spe = substr($specie, 0, 3);
	$reference = "${short_spe}${reference}";
	print "Waiting for the process on $reference\n";
	if (length $reference > 8){
		$reference = substr($reference, 0, 8);
	}
	EVAL:
	my $state_queue = `squeue | grep $reference`; #Capture state
	my $state = 0;
	if ($state_queue && length $state_queue > 0){
		$state = 1;
		sleep(5); #Sleep 60 s.
		goto EVAL;
	}
	$state_queue = ();
	return;
}	

sub check_if_file_exists {
	my ($file) = shift;
	my $answer;
	if (-e $file && !-z $file){ #Exists
		$answer = 1;
	} else { #Non-exists
		$answer = 0;
	}
	return $answer;
}

sub print_error {
	my $text = shift;
	local $Term::ANSIColor::AUTORESET = 1;
	print BOLD RED "[ERROR] $text\n";
	return;
}

sub print_result {
	my $text = shift;
	local $Term::ANSIColor::AUTORESET = 1;
	print BOLD GREEN "[RESULT] $text\n";
	return;
}

sub print_end {
	my $text = shift;
	local $Term::ANSIColor::AUTORESET = 1;
	print BOLD BRIGHT_WHITE"[END] $text\n";
	return;
}

sub print_process {
	my $text = shift;
	local $Term::ANSIColor::AUTORESET = 1;
	print FAINT BLUE "[PROCESS] $text\n";
	return;
}

=head1 load_correspondence_models_rfam 
	Title: load_correspondence_models_rfam
	Usage: load_correspondence_models_rfam(<path_miRNAture_code>);
	Function: Sort all out to validate RFAM or other families into specific precalculated models. The main files
	are located in the ~/miRNAture/Data/ folder: <mirbase_rfam_correspondence.txt> and <user_correspondence_families.txt>
	that could be modified by the user to provide another relations.
	Returns: Change the corresponding validation family into the database table all_RFAM_*_Final.ncRNAs_homology.txt.db \
	    and fasta files.
=cut

sub load_correspondence_models_rfam {
	my $base_path = shift;
	my $user_path = shift;
	my $models_file = "$base_path/mirbase_rfam_correspondence.txt";
	open my $IN, "< $models_file" or die "The file mirbase_rfam_correspondence.txt is required to perform annotation of mature's.\n";
	# <RFAM_ACC> <MIRBASE_ACC_FAMILY>
	my %db;
	while (<$IN>){
		chomp;
		my @all = split /\s+|\t/, $_;
		$db{$all[0]} = $all[-1]; #Rfam->miRBase
		$db{$all[-1]} = $all[-1]; #miRBase->miRBase
		#TODO: include multiple families to test for the same homology model
		#push @{$db{$all[0]}}, $all[-1]; #Rfam->miRBase
		#push @{$db{$all[-1]}}, $all[-1]; #miRBase->miRBase
	}
	close $IN;
	if (length $user_path !~ /^NO$/){
		if (-e "$user_path/user_correspondence_families.txt" && !-z "$user_path/user_correspondence_families.txt"){
			open my $IN3, "< $user_path/user_correspondence_families.txt" or die "The file $user_path/user_correspondence_families.txt is broken, please create it\n<ModelACC> <miRBase_ACC>\n\n";
			while (<$IN3>){
			       	chomp;
			       	my @all = split /\s+|\t/, $_;
			       	$db{$all[0]} = $all[-1];
		       	}
		       	close $IN3;
	       	}
       	}
	return \%db;
}

=head1 create_folders 
    Title: create_folders 
    Usage: create_folders(base_folder, new_folder);
    Function: Creates a new folder based on the given path. 
    Returns: Creates a new folder if not exists. 
=cut 

sub create_folders {
	my ($in, $new) = @_;
	my $dir = "$in/$new";
	if (!-d $dir){ #If DIR doesn't exists
		mkdir($dir, 0755);
	}
	return;
}

=head1 evaluated_family 
    Title: evaluated_family 
    Usage: evaluated_family(mirna_fam);
    Function: Check if folder exists. 
    Returns: '1' if does not exists, otherwise '0'.
=cut 

sub evaluated_family {
	my ($fam, $outpath_eval) = @_;
	my $dir = "$outpath_eval/$fam";
	if (-d $dir){ #If DIR does exists
		return 1;
	}
	return 0;
}
#sub read_data_share {
#	my $data_location = dist_dir('Bio-miRNAture');
#	print "$data_location\n";
#	return $data_location;
#}

sub end_close {
	print_end("\n-¿Olvida usted algo?-\n ¡Ojalá!.\n\n  El emigrante. Luis Felipe Lornelí. (2005)\n");
	return;
}
no Moose::Role;
1;
