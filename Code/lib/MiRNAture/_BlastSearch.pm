package MiRNAture::_BlastSearch;

use Moose;
use Data::Dumper;
use List::MoreUtils;
use MiRNAture::_BlastSearch_ResolveBlastMergings;
with 'MiRNAture::HMMsearch';

has 'blast_str' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'output_folder' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	required => 1,
	coerce => 1,
);

has 'query_folder' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	required => 1,
	coerce => 1,
);

has 'genome_subject' => (
	is => 'ro',
	isa => "Str",
	required => 1,
);

has 'subject_specie' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'data_folder' => (
	is => 'ro',
	isa => 'Path::Class::Dir',	
	required => 1,
	coerce => 1
);

has 'current_directory' => (
	is => 'ro',
	isa => 'Path::Class::Dir',	
	required => 1,
	coerce => 1
);

has 'path_covariance' => (
	is => 'ro',
	isa => 'ArrayRef[Str]',
	required => 1,
);

has 'length_CM' => (
	is => 'rw',
	isa => 'HashRef[Str]',
);

has 'names_CM' => (
	is => 'rw',
	isa => 'HashRef[Str]',
);

has 'bitscores_CM' => (
	is => 'rw',
	isa => 'HashRef[Str]',
);

has 'input_path_db_file' => (
	is => 'rw',
	isa => 'Str',
); 

has 'output_merged' => (
	is => 'rw',
	isa => 'Str',
);
has 'output_special' => (
	is => 'rw',
	isa => 'Str',
);

has 'parallel_running' => (
	is => 'ro',
	isa => 'Int',
	required => 1
);

has 'models_list' => (
	is => 'ro',
	isa => 'Path::Class::File',	
	required => 1,
	coerce => 1
);

has 'blast_program_path' => (
	is => 'ro',
	isa => 'Path::Class::File',	
	required => 1,
	coerce => 1
);

has 'cmsearch_program_path' => (
	is => 'ro',
	isa => 'Path::Class::File',	
	required => 1,
	coerce => 1
);

has 'makeblast_program_path' => (
	is => 'ro',
	isa => 'Path::Class::File',	
	required => 1,
	coerce => 1
);

with 'MiRNAture::ToolBox'; 
with 'MiRNAture::BlastPrepareQueries';
with 'MiRNAture::_BlastSearch_ResolveBlastMergings';
with 'MiRNAture::_BlastSearch_BlockDetection';
with 'MiRNAture::Cleaner';

sub searchHomologySequenceBlast {
	my $shift = shift;
	my @id_running;
	my $param = blast_strategies($shift->blast_str); #Get specific parameters.
	create_folders($shift->output_folder, $shift->subject_specie ); #Create folder specific to spe
	create_folders($shift->current_directory,".blastTemp/");
	create_folders($shift->current_directory,".blastTemp/LOGs");
	# Load both score files, from RFAM and others
	my $families = index_ncRNA_families($shift->data_folder."/Basic_files/all_RFAM_scores.txt", $shift->data_folder."/Basic_files/all_other_scores.txt"); #Create groups of ncRNA families from CMs
	my ($molecules, $query_species, $files_relation) = infer_data_queries($shift->query_folder); #Based on metafile, infer query species and molecules
	foreach my $molecule (@$molecules){ #Each fasta file provided by the user
		foreach my $queryS (@$query_species){
			print_process("Searching $molecule homologs from $queryS");
			my $tag_query_spe = create_tag_specie($queryS);
			my $query_seq;
			if (exists $$files_relation{"${molecule}_${queryS}"}){
				my $file_name = $$files_relation{"${molecule}_${queryS}"};
				$query_seq = $shift->query_folder."/".$file_name;
			}
			if (-z $query_seq || !-e $query_seq){
				next;
			} else {
				writeBlastn($param, $shift->genome_subject, $shift->blast_str, $query_seq, $tag_query_spe, $molecule, $shift->subject_specie, $shift->output_folder."/".$shift->subject_specie, $shift->current_directory, $shift->blast_program_path->stringify);
				my $runId = runBlastn($param, $shift->genome_subject, $shift->blast_str, $query_seq, $tag_query_spe, $molecule,  $shift->subject_specie, $shift->output_folder."/".$shift->subject_specie, $shift->current_directory, $shift->parallel_running);
				push @id_running, $runId;
			}
		}
	}
	return (\@id_running, $molecules, $query_species, $families, $files_relation);
}

sub searchHomologyBlast {
	my ($shift, $molecules, $query_species, $families, $files_relation, $Zscore, $minBitscore) = @_;
	my @id_running; # Save the id of running infernal processes
	my $result_blast_experiment = MiRNAture::_BlastSearch->new(
		blast_str => $shift->blast_str, 
		output_folder => $shift->output_folder,
		query_folder => $shift->query_folder,
		genome_subject => $shift->genome_subject,
		subject_specie => $shift->subject_specie,
		data_folder => $shift->data_folder,
		current_directory => $shift->current_directory,
		path_covariance => $shift->path_covariance,
		length_CM => $shift->length_CM,
		names_CM => $shift->names_CM, 
		bitscores_CM => $shift->bitscores_CM,
		parallel_running => $shift->parallel_running,
		blast_program_path => $shift->blast_program_path,
		cmsearch_program_path => $shift->cmsearch_program_path,	
		makeblast_program_path => $shift->makeblast_program_path,
		models_list => $shift->models_list
	);
	foreach my $molecule (@$molecules){
		foreach my $queryS (@$query_species){
			print_process("Cleaning $molecule homologs from $queryS into ", $shift->subject_specie);
			my $tag_query_spe = create_tag_specie($queryS);
			my $query_seq;
			if (exists $$files_relation{"${molecule}_${queryS}"}){
				my $file_name = $$files_relation{"${molecule}_${queryS}"};
				$query_seq = $shift->query_folder."/".$file_name;
			}
			my $queryDB = makeDatabaseQueryLen("$query_seq.len");
			cleanBlastnOut($shift->genome_subject, $shift->blast_str, $query_seq, $tag_query_spe, $molecule, $shift->subject_specie, $shift->output_folder."/".$shift->subject_specie,$queryDB);
			# Specie and molecule-specific locations
			my $output_merged = $shift->output_folder."/".$shift->subject_specie."/".$shift->subject_specie."_".$shift->blast_str.".".$molecule.".".$tag_query_spe.".tab.db.location";
			my $input_path_db_file = $shift->output_folder."/".$shift->subject_specie."/".$shift->subject_specie."_".$shift->blast_str.".".$molecule.".".$tag_query_spe.".tab.db";
			my $output_special = $shift->output_folder."/".$shift->subject_specie."/".$shift->subject_specie."_".$shift->blast_str.".".$molecule.".".$tag_query_spe.".tab.db.special";
			# Here join all to create location file: *.db -> *db.location
			$result_blast_experiment->resolveBlastMergings($input_path_db_file, $output_merged, $output_special);
		}
		#Concatenate all location files with the same str and different query
		concatenate_locations($shift, $molecule);
		my $output_merged_concatenated = $shift->output_folder."/".$shift->subject_specie."/".$shift->subject_specie."_".$shift->blast_str.".".$molecule.".tab.db.location";
		my $output_blocks = $shift->output_folder."/".$shift->subject_specie."/".$shift->subject_specie."_".$shift->blast_str.".".$molecule.".tab.db.location.blocks.coord";
		$result_blast_experiment->generate_blocks($output_merged_concatenated); #Based on location concatenated, generate blocks
		getSequencesFasta($result_blast_experiment->subject_specie, $result_blast_experiment->genome_subject, "NA", $output_blocks, "6", $result_blast_experiment->length_CM, $result_blast_experiment->names_CM); #Header mode == 6 BLAST blocks
		create_folders($result_blast_experiment->output_folder."/".$result_blast_experiment->subject_specie, "Infernal");#Create folder specific to specie
		create_folders($result_blast_experiment->output_folder."/".$result_blast_experiment->subject_specie."/Infernal","Final");#Create folder to save the classified results
		my $final_fasta_file = $output_blocks.".fasta";
		if (-z $final_fasta_file || !-e $final_fasta_file){
			next;
		} else {
			my $infernal_out_path = $shift->output_folder."/".$shift->subject_specie."/Infernal";
			my $list_file = $result_blast_experiment->models_list->stringify;
			my $all_cm_list = get_list_cms($molecule, $families, $list_file); #list of cms
			write_cmsearch_specific_sequence_group($shift->current_directory, \@$all_cm_list, $shift->path_covariance, $infernal_out_path, $shift->subject_specie, $shift->blast_str, $molecule, $shift->cmsearch_program_path, $Zscore);
			my $runIdCM = runcmSearch($shift->blast_str, $molecule, $shift->subject_specie, $shift->current_directory,$shift->parallel_running);
			push @id_running, $runIdCM;
		}
	}
	wait_processes($shift, \@id_running, $shift->parallel_running); #Wait until complete all processes from Str
	# Iterate over infernal results, evaluate and merge
	foreach my $molecule (@$molecules){
		#foreach my $queryS (@$query_species){
		#    my $tag_query_spe = create_tag_specie($queryS);
		my $infernal_out_path = $shift->output_folder."/".$shift->subject_specie."/Infernal";
		my @result_cmsearch = check_folder_files($infernal_out_path, $shift->subject_specie."\_".$shift->blast_str."\\.$molecule\\.\.*\\.tab"); #Dive_8.miRNA.RF02024.tab
		foreach my $file_out_cmsearch (@result_cmsearch){
			classify_2rd_align_results($result_blast_experiment->subject_specie, "NA", $infernal_out_path, $file_out_cmsearch,"BLAST", $molecule, $result_blast_experiment->bitscores_CM, $result_blast_experiment->length_CM, $result_blast_experiment->names_CM, $minBitscore); #Obtain true candidates
		}
		print_process("Structural evaluation of $molecule complete");
	}
	#Here, concatenate by Str
	if ($shift->blast_str =~ /^\d+$/){
		my $infernal_out_path_final = $shift->output_folder."/".$shift->subject_specie."/Infernal/Final";
		define_final_CMs($infernal_out_path_final, $shift->subject_specie, $shift->blast_str, $shift->length_CM);	
	} else {
		print_error("Missing a correct strategy");
	}
	return;
}

sub concatenate_locations {
	my $shift = shift;
	my $molecule = shift;
	my $path = $shift->output_folder."/".$shift->subject_specie."/".$shift->subject_specie."_".$shift->blast_str.".".$molecule."."."*.tab.db.location";
	my $out = $shift->output_folder."/".$shift->subject_specie."/".$shift->subject_specie."_".$shift->blast_str.".".$molecule.".tab.db.location";
	system "cat $path > $out";
	return;
}

sub join_all_blast_str {
	my $shift = shift;
	my @all_converted_cmsearch = check_folder_files($shift->output_folder."/".$shift->subject_specie."/Infernal/Final", "truetable\\.clean"); #Files generated with define_final_CMs on all strategies
	concatenate_true_cand($shift->subject_specie, $shift->output_folder."/".$shift->subject_specie."/Infernal/Final", \@all_converted_cmsearch, $shift->blast_str);
	resolve_mergings($shift->subject_specie, $shift->output_folder."/".$shift->subject_specie."/Infernal/Final", "4", $shift->blast_str);
	return;
}

sub blast_strategies {
	my $number = shift;
	my $str;
	if ($number =~ /\d+/){
		if ($number == 6){ #Default blastn
			$str = "-outfmt 6 -out";
		} elsif ($number == 1){
			$str = "-dust no -soft_masking false -evalue 0.01 -reward 5 -penalty -4 -gapopen 10 -gapextend 6 -word_size 7 -outfmt 6 -out";
		} elsif ($number == 2){
			$str = "-dust no -soft_masking false -evalue 0.01 -reward 4 -penalty -5 -gapopen 3 -gapextend 5 -word_size 7 -outfmt 6 -out";
		} elsif ($number == 3){
			$str = "-dust no -soft_masking false -evalue 0.01 -reward 5 -penalty -4 -gapopen 25 -gapextend 10 -word_size 7 -outfmt 6 -out";
		} elsif ($number == 4){
			$str = "-dust no -soft_masking false -evalue 0.01 -reward 4 -penalty -5 -gapopen 12 -gapextend 8 -word_size 7 -outfmt 6 -out";
		} elsif ($number == 9){
			$str = "-dust no -soft_masking false -evalue 0.01 -reward 4 -penalty -5 -gapopen 3 -gapextend 5 -word_size 7 -outfmt 6 -out";
		} elsif ($number == 10){
			$str = "-dust no -soft_masking false -evalue 0.01 -reward 5 -penalty -4 -gapopen 25 -gapextend 10 -word_size 7 -outfmt 6 -out";
		} elsif ($number == 7){
			$str = "-dust no -soft_masking false -evalue 0.01 -reward 4 -penalty -5 -gapopen 3 -gapextend 5 -word_size 7 -outfmt 6 -out";
		} elsif ($number == 8){
			$str = "-dust no -soft_masking false -evalue 0.01 -reward 5 -penalty -4 -gapopen 25 -gapextend 10 -word_size 7 -outfmt 6 -out";
		} elsif ($number == 5){ # Hertel parameters: doi:10.3390/life5010905
			$str = "-dust no -soft_masking false -evalue 10e-10 -outfmt 6 -out";
		} else {
			die "Invalid blastn strategy!\n";
		}
	} else { #Text 
		if ($number =~ /^ALL$|^COMPLETE$/){
			$str = "NA";
		} else {
			die "Invalid blastn strategy!\n";
		}
	}
	return $str;
}

sub index_ncRNA_families {
	my ($fileRFAM, $fileOthers) = @_;
	my %families;
	my @score_files = ("$fileRFAM", "$fileOthers");
	foreach my $file (@score_files){
		open my $FAMS, "< $file" or die "Please provide the families_ncRNAs file\n";
		while (<$FAMS>){
			chomp;
			my @all = split /\s+|\t/, $_;
			$all[-1] =~ s/^(.*RNA)(s)$/$1/g;
			if ($all[-1] eq "misc_RNA" || $all[-1] eq "antisense" || $all[-1] eq "antitoxin" || $all[-1] eq "Cis-reg" || $all[-1] eq "Cis-reg_frameshift_element" || $all[-1] eq "Cis-reg_frameshift_element" || $all[-1] eq "Cis-reg_IRES" || $all[-1] eq "Cis-reg_leader" || $all[-1] eq "Cis-reg_riboswitch" || $all[-1] eq "Cis-reg_thermoregulator" || $all[-1] eq "CRISPR" || $all[-1] eq "Intron" || $all[-1] eq "ribozyme" || $all[-1] eq "sRNA"){
				push @{ $families{"others_ncRNA"} }, $all[0];
			} else { #Here the searches are restricted to miRNAs
				push @{ $families{"miRNA"}}, $all[0];
			}
		}
		close $FAMS;
	}
	return \%families;
}

sub infer_data_queries {
	my $blastqueriesfolder = shift;
	open my $METADATA, "< $blastqueriesfolder/queries_description.txt" or die "The Metadata file describing the blast queries is missing\n";	
	my (@molecules, @species);
	my (%mol, %specie, %files);
	while (<$METADATA>){
		#tRNA.fasta	tRNA	Ciona intestinalis
		chomp;
		#print "$_\n";
		my @all = split /\t|\s+/, $_;
		$mol{$all[1]} = 1;
		$specie{"$all[-2] $all[-1]"} = 1;
		$all[0] =~ s/(.*)(\.fa$|\.fasta$|\.fna$)/$1.new.fasta/g;
		$files{"$all[1]_$all[-2] $all[-1]"} = $all[0];	
	}
	foreach my $molecule (sort keys %mol){
		push @molecules, $molecule;
	}
	foreach my $spe (sort keys %specie){
		push @species, $spe;
	}
	return (\@molecules, \@species, \%files);
}

sub makeDatabaseQueryLen {
	my $len_file = shift;
	my %lenquery;
	open my $IN, "< $len_file" or die;
	while (<$IN>){
		chomp;
		my @all = split /\t|\s+/, $_;
		$lenquery{$all[0]} = $all[1];
	}
	return \%lenquery;
}

sub writeBlastn {
	my ($parameters, $genome, $strategy, $query_seq, $query_tag, $ncrna, $specie, $out_path_blast, $dir_now, $blast_path) = @_;	
	$genome =~ s/\"//g;
	open my $OUTTEMP, "> .blastTemp/${specie}_$strategy.$ncrna.$query_tag.sh";
	my $out_file_blast = "${specie}_$strategy.$ncrna.$query_tag.tab";
	print $OUTTEMP "\#!\/bin\/bash\n";
	print $OUTTEMP "$blast_path -db $genome -query $query_seq -num_threads 4 $parameters $out_path_blast/$out_file_blast\n";
	close $OUTTEMP;
	return;
}

sub runBlastn {
	my ($parameters, $genome, $strategy, $query_seq, $query_tag, $ncrna, $specie, $out_path_blast, $dir_now, $parallel) = @_;	
	#"$tag_query_spe-$spe-$molecule"
	my $nameTTO = "$query_tag-$specie-$ncrna";
	my $id; 
	system "chmod 755 $dir_now/.blastTemp/${specie}_$strategy.$ncrna.$query_tag.sh";
	if ($parallel == 1){
		print_process("Going into parallel running!\n");
		$id = `sbatch --job-name=$nameTTO --nodes=1 --ntasks=1 --cpus-per-task=5 --time=96:00:00 --mem=2G --output=$dir_now/.blastTemp/LOGs/out_${specie}_$strategy.$ncrna.$query_tag.out --error=$dir_now/.blastTemp/LOGs/error_${specie}_$strategy.$ncrna.$query_tag.out $dir_now/.blastTemp/${specie}_$strategy.$ncrna.$query_tag.sh`;
		$id =~ /^Submitted batch job (\d+)/; 
		$id = $1;
		#print "$id\n";
	} else {
		system("$dir_now/.blastTemp/${specie}_$strategy.$ncrna.$query_tag.sh &");
		$id = `ps aux | grep ${specie}_$strategy.$ncrna.$query_tag.sh | grep bash`;
		$id = (split /\s+/,(split /\n/, $id)[0])[1];
		#$id = `ps aux | grep ${specie}_$strategy.$ncrna.$query_tag.sh | grep "bash \./" |awk '{print \$2}'`;
		#print "$id\n";
		$id =~ s/(\d+)(\s*)/$1/g;
	}
	####
	return $id;
}

sub wait_processes {
	my ($shift, $processes, $parallel) = @_;
	my $num = length $$processes[0];
	if ($num){
		return if $num == 0;
	}
	if ($parallel == 1){
		foreach my $reference (@$processes){
			chomp $reference;
			print_process("Waiting for the process: $reference"); 
			EVAL:
			my $state_qsub = `squeue -o "%8i %20j %4t %10u" | grep "^$reference"`; #Capture state
			if ($state_qsub && length $state_qsub > 0){
				sleep(3);
				$state_qsub = ();
				goto EVAL;
			} else {
				;
			} 
		}
		return;
	} else {
		foreach my $reference (@$processes){
			print_process("Waiting for the process: $reference"); 
			EVAL2:
			my $exists = kill 0, $reference;
			if ( $exists ){
				sleep(3); 
				$exists = ();
				goto EVAL2;
			} else {
				;
			} 
		}
		return;
	}
}

sub wait_slurm_process {
	my $reference = shift;
	print "Waiting for the process on $reference\n";
	EVAL:
	my $state_qsub = `squeue -o "%8i %20j %4t %10u" | grep $reference`; #Capture state
	my $state = 0;
	if ($state_qsub && length $state_qsub > 0){
		$state = 1;
		sleep(3); #Sleep 60 s.
		goto EVAL;
	}
	$state_qsub = ();
	return;
}

sub cleanBlastnOut {
	my ($genome, $strategy, $query_seq, $query_tag, $ncrna, $specie, $out_path_blast, $lenQueryDB) = @_;	
	#makeIndex($query_tag, $ncrna);
	open my $INFILE, "< $out_path_blast/${specie}_$strategy.$ncrna.$query_tag.tab" or die "Not possible to open the blastn input file: $out_path_blast/${specie}_$strategy.$ncrna.$query_tag.tab\n";
	open my $OUTT, "> $out_path_blast/${specie}_$strategy.$ncrna.$query_tag.tab.db" or die;
	while (<$INFILE>){
		chomp;
		next if $_ =~ m/^\#/;
		my $filter = filter_blast($_, $lenQueryDB);
		if ($filter && length($filter) > 0){
			print $OUTT "$filter";
		} else {
			next;
		}
	}
	close $INFILE; close $OUTT;
	return;
}

sub filter_blast {
	#dre127#0  scaffold15793-size8002#1  65.38#2   130#3     38#4      4#5       1#6       127#7     6186#8    6061#9    0.24#10    38.1#11
	my ($ln, $queryLen) = @_;
	my @temp = split /\s+|\t/, $ln;
	my $Lquery = abs($temp[7] - $temp[6]) + 1;
	my $Lsubject = abs($temp[9] - $temp[8]) + 1;
	my $totalAlign = (($Lsubject + 1) + $temp[5]); #length Subject + Gap_openings
	my $id = $temp[2];
	my $evalue = $temp[10];
	my ($result, $strand, $startT, $endT);	
	if($temp[8] < $temp[9]){
		$strand = 1;
		$startT = $temp[8];
		$endT = $temp[9];
	} else {
		$strand = -1;
		$startT = $temp[9];
		$endT = $temp[8];
	}
	if(($id >= 75) && ($totalAlign >= 20) && ($evalue <= 0.01)){
		#dre127#0	scaffold15793-size8002#1 	6186#2    6061#3	+#4	125#5	65.38#6	1#7       127#8	100#9
		$result = "$temp[0]\t$temp[1]\t$startT\t$endT\t$strand\t$Lsubject\t$id\t$temp[6]\t$temp[7]\t$$queryLen{$temp[0]}\n";
	} else {
		;
	}
	return $result;
}

sub get_list_cms {
	my ($molecule, $families, $list_file) = @_;
	my $final_list;
	my %db;
	open my $IN, "< $list_file" or die "The list file is not available\n"; 
	#Index list file
	while (<$IN>){
		chomp;
		$db{$_} = 1;
	}
	close $IN;
	if ($molecule eq "misc_RNA" || $molecule eq "antisense" || $molecule eq "antitoxin" || $molecule eq "Cis-reg" || $molecule eq "Cis-reg_frameshift_element" || $molecule eq "Cis-reg_frameshift_element" || $molecule eq "Cis-reg_IRES" || $molecule eq "Cis-reg_leader" || $molecule eq "Cis-reg_riboswitch" || $molecule eq "Cis-reg_thermoregulator" || $molecule eq "CRISPR" || $molecule eq "Intron" || $molecule eq "lncRNA" || $molecule eq "ribozyme" || $molecule eq "sRNA"){
		$final_list = $$families{"others_ncRNA"};
	} elsif ($molecule eq "miRNA.hairpin"){
		$final_list = $$families{"miRNA"};
	} else {
		$final_list = $$families{$molecule};
	}
	my $subset;
	foreach my $ln (@$final_list){
		if (exists $db{$ln}){ #Test if Family is included on the current model list
			push @{$subset}, $ln;
		}
	}
	return $subset;
}

sub write_cmsearch_specific_sequence_group {
	my ($dir_now, $cm_group, $cm_models_path, $out_path_infernal, $genome_tag, $str, $molecule_query, $cmsearch_path, $Zvalue) = @_;
	create_folders("$dir_now",".infernalTemp"); #Create folder specific to specie
	create_folders("$dir_now/.infernalTemp","LOGs"); #Create folder specific to specie
	existenceProgram($cmsearch_path);
	#Dive_9.snRNA.pema.tab.db.location.fasta
	open my $OUTTEMPC, "> .infernalTemp/${genome_tag}_$str.$molecule_query.sh";
	print $OUTTEMPC "\#!\/bin\/bash\n";
	# lacht_9.miRNA.ciro.tab.db.location.fasta
	my $query_file = "$out_path_infernal/../${genome_tag}_$str.$molecule_query.tab.db.location.blocks.coord.fasta";
	my $query_file_modified = $query_file;
	$query_file_modified =~ s/(\/.*\/|\.\.\/|.*\/)(.*)(\.tab\.db\.location\.blocks\.coord\.fasta)/$2/g;	
	#Iterate over the rescued cm list
	$Zvalue = $Zvalue/2; #The search is performed only in one strand which is defined by the previous blast search
	foreach my $cm_name (@$cm_group){
		foreach my $cm_models_path_specific (@$cm_models_path){ #Search on available CM folders
			my $cm = "$cm_models_path_specific/${cm_name}.cm"; #Here modify with the correct miRNA name
			if (!-e $cm || -z $cm ){
				next;
			} else {
				my $param = "--cpu 8 --notrunc -Z $Zvalue --noali --nohmmonly --toponly --tblout $out_path_infernal/${query_file_modified}.${cm_name}.tab $cm $query_file";
				print $OUTTEMPC "$cmsearch_path $param\n";
			}
		}
	}
	close $OUTTEMPC;
	return;
}

sub runcmSearch {
	#infernalTemp/${genome_tag}_$str.$molecule_query.$query_name_spe.sh";
	my ($strategy, $ncrna, $specie, $dir_now, $parallel) = @_;	
	my $name = "${ncrna}_${specie}";
	my $nameTTO = "$specie-$ncrna";
	my $id;
	system "chmod 755 $dir_now/.infernalTemp/${specie}_$strategy.$ncrna.sh";
	if ($parallel == 1){
		$id = `sbatch --job-name=$nameTTO --nodes=1 --ntasks=1 --cpus-per-task=4 --time=3-01:00:00 --mem=2G --output=$dir_now/.infernalTemp/LOGs/out_${specie}_$strategy.$ncrna.out --error=$dir_now/.infernalTemp/LOGs/error_${specie}_$strategy.$ncrna.out $dir_now/.infernalTemp/${specie}_$strategy.$ncrna.sh`;
		$id =~ /^Submitted batch job (\d+)/; 
		$id = $1;
	} else {
		system("$dir_now/.infernalTemp/${specie}_$strategy.$ncrna.sh 1> /dev/null &");
		$id = `ps aux | grep ${specie}_$strategy.$ncrna.sh | grep bash`;
		$id = (split /\s+/,(split /\n/, $id)[0])[1];
		#$id = `ps aux | grep ${specie}_$strategy.$ncrna.sh | grep "bash \./" |awk '{print \$2}'`;
		#print "$id\n";
		$id =~ s/(\d+)(\s*)/$1/g;
	}
	return $id;
}

sub clean_empty {
	my $shift = shift;
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_specie); #Blast/Specie
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_specie."/Infernal"); #Blast/Specie/Infernal
	delete_empty_files($shift->output_folder->stringify."/".$shift->subject_specie."/Infernal/Final"); #Blast/Specie/Infernal/Final
	return;
}

no Moose;
__PACKAGE__->meta->make_immutable;
