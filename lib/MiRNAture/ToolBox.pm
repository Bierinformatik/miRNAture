package MiRNAture::ToolBox;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(evaluate_input_flags get_basic_files test_basic_file openF create_folders is_folder_empty copy_files check_folder_files getSequencesFasta extendBlastnCoordinates get_header_name generate_key check_if_exists getSpecieName make_blast_database existenceProgram classify_2rd_align_results infer_data_from_cm infer_list_from_cm cmsearch print_error print_result print_process read_config_file calculate_Z_value calculate_minimum_bitscore getSequencesFasta_final infer_name_database_cm);

use Moose::Role; 
use File::Copy; 
use File::Find;
BEGIN {
	local $SIG{__WARN__} = sub {};
	require Bio::DB::Fasta;
	Bio::DB::->import();
}
use Bio::Seq;
use Bio::SeqIO;
use RNA;
use MiRNAture::ValidationCM;
use Data::Dumper;
use Term::ANSIColor qw(:constants); 
use YAML::Tiny;

with 'MiRNAture::Evaluate'; #Load tool subroutines

=head1 evaluate_input_flags 
    Title: evaluate_input_flags 
    Usage: evaluate_input_flags(<PATH_CURRENT_DIR>);
    Function: Evaluate each flag, based on defined constrains
    Returns: Continues the program unless found a exception.
	 At the same time, creates the defined Output folder.
=cut 

sub evaluate_input_flags {
	#cmlist file, mode run, out folder, parallel mode.
	my ($nameC, $mode, $work_folder, $subject_specie, $parallel, $repetition_threshold, $new_models) = @_;
	# Evaluate the existence of list of CM models
	my $name = $nameC;
	$name =~ s/(.*\/Data\/|\.\/Data\/|\/.*\/|\.\/)(.*)/$2/g;
	if (!-z $nameC && -e $nameC){ #This is a file, must exists and must be non-zero
		print_process("The $name covariance list will be used");
	} else {
		print_process("The default covariance list will be used");
	}
	#Check the selected mode
	if ($mode =~ /blast|rfam|hmm|mirbase|user|final/){ # This is the valid modes
		;
	} else {
		print_error("The mode $mode is not valid!");
	}
	# Check user models folder and the definition of the correct OTHER_CM mode
	if ($mode =~ /user/){
		if ($new_models eq "NA"){
			print_error("The flag -usrM <USER_MODELS_PATH> is missing. It is required to run the OTHER_CM mode.");
		}
	}
	if ($parallel == 1 || $parallel == 0){
		;	
	} else {
		print_error("The mode $parallel mode is not valid!");
	}
	# Check if work folder exists, if not create it.
	if (-d $work_folder){ #This is the outfolder, test if exists or create
		;
	} else {
		create_folders($work_folder, "");
	}
	if (length $subject_specie > 0){
		;
	} else {
		print_error("The subject specie wasn't defined, please be specific with the tag name of your subject genome");
	}
	if ($repetition_threshold =~ /default/){
		my @all = split /,/, $repetition_threshold; 
		if ($all[0] eq "default" && $all[1] == 200 && $all[2] == 100){
			return $repetition_threshold;
		} else {
			print_error("structure of -r is wrong, to modify the parameters activate the 'relax' flag.");
		}
	} elsif ($repetition_threshold =~ /relax/){
		my @all = split /\,/, $repetition_threshold;
		if ($all[0] =~ m/relax/ && $all[1] >= 0 && $all[2] >= 0){
			return $repetition_threshold;
		} else {
			print_error("structure of -r is wrong. Passed parameters are inconsistent.");
		}        
	} elsif (length $repetition_threshold == 0){
		$repetition_threshold = "default,200,100";
		return $repetition_threshold;
	} else {
		print_error("Structure of -rep flag is malformed, did not identify the 'default' or 'relax' modes");
	}
}

=head1 get_basic_files
    Title: get_basic_files
    Usage: get_basic_files(<PATH_CURRENT_DIR>);
    Function: Create basic miRNAture basic folder configuration
	  to run the prediction of miRNAs.
    Returns: Created environment in the current folder:
	 Creates: Data/, Basic_files/ including the basic
	 files to run the pipeline (RFAM scores file).
=cut 

sub get_basic_files {
	my $working_folder = shift;
	my ($scoresCMmodels, $namesCM, $LENGTH_OTHER);
	my $data_folder = "$working_folder";
	my $basic_folder = "$data_folder/Basic_files";
	#create_folders($working_folder, "Data");
	create_folders($data_folder, "Basic_files");
	if (is_folder_empty($basic_folder)){
		print_error("Seems that you do not have the scores files to run miRNAture.")
		#TODO: Consider point the location to download those files
		#print_process("Seems that you do not have the basic file to start running miRNAture, let me copy it for you on:\n $data_folder");
		#copy_files("$working_folder/Default_Data/all_RFAM_scores.txt", $basic_folder);
		#copy_files("$working_folder/Default_Data/genomes.txt", $basic_folder);
	} else { #Test if files are in the right format to be processed
		if (-s "$basic_folder/all_rfam_scores.txt"){
			test_basic_file("$basic_folder/all_rfam_scores.txt", "scores");
		} 
		if (-s "$basic_folder/all_mirbase_scores.txt"){
			test_basic_file("$basic_folder/all_mirbase_scores.txt", "scores");
		} 
		#else {
		#	copy_files("$working_folder/Default_Data/all_RFAM_scores.txt", $basic_folder);
		#	test_basic_file("$basic_folder/all_RFAM_scores.txt", "scores");
		#}
        ##if (-s "$basic_folder/genomes.txt"){
        ##	test_basic_file("$basic_folder/genomes.txt", "genomes");
        ##} else {
        ##	;
			#This file is generated automatically, so it must be on the folder. If not, is an error of miRNAture.
			#print_process("The file: $data_folder/genomes.txt is missing, miRNAture will create it.");
			#copy_files("Default_Data/genomes.txt", $basic_folder);
			#test_basic_file("$basic_folder/genomes.txt", "genomes");
        ##}
	}
	#if (-s "$data_folder/miRNA_RFAM14-1_ACC.lista"){
	#	copy_files("$working_folder/Default_Data/miRNA_RFAM14-1_ACC.lista", $data_folder);
	#} 
	return;
}

=head1 test_basic_file
    Title: test_basic_file
    Usage: test_basic_file(file_to_test, type_of_file);
    Function: Test default files, score and genome, files copied from
	default folder.		
    Returns: undef value, only performs the evaluation of the file. 
=cut 

sub test_basic_file {
	my ($prefile, $mode) = @_;
	open my $IN, "< $prefile" or die "Does not possible to open the basic file\n";
	my $countLine = 1;
	while (<$IN>){
		chomp;
		next if ($_ =~ /^$|^\#/);
		my @all = split /\s+|\t/, $_;
		if ($mode eq "scores"){
			if ($all[0] =~ /^RF\d+|.*/ && $all[1] =~ /^\d+/ && $all[2] =~ /^\d+/ && $all[3] =~ /^\w+/ && $all[3] =~ /^\w+/){
				#RF00001	38.00	119	5S_RNA	rRNA
				$countLine++;
			} else {
				print_error("The $countLine line of your $prefile is not in the right format, please correct it!");
			}
		} elsif ($mode eq "genomes"){
			next if ($_ =~ /^#/);
			my $test_path = (split /\=/, $_)[-1];
			my $spe = (split /\=/, $_)[0];
			$test_path =~ s/"//g;
			if (-s $test_path && -e $test_path){
				$countLine++;
			} else {
				print_error("In the genome file $prefile for the $spe and line $countLine exists an error!");
			}
		} else {
			print_error("Wrong value in test_basic_file subrotine");
		}
	}
	return;
}

=head1 openF 
    Title: openF 
    Usage: openF(tag_time);
    Function: Generate configuration file with all
	  parameters to run miRNAture.
    Returns: Creates the configuration file miRNAture_configuration_${tag}.txt
	 in the working folder.
=cut 

sub openF {
	my $tag = shift;
	my $confile = "miRNAture_configuration_${tag}.txt";
	open my $CONFF, "> $confile" or die "cannot open the configuration file $confile: $!";
	close $CONFF;
	return;
}

=head1 read_config_file 
    Title: read_config_file 
    Usage: read_config_file(Configuration_file);
    Function: Read configuration file generated by miRNAture in YAML format.
    Returns: A YAML object with all variables inside the configuration file.
=cut 

sub read_config_file {
	my $config_file = shift;
	my $final = YAML::Tiny->read($config_file);
	return $final;
}

##sub read_genomes_paths { #File
##	my $name = shift;
##	$name =~ s/(.*\/|\/.*\/)(.*)$/$2/g;
##	open my $genomes_paths, "< Data/Basic_files/genomes.txt" or die "The file <Data/Basic_files/genomes.txt> does not exists\n"; 
##	my $targetGenomes;
##	my %genomes;
##	while (<$genomes_paths>){
##		chomp;
##		next if ($_ =~ /^#|^$/); 
##		my @splitline = split /\=/, $_;
##		$splitline[1] =~ s/"//g;
##		$genomes{$splitline[0]} = $splitline[1];
##	}
##	foreach my $keys (sort keys %genomes){
##		$targetGenomes .= "$keys ";
##	}
##	$targetGenomes =~ s/(.*)(\s+)$/2\.$1/g; #add identification
##	store_conf_file($targetGenomes, $name);
##	return %genomes;
##}

#sub read_genomes_paths {
#    my ($tag, $path) = @_;
#    my %genomes; 
#    $genomes{$tag} = $path;
#    return \%genomes;
#}

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

=head1 is_folder_empty
    Title: is_folder_empty
    Usage: is_folder_empty(test_folder);
    Function: Test if folder is empty or contains elements. 
    Returns: Returns '1' if the folder is empty. 
=cut 

sub is_folder_empty {
	my $dirname = shift;
	opendir(my $dh, $dirname); 
	return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

=head1 copy_files
    Title: copy_files
    Usage: copy_files(file_to_copy, path_to_copy);
    Function: Use copy module to copy selected files into
	described folder.
    Returns: undef value, only performs the copy. 
=cut 

sub copy_files {
	my ($file, $location) = @_;
	if (!-e $file || -z $file){
		print_error("$file file does not exists or is empty");
	} else {
		copy($file, $location) or die "Copy failed: $!";
	}
	return;
}

=head1 check_folder_files
    Title: check_folder_files
    Usage: check_folder_files(directory, prefix);
    Function: Checks the existence of specific files in designed folder.
    Returns: Array with matching files, based on defined prefix. 
=cut 

sub check_folder_files {
	my ($dir, $prefix) = @_;
    my @files;
    find( sub {
            return unless -f; # Test if is a file
            return if /^\.\_/; # Not accept hidden/MacOS-system files
            return unless /$prefix$/;
            push @files, $_; 
        }, $dir);
    my $num = scalar @files;
    if ($num == 0){
        print_error("The $dir does not contain desired files with pattern $prefix");
    }
    return @files;
}

#sub check_folder_files {
#	my ($dir, $prefix) = @_;
#	opendir(DIR, $dir);
#	my @files = grep(/$prefix$/, readdir(DIR));
#	closedir(DIR);
#	return @files;
#}

=head1 calculate_Z_value
    Title: calculate_Z_value
    Usage: calculate_Z_value(genome);
    Function: Calculate the sequence database size in millions of nucleotides
    Returns: Number of nt in the current genome. 
=cut 

sub calculate_Z_value {
	my ($genome, $mode) = @_; #mode: Genome or Region
	if (!-e $genome || -z $genome){
		print_error("Seems that your $genome does not exists or it is empty, fix that please.");
	}
	my ($value, $value2);
	existenceProgram("esl-seqstat");
	$value = `esl-seqstat --dna $genome | grep "#"`;
	$value = (split /\s+/, $value)[-1];
	$value2 = $value;
	if ($mode eq "Genome"){
		$value = (($value * 2) / 1000000);
	} elsif ($mode eq "Region"){ #The search space is restricted only to the provided input sequence and not both strands?
		$value = (($value) / 1000000);
	} else {
		$value = (($value * 2) / 1000000);
	}
	return ($value, $value2); #z-score and genome size
}

=head1 calculate_minimum_bitscore
    Title: calculate_minimum_bitscore
    Usage: calculate_minimum_bitscore(genome);
    Function: Calculate the minimum bitscore for all CMs on specific genome.  
    Based on snoStrip explanation: 'The rule of thumb, which will be applied here, states that a sequence
    is likely to be a true homolog if the bitscore is above log_2 of the genome
    size. Therein, the genome size is the double genome length, since both strands
    have to be searched'.
    Returns: Minimum bitscore to consider a true candidate. 
=cut 

sub calculate_minimum_bitscore {
	my ($genome_size) = @_; #mode: Genome or Region
	if (!$genome_size){
		print_error("No genome size available to calculate the minimum bitscore.");
	}
	my $value = log2($genome_size); #Must be in bp
	return $value;
}

=head1 
    Title: log2
    Usage: log2(genome_size(bp));
    Function: Calculate log to base 2
    Returns: log2 value 
=cut
sub log2 {
	my $n = shift;
	my $log2 = log($n)/log(2);
	return ($log2);
}

=head1 
    Title: getSequencesFasta 
    Usage: (directory, prefix);
    Function: 
    Returns: 
=cut 

sub getSequencesFasta {
	my ($specie_name, $genome, $hmm_query, $out_fasta, $mode, $len_r, $names_r, $file_in, $complete_name) = @_; #mode: 1-> fasta as miRBase, 2-> as HMM, 3-> blast, 4-> final Coord ## $file_in is optional, only with TTO 5!, 6-> blocks
	### First, index the genome
	my $extension_blastn_candidates = 10;
	my $dbCHR;
	# Look if exists index
	if (!-e $genome || -z $genome){
		$dbCHR = Bio::DB::Fasta->new($genome, -reindex=>1);
	} else {
		$dbCHR = Bio::DB::Fasta->new($genome, -reindex=>0);
	}
	#my $dbFasta = Bio::SeqIO->new(-file => $genome, '-format' => 'Fasta');
	my $molecule_name;
	my ($IN, $OUTFILE, $DBFILE);
	if ($mode == 1 || $mode == 2){
		open $DBFILE, ">> $out_fasta/${specie_name}.db" or die; #Names DB
		open $IN, "< $out_fasta/${specie_name}.${hmm_query}.tab.true.table" or die "I need the tabular output from nhmmer\n"; #True table nhmmer 
	} elsif ($mode == 3){
		open $DBFILE, ">> $out_fasta.db" or die; #Names DB
		open $IN, "< $out_fasta" or die "I need the location file, output from blastn\n"; #Dive_9.snRNA.anca.tab.db.location
		$molecule_name = $out_fasta;
		my @tmp_m = split /\./, $molecule_name;
		$molecule_name = $tmp_m[1];
	} elsif ($mode == 4){
		open $DBFILE, ">> $out_fasta/all_RFAM_${specie_name}.truetable.joined.table.db" or die; #Names DB
		open $IN, "< $out_fasta/all_RFAM_${specie_name}.truetable.joined.table" or die "I need the tabular output\n"; #Final Tab
	} elsif ($mode == 5){
		open $DBFILE, ">> ${file_in}.db" or die;
		open ($IN, "<", $file_in) or die; 
	} elsif ($mode == 6){
		open $DBFILE, ">> $out_fasta.db" or die; #Names DB
		open $IN, "< $out_fasta" or die "I need the location file, output from blastn\n"; #Dive_9.snRNA.anca.tab.db.location
		$molecule_name = $out_fasta;
		my @tmp_m = split /\./, $molecule_name;
		$molecule_name = $tmp_m[1];
	}
	my $spe = lc($specie_name);
	my $specie_name_complete;
	if (!$complete_name){
		$specie_name_complete = getSpecieName($spe);
	} else {
		$specie_name_complete = $complete_name;
		$specie_name_complete =~ s/\_/ /g;
	}
	my (%tosearch,%db,%final_seq, %database);
	my ($header, $sta, $end, $lon, $seq);
	my $count = 0; #temporal for names!
	if ( !-e $out_fasta || -z $out_fasta ){
		print_result("The $out_fasta file does not contain homology candidates (from blast or hmm searches)");
		return;
	}
	while (<$IN>){
		chomp;
		my ($ids, $idsDouble);
		if ($mode == 1 || $mode == 5){
			my $strand = (split /\s+|\t/, $_)[1];
			#if ($strand =~ /\,/){
			#    $idsDouble = generate_key_double(); 
			#} else {
			$ids = generate_key();
			#}
		} elsif ($mode == 2 || $mode == 3 || $mode == 4 || $mode == 6){
			$ids = "Temporal$count"; #generate_key(); #ConsiderChange
		}
		my @tmp = split /\s+|\t/, $_;
		my ($frame, $gene_name, $gen_seq);
		if ($mode == 1 || $mode == 2 ){
			my $chr_length = $dbCHR->length("$tmp[0]");
			if (!$chr_length){
				my $backfile = "$genome.len";
				$chr_length = obtain_len_database($backfile, $tmp[0]); 
			}
			my ($exStart, $exEnd); 
			if ($tmp[3] eq "-"){ # In case new HMM without ACC name would be created
				$tmp[3] = $tmp[2];
			}
			if ($tmp[8] < $tmp[9]){ #Verify sense, nhmmer report - swapped. Chr, smaller, greater.
				($exStart, $exEnd) = extendBlastnCoordinates($tmp[8], $tmp[9], $tmp[4], $tmp[5], $$len_r{$tmp[3]}, $chr_length, $extension_blastn_candidates);
				$gen_seq = $dbCHR->seq( $tmp[0], $exStart, $exEnd ); #Coordinates from nhmmer output Chr,S,E
			} else {
				($exStart, $exEnd) = extendBlastnCoordinates($tmp[9], $tmp[8], $tmp[4], $tmp[5],$$len_r{$tmp[3]}, $chr_length, $extension_blastn_candidates);
				#$gen_seq = $dbCHR->seq( $tmp[0], $exEnd, $exStart ); #Coordinates from nhmmer output Chr,S,E
				$gen_seq = $dbCHR->seq( $tmp[0], $exStart, $exEnd ); #Coordinates from nhmmer output Chr,S,E
			}
			if ($mode == 2){ #HMM specific
				$frame = $tmp[11]; #Strand
			} else {
				$frame = $tmp[1];
			}
			my $other_name = $$names_r{$tmp[3]};
			### Modes:
			#>hsa-let-7a-1 MI0000060 Homo sapiens let-7a-1 stem-loop #1
			# OR
			#>dvexscf69170.-.49.191.RF00001.3.117.119 #2
			##
			$tmp[8] = $exStart;
			$tmp[9] = $exEnd;
			$gene_name = get_header_name(\@tmp, $mode, $ids, $count, $other_name, $spe, $specie_name_complete, $len_r);
		} elsif ($mode == 3){
			#dre1	scaffold2992-size13932	+	1013	1127	46	161	164
			my $chr_length = $dbCHR->length("$tmp[1]");
			if (!$chr_length){
				my $backfile = "$genome.len";
				$chr_length = obtain_len_database($backfile, $tmp[1]); 
			}
			my ($exStart, $exEnd) = extendBlastnCoordinates($tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7], $chr_length, $extension_blastn_candidates);
			if ($tmp[4] > $tmp[3]){ #Verify sense, nhmmer report - swapped. Chr, smaller, greater.
				$gen_seq = $dbCHR->seq( $tmp[1], $exStart, $exEnd); #Coordinates from nhmmer output Chr,S,E
			} else {
				$gen_seq = $dbCHR->seq( $tmp[1], $exStart, $exEnd); #Coordinates from nhmmer output Chr,S,E
			}
			$frame = $tmp[2]; #Strand
			$tmp[3] = $exStart;
			$tmp[4] = $exEnd;
			$gene_name = get_header_name(\@tmp, $mode, $ids, $count, "NA", $spe, $specie_name_complete, $len_r);
		} elsif ($mode == 4 || $mode == 5){ #Final Tables
            if ($mode == 5){
                my $extension_homology = 15; # Extension to the homology candidates
                my $chr_length = $dbCHR->length("$tmp[0]");
                ($tmp[5], $tmp[6]) = extendBlastnCoordinates($tmp[5], $tmp[6], 1, 0, 0, $chr_length, $extension_homology);
            }
			if ($tmp[6] > $tmp[5]){ #Verify sense, nhmmer report - swapped. Chr, smaller, greater.
				$gen_seq = $dbCHR->seq( $tmp[0], $tmp[5], $tmp[6]); #Coordinates from nhmmer output Chr,S,E
			} else {
				$gen_seq = $dbCHR->seq( $tmp[0], $tmp[6], $tmp[5]); #Coordinates from nhmmer output Chr,S,E
			}
			$frame = $tmp[1]; #Strand
			my $other_name = $$names_r{$tmp[10]};
			push @tmp, $other_name;
			$gene_name = get_header_name(\@tmp, $mode, $ids, $count, "NA", $spe, $specie_name_complete, $len_r);
		} elsif ($mode == 6){ # Blocks
			#dre1	scaffold2992-size13932	+	1013	1127	46	161	164
			# 20	JH126831.1	-	1045134	1045213
			#$tmp[1] = modify_chr($tmp[1], $spe);

			my $chr_length = $dbCHR->length("$tmp[1]");
			if (!$chr_length){
				my $backfile = "$genome.len";
				$chr_length = obtain_len_database($backfile, $tmp[1]); 
			}
			# Here the coordinates included the difference respect to the query 
			# on the block
			my ($exStart, $exEnd) = extendBlastnCoordinates($tmp[3], $tmp[4], 1, 0, 0, $chr_length, $extension_blastn_candidates);
			$gen_seq = $dbCHR->seq($tmp[1], $exStart, $exEnd); #Coordinates from nhmmer output Chr,S,E
			$frame = $tmp[2]; #Strand
			$tmp[3] = $exStart;
			$tmp[4] = $exEnd;
			$gene_name = get_header_name(\@tmp, $mode, $ids, $count, "NA", $spe, $specie_name_complete, $len_r);
		}
		my $output_nucleotide = Bio::Seq->new(
			-seq        => $gen_seq,
			-id         => $gene_name,
			-display_id => $gene_name,
			-alphabet   => 'dna',
		);
		my $output_nucleotideB;
		if ($frame =~ /^\-$/) {
			$output_nucleotide = $output_nucleotide->revcom();
		} elsif ($frame =~ /,/){
			$output_nucleotideB = $output_nucleotide->revcom(); #Here B is the reversed one
		}
		if ($mode == 1 || $mode == 2){
			push @{ $final_seq{$tmp[3]}{$output_nucleotide->id}}, $output_nucleotide->seq;
		} elsif ($mode == 3 || $mode == 6){
			push @{ $final_seq{$molecule_name}{$output_nucleotide->id}}, $output_nucleotide->seq;
		} elsif ($mode == 4){
			push @{ $final_seq{$tmp[2]}{$output_nucleotide->id}}, $output_nucleotide->seq;
		} elsif ($mode == 5){
			# Only at the end table of homology, evaluate if exists strand overlapping candidates. 
			if ($tmp[2] =~ /\-/){ #Some CMs has this additional thing, just report the name: MIPF0000001-100-Mammalia-F -> MIPF0000001, Only modify mirbase references with additions
				if ($tmp[2] =~ /^MIPF/){
					my @temp = split /\-/, $tmp[2];
					$tmp[2] = $temp[0];
				}
			}
			$_ = join "\t", @tmp; #Modified id of the family withouth dashes
			if ($frame =~ /^\-$|^\+$/){
				push @{ $final_seq{$tmp[2]}{$output_nucleotide->id}}, $output_nucleotide->seq;
			} elsif ($frame =~ /,/){ #Here include both to be evaluated by miRNAnchor
				my $header = $output_nucleotide->id;
				my @modification = split /\s+|\t/, $header; 
				my $idA = $modification[1]."A"; #Forward
				my $idB = $modification[1]."B"; #Reverse
				my $complement = join " ", @modification[2..$#modification];
				my $headerA = "$modification[0] $idA $complement"; #Forward
				my $headerB = "$modification[0] $idB $complement"; #Reverse
				push @{ $final_seq{$tmp[2]}{$headerA}}, $output_nucleotide->seq; #Forward 
				push @{ $final_seq{$tmp[2]}{$headerB}}, $output_nucleotideB->seq; #Reverse
			}
		}
		if ($frame =~ /,/ && $mode == 5){
			my @allA = split /\t|\s+/, $_;
			my @allB = split /\t|\s+/, $_;
			$allA[1] = "+"; #Reeplace combined strand to forward
			$allB[1] = "-"; #Reeplace combined strand to reverse
			my $lineA = join "\t", @allA;
			my $lineB = join "\t", @allB;
			$database{$ids."A"} = $lineA;
			$database{$ids."B"} = $lineB;            
		} else {
			$database{$ids} = $_;
		}
		$count++;
	}
	foreach my $acc (sort keys %final_seq){
		if ($mode == 1 || $mode == 2){
			open $OUTFILE, ">> $out_fasta/${specie_name}.${acc}.tab.true.table.fasta" or die;
		} elsif ($mode == 3 || $mode == 6){
			open $OUTFILE, ">> $out_fasta.fasta" or die;
		} elsif ($mode == 4 || $mode == 5){
			open $OUTFILE, ">> $out_fasta/all_RFAM_${specie_name}.${acc}.fasta" or die;
		}
		foreach my $ids (sort keys %{$final_seq{$acc} }){
			my $idsm = $ids;
			#$idsm =~ s/\#/ /g;
			print $OUTFILE ">$idsm\n";
			my $all = $final_seq{$acc}{$ids}; #Here is a problem
			foreach my $sq (@$all){
				if (!$sq){
					print_error("The sequence for $acc and $ids is empty");
				}
				$sq = uc($sq);
				$sq =~ s/T/U/g;
				$sq =~ s/(.{60})/$1\n/g;
				print $OUTFILE "$sq\n";
			}
		}
		close $OUTFILE;
	}
	foreach my $reg (keys %database){
		print $DBFILE "H$reg\t$database{$reg}\n";
	}
	close $IN; close $DBFILE;
    if ($mode == 1 || $mode == 5){
        system("rm .used_ids.txt");
    }
	return;
}

##############

sub getSequencesFastaSubGenome {
	my ($specie_name, $genome, $out_fasta, $input_table) = @_;  #mode:validatedStr, validatedNoStr, discarded
	my %final_seq;
	### First, index the genome
	my $dbCHR;
	# Look if exists index
	if (!-e "$genome" || -z "$genome"){
		$dbCHR = Bio::DB::Fasta->new($genome, -reindex=>1);
	} else {
		$dbCHR = Bio::DB::Fasta->new($genome, -reindex=>0);
	}
	##
    # Here read the .db file which already included +-15nt:  <30-08-21, cavelandiah> #
	open my $IN, "< $input_table.db" or die;
	my ($gen_seq);
    my $extension_genome = 250; # Extension to the homology candidates
	while (<$IN>){
		chomp;
		my $ids;
        #H1625790349	S104	+	MIPF0000079	1	88	271598	271735	33.8	0.00099	no	mir-145	88	0	Infernal	Infernal	mir-145	miRNA	RF00675
		my @tmp = split /\s+|\t/, $_;
        my $chr_length = $dbCHR->length("$tmp[1]");
        ($tmp[6], $tmp[7]) = extendBlastnCoordinates($tmp[6], $tmp[7], 1, 0, 0, $chr_length, $extension_genome);
		if ($tmp[6] < $tmp[7]){ 
			$gen_seq = $dbCHR->seq( $tmp[1], $tmp[6], $tmp[7]); 
        } else {
			$gen_seq = $dbCHR->seq( $tmp[1], $tmp[7], $tmp[6]); 
        }
		my $frame = $tmp[2]; #Strand
		my $gene_name = "$tmp[0] $specie_name $tmp[11] $tmp[1]#$tmp[6]#$tmp[7]"; #>H1595458263 Latimeria chalumnae mir-26 stem-loop 
		my $output_nucleotide2 = Bio::Seq->new(
			-seq        => $gen_seq,
			-id         => $gene_name,
			#-display_id => $gene_name,
			-alphabet   => 'dna'
		);
        # In case was detected in forward strand, just report as a genome in 5'->3'
        # because MIRfix will search in both strands
        #if ($frame eq "-") {
        #	$output_nucleotide2 = $output_nucleotide2->revcom();
        #}
		push @{ $final_seq{$output_nucleotide2->id}}, $output_nucleotide2->seq;
	}
	close $IN;
	my $outfileF = "$out_fasta/${specie_name}_subgenome.fasta";
	open my $OUTFILE, ">> $outfileF";
	foreach my $ids (sort keys %final_seq){
		my $idsm = $ids;
		print $OUTFILE ">$idsm\n";
		my $all = $final_seq{$ids};
		foreach my $sq (@$all){
			if (!$sq){
				print_error("The sequence for $ids is empty");
			}
			$sq = uc($sq);
			$sq =~ s/(.{60})/$1\n/g;
			print $OUTFILE "$sq\n";
		}
	}
	close $OUTFILE;
	return $outfileF;
}

#############

sub getSequencesFasta_final {
	my ($specie_name, $genome, $out_fasta, $input_table, $mode) = @_;  #mode:validatedStr, validatedNoStr, discarded
	my %final_seq;
	### First, index the genome
	my $dbCHR;
	# Look if exists index
	if (!-e "$genome" || -z "$genome"){
		$dbCHR = Bio::DB::Fasta->new($genome, -reindex=>1);
	} else {
		$dbCHR = Bio::DB::Fasta->new($genome, -reindex=>0);
	}
	##
	open my $IN, "< $input_table" or die;
	my ($gen_seq);
	while (<$IN>){
		chomp;
		my $ids;
		#H1595452937.0	JH127654.1	+.2	RF00929.3	1.4	77.5	486168.6	486267.7	33.8.8	0.0034.9	no.10	mir-574.1	77.2	9.3	Blast.4	ANCA.5	anca148294,anca148390.6	miRNA.7	2.8	18.9	-8.2	100.1
		my @tmp = split /\s+|\t/, $_;
		if ($tmp[6] < $tmp[7]){ 
			$gen_seq = $dbCHR->seq( $tmp[1], $tmp[6], $tmp[7]); 
        } else {
			$gen_seq = $dbCHR->seq( $tmp[1], $tmp[7], $tmp[6]); 
        }
		my $frame = $tmp[2]; #Strand
		my $gene_name = "$tmp[0] $specie_name $tmp[11] stem-loop"; #>H1595458263 Latimeria chalumnae mir-26 stem-loop 
		my $output_nucleotide2 = Bio::Seq->new(
			-seq        => $gen_seq,
			-id         => $gene_name,
			#-display_id => $gene_name,
			-alphabet   => 'dna'
		);
		if ($frame eq "-") {
			$output_nucleotide2 = $output_nucleotide2->revcom();
		}
		push @{ $final_seq{$output_nucleotide2->id}}, $output_nucleotide2->seq;
	}
	close $IN;
	my $outfileF = "$out_fasta/${specie_name}_miRNAs_$mode.fasta";
	open my $OUTFILE, "> $outfileF";
	foreach my $ids (sort keys %final_seq){
		my $idsm = $ids;
		print $OUTFILE ">$idsm\n";
		my $all = $final_seq{$ids};
		foreach my $sq (@$all){
			if (!$sq){
				print_error("The sequence for $ids is empty");
			}
			$sq = uc($sq);
			$sq =~ s/T/U/g;
            $sq =~ s/(.{60})/$1\n/g;
			print $OUTFILE "$sq\n";
		}
	}
	close $OUTFILE;
	return $outfileF;
}

sub perform_measure_MFE {
	my $file = shift;
	my %sequences;
	if (!-e $file || -z $file){
		return "NA";
	} else { #Perform evaluation
		open my $IN, "< $file";
		my $head;
		while (<$IN>){
			chomp;
			if ($_ =~ /^>/){
				$head = $_;
			} else {
				$sequences{$head} .= $_;
			}
		}
		close $IN;
		open my $OUT, "> $file.mfe";
		foreach my $head (sort keys %sequences){
			my $sequence = $sequences{$head};
			my $len = length $sequence;
			my ($ss, $mfe) = RNA::fold($sequence);
			my $code = (split /\s+|\t/, $head)[0];
			$code =~ s/(\>)(.*)/$2/g;
			print $OUT "$code\t$mfe\t$len\n";
		}
		close $OUT;
		return "$file.mfe";
	}
}

sub cross_table_final {
	my ($mfe, $allData, $output_file) = @_;
	open my $MFE, "< $mfe";
	open my $TAB, "< $allData";
	my (%energy, %data);
	while (<$MFE>){
		chomp;
		my @all = split /\s+|\t/, $_;
		$energy{$all[0]} = "$all[1] $all[2]";
	} 
	while (<$TAB>){
		chomp;
		my @all = split /\s+|\t/, $_;
		$data{$all[0]} = $_;
	}
	open my $FINAL, "> $output_file" or die;
	foreach my $code (sort keys %energy){
		if (exists $data{$code}){
			print $FINAL "$data{$code}\t$energy{$code}\n";
		} else {
			print_error("Something is broken on the MFE assignment");
		}
	}
	close $FINAL;
	return;
}

sub obtain_len_database {
	my ($backfile, $namechr) = @_;
	open my $IN, "< $backfile" or die;
	while (<$IN>){
		chomp;
		my $id = (split /\s+|\t/, $_)[0];
		my $len = (split /\s+|\t/, $_)[1];
		if ($id =~ /^$namechr$/){
			return $len;
		}
	}
	print_error("Something is wrong with the database references of $namechr sequence, fatal error");
}

sub extendBlastnCoordinates {
	my ($start_s, $end_s, $start_q, $end_q, $len_q, $len_chr, $extension) = @_;
	my ($ex_start, $ex_end);
	my $diffS = abs($start_q - 1); #|--.
	my $diffE = abs($len_q - $end_q); #.---|
	$ex_start = (($start_s - $diffS) - $extension);
	$ex_end = ($end_s + $diffE + $extension);
	if ($ex_start < 1){
		$ex_start = 1;
	}
	if ($ex_end > $len_chr){
		$ex_end = $len_chr;
	}
	return ($ex_start, $ex_end);
}

sub get_header_name {
	my ($info, $mode, $id, $count, $other_name, $spe, $specie_nm_all, $len) = @_;
	my $name;
	if ($mode == 1){ #Header as miRBase
		$name = "$spe-$other_name-${count}#A$id#$specie_nm_all#$$info[3]#ncRNA"; #New ID is 'A'
	} elsif ($mode == 2){ #Header as HMM requirements
		#dvexscf69170.-.49.191.RF00001.3.117.119
		if ($$info[8] <= $$info[9]){
			$name = "$$info[0].$$info[11].$$info[8].$$info[9].$$info[3].$$info[4].$$info[5].$$len{$$info[3]}";
		} else {
			$name = "$$info[0].$$info[11].$$info[9].$$info[8].$$info[3].$$info[4].$$info[5].$$len{$$info[3]}";
		}
	} elsif ($mode == 3){
		#>dvexscf46521.+.7843.8520.sce1.46.437.483
		$name = "$$info[1].$$info[2].$$info[3].$$info[4].$$info[0].$$info[5].$$info[6].$$info[7]";
	} elsif ($mode == 4){
		$name = "$$info[0].$$info[1].$$info[5].$$info[6].$$info[2].$$info[3].$$info[4].$$info[7]";
	} elsif ($mode == 5){
		#JH126831	-	mir-942	1	86	7	67	29.2	3.7e-05	no	mir-942	86	0 RF00997
		#ENSLACT00000020000.1	+	RF01014	1	74	50	123	99.5	1.4e-28	no	mir-1306	74 10,9	Blast	ANCA,DARE,XELA,XETR	anca148572,dare191,xela36937,xetr1950	miRNA RF01014
		$name = "$spe-$$info[-1]-${count} H$id $specie_nm_all $$info[10] stem-loop"; 
	} elsif ($mode == 6){
		# 20    JH126831.1  -   1045134 1045213
		$name = "$$info[1].$$info[2].$$info[3].$$info[4].$$info[0]";
	}
	return $name;
}

sub generate_key {
    my $folder = shift;
	FLAG:	
	my $id = time() + int rand(10000); #entropy_source->get_int(12345);#time(); #+ rand(1000);
	$id =~ s/(.*)(\..*)/$1/g;
	my $exist = check_if_exists($id);
	if ($exist == 1){
		goto FLAG;
	} elsif ($exist == 0){
		;
	}
	return $id;
}

sub generate_key_double {
	my $folder = shift;
    FLAG:	
	my $id = time() + int rand(10000); #entropy_source->get_int(12345);#time(); #+ rand(1000);
	$id =~ s/(.*)(\..*)/$1/g;
	my $exist = check_if_exists($id);
	if ($exist == 1){
		goto FLAG;
	} elsif ($exist == 0){
		;
	}
	my $id1 = $id."A";
	my $id2 = $id."B";
	my @ids;
	push @ids, $id1;
	push @ids, $id2;
	return \@ids;
}

sub check_if_exists {
	my $test = shift;
	my $OUTID;
	my $IDS;
	my %usedIds;
	if (!-f ".used_ids.txt"){
		open $OUTID, ">> .used_ids.txt";
	} else {
		open $IDS, "< .used_ids.txt";
		open $OUTID, ">> .used_ids.txt";
		while (my $ln = <$IDS>){
			chomp $ln;
			$usedIds{$ln} = 0;
		}
		close $IDS;
	}
	if (exists $usedIds{$test}){
		return "1";
	} else {
		print $OUTID "$test\n";
		close $OUTID;
		$usedIds{$test} = 0;
		return "0";
	}
}

sub getSpecieName {
	my $abb = shift;
	my $specie_name;
	if ($abb eq "bole"){
		$specie_name = "Botrylloides leachii";
	} elsif ($abb eq "bosc"){
		$specie_name = "Botryllus schlosseri";
	} elsif ($abb eq "brbe"){
		$specie_name = "Branchiostoma belcheri";
	} elsif ($abb eq "brfl"){
		$specie_name = "Branchiostoma floridae";
	} elsif ($abb eq "ciin"){
		$specie_name = "Ciona intestinalis";
	} elsif ($abb eq "ciro"){
		$specie_name = "Ciona robusta";
	} elsif ($abb eq "cisa"){
		$specie_name = "Ciona savignyi";
	} elsif ($abb eq "clob"){
		$specie_name = "Clavelina oblonga";
	} elsif ($abb eq "dare"){
		$specie_name = "Danio rerio";
	} elsif ($abb eq "dive"){
		$specie_name = "Didemnum vexillum";
	} elsif ($abb eq "haau"){
		$specie_name = "Halocynthia aurantium";
	} elsif ($abb eq "haro"){
		$specie_name = "Halocynthia roretzi";
	} elsif ($abb eq "lach"){
		$specie_name = "Latimeria chalumnae";
	} elsif ($abb eq "mata"){
		$specie_name = "Molgula oculata";
	} elsif ($abb eq "mlis"){
		$specie_name = "Molgula occidentalis";
	} elsif ($abb eq "mlta"){
		$specie_name = "Molgula occulta";
	} elsif ($abb eq "oidi"){
		$specie_name = "Oikopleura dioica";
	} elsif ($abb eq "pema"){
		$specie_name = "Petromyzon marinus";
	} elsif ($abb eq "pami"){
		$specie_name = "Patiria miniata";
	} elsif ($abb eq "pevi"){
		$specie_name = "Perophora viridis";
	} elsif ($abb eq "phfu"){
		$specie_name = "Phallusia fumigata";
	} elsif ($abb eq "phma"){
		$specie_name = "Phallusia mammillata";
	} elsif ($abb eq "sath"){
		$specie_name = "Salpa thompsoni";
	} elsif ($abb eq "stpu"){
		$specie_name = "Strongylocentrotus purpuratus";
	} elsif ($abb eq "saco"){
		$specie_name = "Saccoglossus kowalewski";
	} else {
		#print "$abb doesn´t recognized!\n";
		$specie_name = "Specie unrecognized";
	}
	return $specie_name;
}

sub make_blast_database {
	my ($genome, $makeblast_path) = @_;
	my $param = "-in $genome -dbtype nucl -out $genome";
	existenceProgram($makeblast_path);
	system "$makeblast_path $param 1> /dev/null";
	return;
}


sub existenceProgram {
	my $tool = shift;
	my $tool_path = `which $tool 2> /dev/null | tr -d '\n'`;
	unless ($tool_path){
		print_error("No $tool command is available in your system!");
	}
	return;
}


sub classify_2rd_align_results {
	my ($spe, $cm, $folder_in, $file, $mode, $molecule, $cm_scores, $len_scores, $names_cms, $minBitscore, $maxthreshold) = @_;
	my $filein;
	if ($cm ne "NA"){
		#$filein = "$folder_in/$spe.$cm.tab";
		if ($mode eq "INFERNAL" || $mode eq "Infernal" || $mode eq "OTHER_CM"){
			$filein = $file; 
		} elsif ($mode eq "HMM"){
			$filein = $file;
		}
	} else { #Here is the Blast/Infernal input data, not too much input data, then reeplace.
		#Ncov_2.snRNA.ciin.RF00438.tab
		$filein = $file; #$cm_scores; #Here we provide the file
		$filein = "$folder_in/$filein";
	}
	if (-e $filein && !-z $filein){
		if ($molecule =~ /^miRNA/){ #Based on evidence, miRNAs are evaluated with 32% of Bitscore.
			cleancmsearch($filein, $maxthreshold, 1, $cm_scores, $len_scores, $names_cms ,$minBitscore); #filein, threshold bitscore, mode GA score miRNAture
			#cleancmsearch($filein, 1, 1, $cm_scores, $len_scores, $names_cms, $minBitscore); #filein, threshold bitscore, mode GA score
		} elsif ($molecule =~ /^NA$/){ #other CMs without bitscore
			cleancmsearch($filein, $maxthreshold, 2, $cm_scores, $len_scores, $names_cms, $minBitscore);
		} else { # Other RNA families
			cleancmsearch($filein, 1, 1, $cm_scores, $len_scores, $names_cms, $minBitscore); 
		}
	}
	return;
}

sub infer_name_database_cm {
	my ($new_cm_folder) = @_;
    my @all_models;
	foreach my $path (@$new_cm_folder){
        my @all_cm_models = check_folder_files($path, "\.cm"); #All files must end at '.cm'
        foreach my $file (@all_cm_models){
            my $complete =  "$path/$file";
            push @all_models, $complete;
        }
    }
	my %database;
	my $IN;
    # Read each CMs and collect its name, build DB
	foreach my $file (@all_models){
		next if $file !~ /\.cm$/;
		open $IN, "< $file" or die;
		my ($name, $acc);
		while (<$IN>){
			chomp;
			if ($_ =~ /^NAME\s/){
				$name = (split /\s+|\t/)[1];
                $database{$name} = $file;
			} elsif ($_ =~ /^ACC\s/){
                $acc = (split /\s+|\t/)[1];
                $database{$acc} = $file;
            } else{
				next;
			}
		}
        close $IN;
	}
	return \%database;
}

sub infer_data_from_cm {
	my ($new_cm_folder, $basic_files, $name) = @_;
	my @all_cm_models = check_folder_files($new_cm_folder, "\.cm"); #All files must end at '.cm'
	my %complete_cm_information;
	my $IN;
	open my $OUT, "> $basic_files/$name" or die "The table with all CMs scores, cannot be created at $new_cm_folder\n";
	#RF00006    34.00   96  Vault   misc_RNA
	foreach my $file (@all_cm_models){
		my $complete_path = $new_cm_folder."/".$file;
		next if $file !~ /\.cm$/;
		open $IN, "< $complete_path" or die;
		my ($length, $name);
		my $family = "miRNA"; #This family does not have identity
		my $acc = "NA";
		my $score = 0;
		while (<$IN>){
			chomp;
			if ($_ =~ /^NAME\s/){
				$name = (split /\s+|\t/)[1];
			} elsif ($_ =~ /^CLEN\s/){
				$length = (split /\s+|\t/)[1];
			} elsif ($_ =~ /^ACC\s/){
				$acc = (split /\s+|\t/)[1];
			} elsif ($_ =~ /^GA\s/){
				$score = (split /\s+|\t/)[1];
			} else {
				next;
			}
		}
		if ($acc eq "NA"){
			$acc = $name;
		}
		my $new_line_cm_families = "$acc\t$score\t$length\t$name\t$family";
		print $OUT "$new_line_cm_families\n";
        close $IN;
	}
	close $OUT;
	return;
}

sub infer_list_from_cm {
	my ($cm_model_others, $data_folder) = @_;
	my @all_cm_models = check_folder_files($cm_model_others, "\.cm"); #All files must end at '.cm'
	open my $OUT, "> $data_folder/new_cm.list" or die "The table with all CMs scores, cannot be created at $cm_model_others\n";
	foreach my $file (@all_cm_models){
		$file =~ s/(.*)(\.cm)/$1/g;
		print $OUT "$file\n";
	}
	close $OUT;
	return;
}

sub cmsearch {
	#In: out, CM model, genome
	my ($nameCM,$genome,$genomeTag,$outFolder, $path_cm, $mode, $cmsearch_path, $zscore) = @_; #modes: 1 X.cm 0 other name
	existenceProgram($cmsearch_path);
	my $nameCMFinal;
	if($mode == 1){ #From other source
		$nameCMFinal = "${nameCM}.cm";
	} elsif ($mode == 0){ #From metazoa
		$nameCMFinal = "${nameCM}_metazoa.cm";
	}
	$genome =~ s/"//g;
	if (-e "$path_cm/${nameCMFinal}" && !-z "$path_cm/${nameCMFinal}"){
		my $param = "--cpu 5 -E 0.015 --notrunc -Z $zscore --nohmmonly --tblout $outFolder/${nameCM}_$genomeTag.tab -o $outFolder/${nameCM}_$genomeTag.out $path_cm/${nameCMFinal} $genome";
		system "$cmsearch_path $param 1> /dev/null";
	} else {
		;
	}
	return;
}


sub print_error {
	my $text = shift;
	local $Term::ANSIColor::AUTORESET = 1;
	print BOLD RED "[ERROR] $text\n";
	die "miRNAture failed to continue working due derected errors\nExiting...\n";
	return;
}

sub print_result {
	my $text = shift;
	local $Term::ANSIColor::AUTORESET = 1;
	print BOLD GREEN "[RESULT] $text\n";
	return;
}

sub print_process {
	my $text = shift;
	local $Term::ANSIColor::AUTORESET = 1;
	print FAINT BLUE "[PROCESS] $text\n";
	return;
}

no Moose::Role;
1;
