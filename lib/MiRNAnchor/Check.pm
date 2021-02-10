package MiRNAnchor::Check;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(process_fasta_sequences check_candidate load_database_miRNAs describe_candidate process_old_line update_coordinates change_to_dna run_blast get_blast_information test_no_empty_variables check_validity_mature_coordinates check_blast_database make_blast_database check_output_files discard_candidate check_if_valid read_fasta_results concatenate_results validate_secondary_structure_alignment get_mature_map_location get_most_frequent);

use Moose;
use MooseX::Types::Path::Class;
use Data::Dumper;
use Bio::SearchIO;
use Bio::Tools::Run::StandAloneBlast;

my %mappedCoordinates; #Mapped coordinates back to genome

with 'MiRNAnchor::Tools';

has 'source_miRNAs_fasta' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1
);

has 'final_miRNA_database' => (
	is => 'rw',
	isa => 'Path::Class::File',
	coerce => 1,
);

has 'output_folder' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1
);

has 'current_dir' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1
);

has 'accession_number_RFAM' => (
	is => 'ro',
	isa => 'Str',
	required => 1
);

has 'subject_tag' => (
	is => 'ro',
	isa => 'Str',
	required => 1
);

has 'subject_genome_path' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'subject_genome_path_original' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'accepted_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'discarded_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'fasta_results' => (
	is => 'rw',
	isa => 'HashRef[Str]',
);

has 'db_codes' => (
	is => 'rw',
	isa => 'HashRef[Str]',
);

has 'mature_description' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

sub process_fasta_sequences {
	my $shift = shift;
	my $fasta_query = (check_folder_files($shift->source_miRNAs_fasta, ".*\_".$shift->subject_tag.".".$shift->accession_number_RFAM."\.*\.fasta"))[0];
	my $fasta_candidates = read_fasta_results($shift->source_miRNAs_fasta."/".$fasta_query);
	$shift->fasta_results($fasta_candidates);
	my $accepted = $shift->accepted_file->stringify;
	my $discarded = $shift->discarded_file->stringify;
	open (my $OUT1, ">", $accepted);
	open (my $OUT2, ">", $discarded);
	my $database_codes = load_database_miRNAs($shift->final_miRNA_database->stringify);
	$shift->db_codes($database_codes);
	return;
}

sub detect_best_group {
	my ($shift, $work_path, $id, $RFAM_ACC) =  @_;
	# RF00711_1-mirstar-map.txt
	my $best;
	my $path = "$work_path/$id/Output/$RFAM_ACC.out/$RFAM_ACC-mirstar-map.txt";
	my $path0 = "$work_path/$id/Output/$RFAM_ACC.out/${RFAM_ACC}_0-mirstar-map.txt";
	my $path1 = "$work_path/$id/Output/$RFAM_ACC.out/${RFAM_ACC}_1-mirstar-map.txt";
	if (-e $path && !-z $path){  #Default case, there is no groups
		$best = 2;
	} else {
		if (-e $path0 && !-z $path0){  #Selected group was 0
			$best = 0;
		} else {
			if (-e $path1 && !-z $path1){  #Selected group was 1
				$best = 1;
			} else {
				print_error("Corrected family files failed to be predicted");
			}
		}
	}
	return $best;
}

sub check_candidate {
	my $shift = shift;
	foreach my $id (sort keys %{ $shift->fasta_results }){
		my $pass1 = check_output_files($shift, $shift->output_folder, $id, $shift->accession_number_RFAM); #Check if exist correct output file from MIRfix
		if ($pass1 == 0){
			next;
		}
		my $pass2 = check_validity_mature_coordinates($shift, $shift->output_folder, $id, $shift->accession_number_RFAM); #Check if reported valid mature coordinates
		if ($pass2 == 0){
			next;
		}
		my $new_line = describe_candidate($shift, $shift->output_folder, $id, $shift->db_codes, $shift->accession_number_RFAM); #Mapp to genome
		my $bestGroup = detect_best_group($shift, $shift->output_folder, $id, $shift->accession_number_RFAM); #Detect the best group in special cases: 0, 1 or NA.
		if ($new_line eq "NA"){
			next;
		} else {
			open (my $ACCEPTED, ">>", $shift->accepted_file->stringify) or die "File with accepted candidates couldn't be open\n";
			print $ACCEPTED "$new_line\t$bestGroup\n";
			close $ACCEPTED;
		}
	}
	return;
}

sub validate_secondary_structure_alignment {
	my ($sto, $current_dir, $positionsBest) = @_;
	my $result = system("evaluate_conserved_str.py $sto $positionsBest 2>/dev/null 1>/dev/null");
	my $evaluation_result;
	if ($result == 0){ #Success!
		$evaluation_result = `evaluate_conserved_str.py $sto $positionsBest 2>&1`;
		chomp $evaluation_result;
	} else { #Failed to read the file, didn't generated
		$evaluation_result = "NA";
	}
	return $evaluation_result;
}

sub load_database_miRNAs {
	my $file_database = shift;
	open (my $DB, "<", $file_database) or die "Database file couldn't be openned\n";
	my %in;
	while (<$DB>){
		chomp;
		my $code = (split /\s+|\t/, $_)[0];
		$in{$code} = $_;
	}
	return \%in;
}

sub describe_candidate {
	#Describe the accepted candidate, in this case report onto the table but with updated genome coordinates
	my ($shift, $work_path, $id, $info, $RFAM_ACC) =  @_;
	open my $IN, "< $work_path/$id/Output/$RFAM_ACC.out/$RFAM_ACC-FinalCoor.txt" or die "Final Coordinates file does not exists!\n";
	my $line;
	while (<$IN>){
		chomp;
		if ($_ =~ m/^H/){
			$line = $_; #Information about the mature position
		}
	}
	my $fasta_file = "$work_path/$id/Output/$RFAM_ACC.out/$RFAM_ACC-Final.fasta" or die "Fasta File does not exist!\n";	
	my $fastaCorrected_file = "$work_path/$id/Output/$RFAM_ACC.out/$RFAM_ACC-FinalCorrected.fasta";
	if (-e $fastaCorrected_file && !-z $fastaCorrected_file){
		$fasta_file = $fastaCorrected_file;
	}
	#H1608041296	9	-	MIPF0000483	1	78	71626963	71627060	37.8	0.00053	no	MIPF0000483	78	HMM	HMM	HMM MIPF0000483	miRNA	MIPF0000483
	my $old_info_miRNA = $$info{$id};
	my $targetChr = (split /\s+|\t/, $old_info_miRNA)[1]; #Get Name chr
	my $targetStart = (split /\s+|\t/, $old_info_miRNA)[6]; #Get reported homology start
	my $targetEnd = (split /\s+|\t/, $old_info_miRNA)[7]; #Get reported homology end
	#Based on the reported new fasta seq by MIRfix, map onto genome and get new coordinates
	my ($Nchr, $Nstart, $Nend, $Nstrand, $queryLen) = update_coordinates($shift, $fasta_file, $line, $old_info_miRNA, $targetChr); #Changed genome_subject to Chr name
	# Include mapped coordinates
	$mappedCoordinates{"$targetChr $Nstrand $Nstart $Nend"} = 1;
	# Update annotated line
	if ($Nchr eq "NA" || $Nstart eq "NA" || $Nend eq "NA" || $Nstrand eq "NA" || $queryLen < 20){
		discard_candidate($shift, $work_path, $id, "SHORT_SEQ_MIRFIX");		
		return "NA";
	} else {
		my $updated_line = process_old_line($old_info_miRNA, $Nchr, $Nstart, $Nend, $Nstrand);
		return $updated_line;
	}
}

sub process_old_line {
	my ($old_info_miRNA, $Nchr, $Nstart, $Nend, $Nstrand) = @_;
	my @split = split /\s|\t/, $old_info_miRNA;
	$split[1] = $Nchr;
	$split[2] = $Nstrand;
	if ($Nstart < $Nend){
		$split[6] = $Nstart;
		$split[7] = $Nend;
	} else {
		$split[6] = $Nend;
		$split[7] = $Nstart;
	}
	my $updated_line = join "\t", @split;
	return $updated_line;
}

sub infer_length_query {
	my ($fasta, $code) = @_;
	open my $IN, "< $fasta" or die; 
	open my $OUT, "> $fasta.2";
	my %select;
	my $head;
	# Load all the db and change to DNA nt
	while (<$IN>){
		chomp;
		if ($_ =~ /^>/){
			$head = $_;
		} else {
			$_ =~ tr/Uu/Tt/;
			push @{$select{$head}}, $_;
		}
	}
	# Select only the interest sequence
	my $length;
	my $tmp_seq;
	foreach my $seq (sort keys %select){
		if ($seq =~ m/$code/){
			print $OUT "$seq\n";
			my $sequence = $select{$seq};
			foreach my $ln (@$sequence){
				$tmp_seq .= "$ln";
				print $OUT "$ln\n";
			}
			$length = length $tmp_seq;
		}
	}
	close $OUT;
	close $IN;
	return ("$fasta.2", $length);
}

sub update_coordinates {
	my ($shift, $fasta_query, $matureInfoLine, $old_info_miRNA, $chrName) = @_;
	#my $genome = "$splitGenome/$chrName.fa"; #Pointer to the splitted Chromosome
	my $genome = $shift->subject_genome_path_original->stringify;
	check_blast_database($genome);
	my $sequence_code = (split /\s+|\t/, $old_info_miRNA)[0];
	my ($selected_fasta_query, $subject_len) = infer_length_query($fasta_query, $sequence_code); # Selects only the subject sequence and infer length
	my ($chr, $start, $end, $strand, $queryLen) = run_blast($genome, $selected_fasta_query, $subject_len, $chrName);
	if ($strand =~ /^-1$/){
		$strand = "-";
	} else {
		$strand = "+";
	}
	return ($chr, $start, $end, $strand, $queryLen);
}

sub change_to_dna {
	my $original = shift;
	open my $IN, "< $original";
	open my $OUT, "> $original.2";
	while (<$IN>){
		chomp;
		if ($_ =~ /^>/){
			print $OUT "$_\n";
		} else {
			$_ =~ tr/Uu/Tt/;
			print $OUT "$_\n";
		}
	}
	close $OUT;
	close $IN;
	return "$original.2";
}

sub run_blast {
	my ($genome, $query, $query_len, $chrName) = @_;
	#my $query2 = change_to_dna($query);
	my @params = (  -program => 'blastn',
		-outfile => "$query.out",
		-database => $genome,
		-F => "F", #Disable low complexity
		-W => 30);
	my $factory =  Bio::Tools::Run::StandAloneBlast->new(@params);
	my $blast_report = $factory->blastall($query);
	my $output = Bio::SearchIO->new(-format => "blast", 
		-file => "$query.out",
		-best_hit_only => 1,
		-hit_filter => sub { my $hit = shift;
			$hit->name eq "$chrName"; 
		}
	);
	my ($chr, $start, $end, $strand, $queryLen, $identity) = get_blast_information($output, $query_len, $chrName);
	system "rm $query"; #Delete temp fasta
	system "rm $query.out"; #Delete BLAST output 
	return ($chr, $start, $end, $strand, $queryLen);
}

sub get_blast_information {
	my ($output, $query_len, $chrName) = @_;
	my ($chr, $start, $end, $strand, $queryLen, $identity);
	my $indicator = 0;
	while (my $result = $output->next_result){
		if ($indicator == 0){
			while (my $hit = $result->next_hit){
				while (my $hsp = $hit->next_hsp){
					if ($indicator == 0){
						$chr = $hit->name;
						$start = $hsp->start('hit');
						$end = $hsp->end('hit');
						$strand = $hsp->strand('hit');
						$queryLen = $hsp->length('query');
						$identity = $hsp->percent_identity;
						my $nstrand;
						if ($strand == 1){
							$nstrand = "+";
						} elsif ($strand == -1){
							$nstrand = "-";
						}
						if (($identity >= 95) && ($chr eq $chrName) && (!exists $mappedCoordinates{"$chr $nstrand $start $end"})){
							$indicator = 1;
							last;
						} else {
							;
						}
					}
				}
			}
		}
	}
	$chr = test_no_empty_variables($chr);
	$start = test_no_empty_variables($start);
	$end = test_no_empty_variables($end);
	$strand = test_no_empty_variables($strand);
	$queryLen = test_no_empty_variables($queryLen);
	$identity = test_no_empty_variables($identity);
	return ($chr, $start, $end, $strand, $queryLen, $identity);
}

sub test_no_empty_variables {
	my $test_variable = shift;
	if (length $test_variable){
		return $test_variable;
	} else {
		return "NA";
	}
}

sub check_validity_mature_coordinates {
	my ($shift, $work_path, $id, $RFAM_ACC) =  @_;
	open my $IN, "< $work_path/$id/Output/$RFAM_ACC.out/$RFAM_ACC-FinalCoor.txt";
	my $Nlines = `wc -l $work_path/$id/Output/$RFAM_ACC.out/$RFAM_ACC-FinalCoor.txt | cut -d \" \" -f 1`;
	chomp $Nlines;
	my $countLines = 0;
	my $des2 = 0;
	while (<$IN>){
		chomp;
		if ($_ =~ /^H/){ #The file contains a specie candidate
			if ($_ =~ m/\s-1\s|\s0\s/){ # Removed \sNULL\s, here I allow to annotate only one sequence
				$des2 = 0;
				print_result("Seems that $id does not fit with the alignment");
				discard_candidate($shift, $work_path, $id, "MISMATCH_ON_ALIGN:-1");
			} else {
				$des2 = 1;
				last;
			}	
		} else {
			if ($_ =~ /^B|^MI/ && $countLines < $Nlines-1){
				$countLines++;
			} elsif ($_ =~ /^B|^MI/ && $countLines == $Nlines-1){ #Until the end of the file there is not ^H lines
				$des2 = 0;
				print_result("Seems that $id was not reported on the Final Coordinates file");
				discard_candidate($shift, $work_path, $id, "NOT_REPORTED_ON_COOR_FILE");
			} 
		}
	}
	close $IN;
	return $des2;
}

sub check_blast_database {
	my $file = shift;
	my $database_file = "$file.nhr"; #Check if exists nhr file, an output from makeblastdb program
	if (-e $database_file && !-z $database_file){
		return;
	} else {
		make_blast_database($file);
	}
	return;
}

sub make_blast_database {
	my $genome = shift;
	my $param = "-in $genome -dbtype nucl -out $genome";
	system "makeblastdb $param 2>/dev/null 1>/dev/null";
	return;
}

sub check_output_files {
	my ($shift, $work_path, $id, $RFAM_ACC) =  @_;
	#Check if exists RF0*-FinalCoor.txt
	my $des1 = check_if_valid("$work_path/$id/Output/$RFAM_ACC.out/$RFAM_ACC-FinalCoor.txt");	
	my $pass;
	if ($des1 == 0){ #Print to the discarted file because MIRfix failed
		print_result("The sequence $id failed the mature annotation");
		discard_candidate($shift, $work_path, $id, "NO_OUTPUT_FILE");
		$pass = 0;
	} else {
		$pass = 1;
	}
	return $pass;
}

sub discard_candidate {
	my ($shift, $path, $id, $reason) = @_;
	open (my $DISCARDED, ">>", $shift->discarded_file->stringify) or die "Couldn't open file to discard all bad candidates\n";
	print $DISCARDED $id."\t".$shift->fasta_results->{$id}."\t".$reason."\n";
	close $DISCARDED;
	return;
}

sub check_if_valid {
	my $file = shift;
	my $recall = 0; #If 0 is not valid, if 1: valid
	if (-e $file && !-z $file){
		$recall = 1;
	} else {
		$recall = 0;
	}
	return $recall;
}

sub read_fasta_results {
	my $location = shift;
	open my $FASTA, "< $location";
	my %dbHeaders;
	while (<$FASTA>){
		chomp;
		if ($_ =~ /^\>/){
			my $code = (split /\s+|\t/, $_)[1];
			$dbHeaders{$code} = $_;
		} else {
			next;
		}	
	}
	return \%dbHeaders;
}

sub generate_gff {
	my ($out, $positive, $info) = @_;
	open (my $GFF, ">>", $out) or die "GFF file fails to be created\n";
	if (-e $out && -z $out){ #Print header
		print $GFF "\#\#gff-version 3\n";
	}
	foreach my $id (sort keys %{$positive}){
		# Check if exist in positive
		print $GFF "$$positive{$id}[0]\tmiRNAture\tpre_miRNA\t$$positive{$id}[2]\t$$positive{$id}[3]\t$$positive{$id}[-1]\t$$positive{$id}[1]\t\.\tID=$id;NAME=$$positive{$id}[4];Alias=$$positive{$id}[5]\n"; 
		if (exists $$info{$id}){
			my $start5p = $$positive{$id}[2] + $$info{$id}[2]; 
			my $end5p = $$positive{$id}[2] + $$info{$id}[3];
			my $start3p = $$positive{$id}[2] + $$info{$id}[4];
			my $end3p = $$positive{$id}[2] + $$info{$id}[5];
			print $GFF "$$positive{$id}[0]\tmiRNAture\tmiRNA\t$start5p\t$end5p\t\.\t$$positive{$id}[1]\t\.\tID=${id}_5p;NAME=$$positive{$id}[4];Alias=$$positive{$id}[5];Parent=$id\n"; 
			print $GFF "$$positive{$id}[0]\tmiRNAture\tmiRNA\t$start3p\t$end3p\t\.\t$$positive{$id}[1]\t\.\tID=${id}_3p;NAME=$$positive{$id}[4];Alias=$$positive{$id}[5];Parent=$id\n";
		} else {
			print_error("The database of results on the validation step is corrupted");
		}
	}
	close $GFF;
	return;
}

sub generate_bed {
	my ($out, $positive, $info) = @_;
	open (my $BED, ">>", $out) or die "GFF file fails to be created\n";
	foreach my $id (sort keys %{$positive}){
		# Check if exist in positive
		print $BED "$$positive{$id}[0]\t$$positive{$id}[2]\t$$positive{$id}[3]\t$$positive{$id}[4]\t0\t$$positive{$id}[1]\n"; 
		if (exists $$info{$id}){
			my $start5p = $$positive{$id}[2] + $$info{$id}[2]; 
			my $end5p = $$positive{$id}[2] + $$info{$id}[3];
			my $start3p = $$positive{$id}[2] + $$info{$id}[4];
			my $end3p = $$positive{$id}[2] + $$info{$id}[5];
			print $BED "$$positive{$id}[0]\t$start5p\t$end5p\t$$positive{$id}[4]_5p\t0\t$$positive{$id}[1]\n"; 
			print $BED "$$positive{$id}[0]\t$start3p\t$end3p\t$$positive{$id}[4]_3p\t0\t$$positive{$id}[1]\n";
		} else {
			print_error("The database of results on the validation step is corrupted");
		}
	}
	close $BED;
	return;
}

sub concatenate_results {
	my $shift = shift;
	my $working_path = shift;
	my $tag = shift;
	open (my $ACCEPTED, "<", $shift->accepted_file->stringify);
	open (my $DISCARDED, "<", $shift->discarded_file->stringify);
	open (my $DESCRIPTION, "<", $shift->mature_description->stringify);
	open (my $CONCATENTATED_ACCEPTED, ">>", $working_path."/all_accepted_${tag}_miRNAs.txt");
	open (my $CONCATENTATED_DISCARDED, ">>", $working_path."/all_discarded_${tag}_miRNAs.txt");
	open (my $CONCATENTATED_MATUREDES, ">>", $working_path."/all_mature_${tag}_miRNAs_description.txt");
	while (<$ACCEPTED>){
		chomp;
		print $CONCATENTATED_ACCEPTED "$_\n";
	}
	while (<$DISCARDED>){
		chomp;
		print $CONCATENTATED_DISCARDED "$_\n";
	}
	while (<$DESCRIPTION>){
		chomp;
		print $CONCATENTATED_MATUREDES "$_\n";
	}
	close $CONCATENTATED_ACCEPTED;
	close $CONCATENTATED_DISCARDED;
	close $CONCATENTATED_MATUREDES;
	close $ACCEPTED;
	close $DISCARDED;
	close $DESCRIPTION;
	return;
}

sub get_mature_map_location {
	my $coordFile = shift;
	# MI0005452 MIMAT0003844_2 MIMAT0004330 10 31 58 78 
	if (!-e $coordFile || -z $coordFile){
		print_error("Seems that $coordFile were not generated. Not possible to analyse RNA consensus structure");
	}
	open my $IN, "< $coordFile";
	my $total = 0; # Number of sequences
	my %ranges5p; #Count ranges 5p
	my %ranges3p; #Count ranges 5p
	while (<$IN>){
		chomp;
		$total++;
		my $fiveP = (split /\s+|\t/, $_)[3];
		my $threeP = (split /\s+|\t/, $_)[6];
		$ranges5p{$fiveP} += 1;
		$ranges3p{$threeP} += 1;	
	}
	# Get the most frequent
	my $start = get_most_frequent(\%ranges5p);
	my $end = get_most_frequent(\%ranges3p);
	my $final = "$start,$end";
	return ($final);
}

sub get_most_frequent {
	my $values = shift;
	my $min = 0;
	my $mostFreq;
	foreach my $key (sort keys %{$values}){
		my $number = $$values{$key};
		if ($number > $min){
			$min = $number;
			$mostFreq = $key;
		}
	}
	return $mostFreq;
}


no Moose;
__PACKAGE__->meta->make_immutable;
