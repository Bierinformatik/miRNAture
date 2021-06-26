package MiRNAnchor::Classify;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(process_all_candidates, concatenate_cms_paths);

use Moose;
use MooseX::Types::Path::Class;
use Data::Dumper;
use RNA;
use File::Copy;
use Bio::AlignIO;
use Bio::SimpleAlign;
use MiRNAnchor::Check;

use MiRNAture::ToolBox;

with 'MiRNAnchor::Tools';

has 'database_mirnas' => (
	is => 'rw',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'accepted_mirnas' => (
	is => 'rw',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'discarted_mirnas' => (
	is => 'rw',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
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

has 'precalculated_path' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1
);

has 'tag_spe_query' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'genome_specie' => (
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

has 'accepted_noStr_file' => (
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

has 'mature_description' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'gff_high_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'bed_high_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'gff_med_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'bed_med_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'gff_NO_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'bed_NO_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'gff_ACCEPTED_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'bed_ACCEPTED_file' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);
sub load_all_results {
	my $shift = shift;
	my $db_codes = load_database_miRNAs($shift->database_mirnas); #Complete miRNA table
	my $acceptedDB = load_accepted($shift->accepted_mirnas); #Partially accepted miRNAs
	my $discardedDB = load_discarted($shift->discarted_mirnas); # Discarted miRNAs
	return ($db_codes, $acceptedDB, $discardedDB);
}

sub process_all_candidates {
	my $shift = shift;
    my $variables = shift;
    my $name_cm_db = shift;
	my ($db_codes, $acceptedDB, $discardedDB) = load_all_results($shift); 
	my $finalValidationPath = $shift->output_folder;
    my $specie_tag = $shift->tag_spe_query;
    my $genome_path_complete = $variables->[3]->{Specie_data}{Old_Genome};
    my ($Zvalue, $genome_size_bp) = calculate_Z_value($genome_path_complete, "Genome");
    my $minBitscore = calculate_minimum_bitscore($genome_size_bp);
	# Output files
	open (my $ACCEPTED, ">", $shift->accepted_file) or die;
	open (my $ACCEPTED_NOALIGN, ">", $shift->accepted_noStr_file) or die; 
    open (my $DISCARTED, ">", $shift->discarded_file) or die; 
	##
	foreach my $code (sort keys %{$db_codes}){
		if (exists $$acceptedDB{$code}){ 
			#H1595270684	ENSLACT00000004505.1	+	RF00143	1	77	3309	3391	37.8	0.00066	no	mir-6	77	0	Infernal	Infernal	mir-6	miRNA
			my $family = (split /\s+|\t/, $$acceptedDB{$code})[3]; #RFAM family
			my $processGroup = (split /\s+|\t/, $$acceptedDB{$code})[-1]; #RFAM family process group: 0, 1, 2
			my $stockholmFile = "$finalValidationPath/$family/$code/Output/$family.out/$family.stk";	
			my $maturePosFile = "$finalValidationPath/$family/$code/Output/$family.out/$family-FinalCoor.txt";
			my $stockholmFileAlternative = "$finalValidationPath/$family/$code/Output/$family.out/${family}corrected.stk";
            my $fastaFinal = "$finalValidationPath/$family/$code/Output/$family.out/${family}-Final.fasta";
			# In cases where alterantive folding were created
			if (-e $stockholmFileAlternative && !-z $stockholmFileAlternative){
				$stockholmFile = $stockholmFileAlternative;
			}
			my ($stockholmFileOld, $stockholmFileOldAlternative);
			# In case exists a multifamily curated from RFAM
			if ($processGroup == 2){ #Usual family no groups
				$stockholmFileOld = $shift->precalculated_path."/".$family."/Output/".$family.".out/".$family.".stk";
				$stockholmFileOldAlternative = $shift->precalculated_path."/".$family."/Output/".$family.".out/".$family."corrected.stk";
			} elsif ($processGroup == 0){ #Selected group was 0 
				$stockholmFileOld = $shift->precalculated_path."/".$family."/Output_0/".$family.".out/".$family.".stk";
				$stockholmFileOldAlternative = $shift->precalculated_path."/".$family."/Output_0/".$family.".out/".$family."corrected.stk";
			} elsif ($processGroup == 1){ #Selected group was 1
				$stockholmFileOld = $shift->precalculated_path."/".$family."/Output_1/".$family.".out/".$family.".stk";
				$stockholmFileOldAlternative = $shift->precalculated_path."/".$family."/Output_1/".$family.".out/".$family."corrected.stk";
			} else {
				print_error("No valid process group, check Classify.pm:process_all_candidates.");
			}

			if (-e $stockholmFileOldAlternative && !-z $stockholmFileOldAlternative){
				$stockholmFileOld = $stockholmFileOldAlternative;
			}
            #
            ##### Extract final fasta and evaluate if truncated with global
            my $line = $$acceptedDB{$code};
            my $fasta_sequence = get_fasta_to_global($fastaFinal, $code, $finalValidationPath, $specie_tag);
            my $fam_sequence = (split /\s+|\t/, $line)[-2];
            my $truncated_scores_global = perform_global_evaluation_covariance($fasta_sequence, $name_cm_db, $finalValidationPath, $specie_tag, $code, $fam_sequence); 
            my ($truncated_value, $bitscore_global_value);
            #Truncated value
            if (exists $$truncated_scores_global{$code}{0}{Truncated}){
                $truncated_value = $$truncated_scores_global{$code}{0}{Truncated}; # Get the best result in position 0.
            } else {
                $truncated_value = "NotCalculated";
            }
            # Bitscore corrected
            if (exists $$truncated_scores_global{$code}{0}{Score}){
                $bitscore_global_value = $$truncated_scores_global{$code}{0}{Score}; # Get the best result in position 0.
            } else {
                $bitscore_global_value = -1;
            }
            ###
			# Check the conserved SS is valid
			my $conservedMatureLocations = get_mature_map_location($maturePosFile); #Get S5p,E5p,S3p,E3p from fam.
			my $validStructure = validate_secondary_structure_alignment($stockholmFile, $shift->current_dir, $conservedMatureLocations); #Evaluate generated alignment
			if ($validStructure =~ m/^Valid/){
				my $newStrLine = obtain_str_line($stockholmFile); #Get alignment structure line including new sequence
				my $oldStrLine = obtain_str_line($stockholmFileOld); #Structure line without new sequence
				#Obtain strings: Consensus, matches and structure
				my ($consensusOLD, $matchOLD, $structureOLD) = analysis_stockholm_file($stockholmFileOld); #Features OLD
				my ($consensusNEW, $matchNEW, $structureNEW) = analysis_stockholm_file($stockholmFile); #Features NEW
				#Normalize lengths conserved line
				$matchOLD = process_length($consensusOLD, $matchOLD); #Add characters if missing to match line
				$matchNEW = process_length($consensusNEW, $matchNEW);	
				my $strDistance = calculate_edit_distance_trees($newStrLine, $oldStrLine);
				my $new_data = $$acceptedDB{$code};
				my $old_data = $$db_codes{$code};
				# Compare distances from annotated and corrected
				my ($distanceStart, $distanceEnd, $distanceFinal) = compare_positions_miRNAs($new_data, $old_data);
                if ($truncated_value eq "no"){ # Global align evaluation, not truncated
                    if ($bitscore_global_value > $minBitscore){ #Look global score
                        if ($strDistance <= 7){
                            print $ACCEPTED "$$acceptedDB{$code}\t$validStructure\t$strDistance\t$distanceStart\t$distanceEnd\t$distanceFinal\t$truncated_value\t$bitscore_global_value\n";	
                        } else {
                            print $ACCEPTED_NOALIGN "$$acceptedDB{$code}\t$validStructure\t$strDistance\t$distanceStart\t$distanceEnd\t$distanceFinal\t$truncated_value\t$bitscore_global_value\n";	
                        }
                    } else { #Discard if global score is low
                        print $DISCARTED "$$db_codes{$code}\tNA\tDISCARDED\tNA\tNA\tNA\tNA\t$truncated_value\t$bitscore_global_value\n";	
                    }
                } else { # Truncated, in relation to the CM model
                    # Evaluate the bitscore, have to be > log2 2N
                    if ($bitscore_global_value > $minBitscore){ #Look global score
                        print $ACCEPTED_NOALIGN "$$acceptedDB{$code}\t$validStructure\t$strDistance\t$distanceStart\t$distanceEnd\t$distanceFinal\t$truncated_value\t$bitscore_global_value\n";	
                    } else {
                        print $DISCARTED "$$db_codes{$code}\tNA\tDISCARDED\tNA\tNA\tNA\tNA\t$truncated_value\t$bitscore_global_value\n";	
                    }
                }
			} elsif ($validStructure =~ m/^No_valid|^NA|^No_valid_No_match/){
                if ($truncated_value eq "no"){ # Global align evaluation, not truncated
                    if ($bitscore_global_value > $minBitscore){ #Look global score
                        print $ACCEPTED_NOALIGN "$$acceptedDB{$code}\t$validStructure\tNA\tNA\tNA\tNA\t$truncated_value\t$bitscore_global_value\n"; 
                    } else {
                        print $DISCARTED "$$db_codes{$code}\tNA\tDISCARDED\tNA\tNA\tNA\tNA\t$truncated_value\t$bitscore_global_value\n";	
                    }
                } else {
                    if ($bitscore_global_value > $minBitscore){ #Discuss this one. No valid str + truncated, but valid score
                        print $ACCEPTED_NOALIGN "$$acceptedDB{$code}\t$validStructure\tNA\tNA\tNA\tNA\t$truncated_value\t$bitscore_global_value\n"; 
                    } else {
                        print $DISCARTED "$$db_codes{$code}\tNA\tDISCARDED\tNA\tNA\tNA\tNA\t$truncated_value\t$bitscore_global_value\n";	
                    }
                }
			} else {
				print_error("There is an error in the syntax code $validStructure and $code");
			}
		} elsif (exists $$discardedDB{$code}){ #This miRNA was classified as discarted
			print $DISCARTED "$$db_codes{$code}\tNA\tDISCARDED\tNA\tNA\tNA\tNA\tNA\tNA\n";	
		} else {
			print $DISCARTED "$$db_codes{$code}\tNA\tCM_MODEL_DISCARDED\tNA\tNA\tNA\tNA\tNA\tNA\n";	
		}
	}
	close $ACCEPTED;
	close $ACCEPTED_NOALIGN;
	close $DISCARTED;
	return;
}

sub perform_global_evaluation_covariance {
    my ($file, $pathsCM, $outFolder, $specieT, $code, $family) = @_;
    my $outfile;
	if (!-e $file || -z $file){
		print_error("Sequence file to global evaluation not created\n");
	} 
    if (exists $$pathsCM{$family}){
        my $cm = $$pathsCM{$family};
        $outfile = cmsearch_global($family, $cm, $file, $specieT, $outFolder, "cmsearch", $code); #, $zscore); 
    } else {
        print_error("CM is $family not available for global validation\n");
    }
    # Here concatenate and classify the results if they are complete or not
    #my $outfile = concatenate_tab_files($outFolder, $specieT);
    # Delete individual files
    #delete_individual_files($outFile, "\\_global\.tab");
    #delete_individual_files($outFile, "\\_global\.out");
    # Build hash with resulting data
    my $truncated_global = get_result_truncated_global($outfile);
    return $truncated_global;
}

sub concatenate_tab_files {
    my ($outFolder, $specie) = @_;
    my $concatenatedFile = "$outFolder/${specie}_concatenated_evaluation.tab";
    open my $OUT, "> $concatenatedFile" or die;
    my @filesToConcatenate = check_folder_files($outFolder, "\\_global\.tab");
    foreach my $out (@filesToConcatenate){
        open my $IN, "< $out" or die;
        while (<$IN>){
            chomp;
            print $OUT "$_\n";
        }
    }
    close $OUT;
    return $concatenatedFile; 
}

sub cmsearch_global {
	#In: out, CM model, genome
	my ($nameCMFinal, $cm, $sequence, $speTag, $outFolder, $cmsearch_path, $code) = @_; #modes: 1 X.cm 0 other name
	existenceProgram($cmsearch_path);
	$sequence =~ s/"//g;
	if (-e "${cm}" && !-z "${cm}"){
        # Here detect truncated sequences
        # Z-score?
        my $param = "-g --cpu 5 --toponly --nohmmonly --tblout $outFolder/${nameCMFinal}_${speTag}_${code}_global.tab -o $outFolder/${nameCMFinal}_${speTag}_${code}_global.out $cm $sequence";
		system "$cmsearch_path $param 1> /dev/null";
	} else {
		print_error("The CM to global evaluation is not available\n");
	}
    my $outfile = "$outFolder/${nameCMFinal}_${speTag}_${code}_global.tab";
	return $outfile;
}

sub generate_temporal_fasta_family {
    my ($dir, $sequences, $family, $specie, $code) = @_;
    my $out_name = "$dir/${family}_${specie}_${code}_temporal_evaluate.fa";
    open my $OUT, "> " or die;
    foreach my $head (sort keys %{$sequences}){
        print $OUT "$head\n";
        print $OUT "$$sequences{$head}\n";
    }
    close $OUT;
    return $out_name;
}

sub get_result_truncated_global {
    my ($file) = @_;
    open my $IN, "< $file";
    my %data;
    my $count = 0; #Track results
    while (<$IN>){
        next if ($_ =~ /^#/);
        my @complete = split /\s+|\t/, $_;
        $data{$complete[0]}{$count}{Truncated} = $complete[10]; #Truncated result
        my $scoreB = $complete[14] - $complete[13]; # Bitscore - bias
        $data{$complete[0]}{$count}{Score} = $scoreB; # Bitscore
        $count++;
    }
    return \%data;
}

sub concatenate_cms_paths {
    my ($default, $user, $other) = @_;
    my @paths;
    my %final;
    if (-d $default) {
        push @paths, $default;
    }
    if (-d $other) {
        push @paths, $other;
    }
    if (-d $user) {
        push @paths, $user;
    }
    return \@paths;
}

sub concatenate_cms_paths_old {
    my ($default, $user, $other, $fam) = @_;
    my @paths;
    my %final;
    push @paths, $default;
    push @paths, $other;
    push @paths, $user;
    foreach my $in (@paths){
        my $cm = "$in/$fam.cm";
        if (-e $cm && !-z $cm){
            $final{$fam} = $cm; #Only take one model. If repeated, take last one.
        }
    }
    return \%final;
}

sub get_fasta_to_global {
    my ($fasta, $codeStart, $finalValidationPath, $specieTag) = @_;
    open my $IN, "< $fasta" or die;
    my $head;
    my $code;
    my %sequences;
    while (<$IN>){
        chomp;
        if ($_ =~ /^>/){
            $head = $_;
            my @headm = split /\s+|\t/, $head;
            $code = (split /\s+|\t/,$head)[1];
            $head = "\>$headm[1] $headm[2] $headm[3] $headm[4] $headm[5]";
        } elsif ($_ !~ /^>/ ) {
            if ($code eq $codeStart){
                $sequences{$head} .= $_;
            } else {
                next;
            }
        }
    }
    close $IN;
    my $outfile = "$finalValidationPath/${codeStart}_${specieTag}_temporal_file.fa";
    open my $TEMP, "> $outfile" or die;
    foreach my $head (sort keys %sequences){
        print $TEMP "$head\n$sequences{$head}\n";
    }
    close $TEMP;
    return $outfile;
}

sub edit_pairings {
	my ($group_elements, $delete, $sense, $symbol) = @_;
	if ($sense eq "start"){ #Remove symbol from start of the block
		my $x = 1;
		for (my $p = 0; $p <= $#$group_elements; $p++){
			#foreach my $ln (@$group_elements){
			if ($x <= $delete - 1){ 
				if ($$group_elements[0] eq $symbol){
					shift @$group_elements; #Remove the first element
					$x++;
				} else {
					shift @$group_elements; #Remove the first element
				}
			} else {
				last;
			}
		}
	} elsif ($sense eq "end"){ #Remove from the end
		my $x = 1;
		my $reversed = reverse @$group_elements;
		my @reversedAll =  split //, $reversed;
		my $numberElements = scalar @reversedAll;
		for (my $p = 0; $p <= $#reversedAll; $p++){
			if ($x <= $delete -1){
				#while ($numberElements && $x <= $delete){
				if ($reversedAll[0] eq $symbol){
					shift @reversedAll; #Remove the first element
					$x++;
				} else {
					shift @reversedAll; #Remove the first element
				}
			} else {
				last;
			}
		}
		$group_elements = reverse @reversedAll;
		return $group_elements;
	} else {
		die "Not right direction!\n";
	}
	my $cleaned_block = join "", @$group_elements;
	return $cleaned_block;
}

sub count_elements {
	my ($group_elements, $sign) = @_;
	my $count = 0;
	for ($a = 0; $a <= scalar @$group_elements - 1; $a++){
		if ($$group_elements[$a] eq "$sign"){
			$count++;
		} else {
			next;
		}
	}
	return $count;
}

sub clean_unpaired_ends {
	my ($startBlock, $endBlock) = @_;
	my @array1 = split //, $startBlock;
	my @array2 = split //, $endBlock;
	my $numA = count_elements(\@array1, "\(");
	my $numB = count_elements(\@array2, "\)");
	if ($numA == $numB){
		return ($startBlock, $endBlock);
	} elsif ($numA < $numB){
		my $diff = abs($numA - $numB) + 1;
		$endBlock = edit_pairings(\@array2, $diff, "end", "\)");
		return ($startBlock, $endBlock);
	} elsif ($numA > $numB){
		my $diff = abs($numA - $numB) + 1;
		$startBlock = edit_pairings(\@array1, $diff, "start", "\(");
		return ($startBlock, $endBlock);
	}
}

sub generate_conserved_blocks {
	my ($consensus, $match, $structure, $symbolStart, $symbolEnd, $code) = @_; 
	#InferBlocks: Slice numbers
	my ($start_block1, $end_block1, $start_block2, $end_block2) = infer_blocks($consensus, $match, $structure, "\(", "\)",$code);
	if ($start_block1 eq "NA"){
		return "NA";
	}
	#IsolateBlocks: Cut each string with slice numbers
	my $block1consensus = isolate_block($consensus, $start_block1, $end_block1);
	my $block2consensus = isolate_block($consensus, $start_block2, $end_block2);
	my $block1match = isolate_block($match, $start_block1, $end_block1);
	my $block2match = isolate_block($match, $start_block2, $end_block2);
	my $block1structure = isolate_block($structure, $start_block1, $end_block1);
	my $block2structure = isolate_block($structure, $start_block2, $end_block2);	
	#PolishBlocks: Refine restincting sliced strings to conserved ends
	my ($consensusS, $matchS, $structureS) = polish_block($block1consensus, $block1match, $block1structure, "\*");
	my ($consensusE, $matchE, $structureE) = polish_block($block2consensus, $block2match, $block2structure, "\*");
	return ($consensusS, $consensusE, $matchS, $matchE, $structureS, $structureE);
}

sub process_length {
	my ($consensus, $match) = @_;
	$match =~ s/\s/N/g;
	if (length $consensus > length $match){
		my $len1 = length $consensus;
		my $len2 = length $match;
		my $diff = $len1 - $len2;
		$match = add_characters($match, $diff);
	} elsif (length $consensus < length $match){
		print_error("Match is larger than Consensus, this not possible: $consensus and $match");
	}
	return $match;
}

sub polish_block {
	my ($consensus, $match, $structure, $symbol) = @_;
	my @mat = split //, $match; #Conserved columns
	my $len4 = scalar @mat;
	my $switchA = 0; #Turns 1 if found symbolStart. Only to save the first one.
	my $start_block1 = 0;
	my $end_block1 = 0;
	for ($a = 0; $a <= $len4 - 1; $a++){
		if ($mat[$a] eq $symbol && $switchA == 0){      
			$start_block1 = $a;    
			$switchA = 1;             
		}
		if ($mat[$a] eq $symbol){
			$end_block1 = $a;
		} else {
			next;
		}    
	}
	my $distance = abs($start_block1 - $end_block1) + 1;
	if ($distance == 1){
		$consensus = "";
		$match = "";
		$structure = "";
	} else {
		$consensus = substr ($consensus, $start_block1, $distance);
		$match = substr ($match, $start_block1, $distance);
		$structure = substr ($structure, $start_block1, $distance);
	}
	return ($consensus, $match, $structure);
}   

sub isolate_block {
	my ($string, $start, $end) = @_;
	my $distance = abs($end - $start) + 1;
	my $block = substr($string, $start, $distance);
	return $block;
}

sub add_characters {
	my ($string, $number) = @_;
	my $complement = "N" x $number;
	$string = $string.$complement;
	return $string;
}

sub infer_blocks {
	my ($consensus, $match, $structure, $symbolStart, $symbolEnd, $code) = @_;
	my @con = split //, $consensus; #Nucleotides
	my @mat = split //, $match; #Conserved columns
	my @struc = split //, $structure; #SS structure
	my $len1 = scalar @con;
	my $len2 = scalar @mat;
	my $len3 = scalar @struc;
	my $oper = ($len1 - $len2) - $len3; #a - a - a = -a
	my ($start_block1, $start_block2);
	my $end_block1 = 0; my $end_block2 = 0;
	if ($oper == (-1 * $len1)){ #Test if all variables have the same length: Have to!
		my $switchA = 0; #Turns 1 if found symbolStart. Only to save the first one.
		my $switchB = 0; #Turns 1 if found a symbolEnd. Only to save the first one.
		for ($a = 0; $a <= $len3 - 1; $a++){
			if ($struc[$a] eq $symbolStart && $switchA == 0){
				$start_block1 = $a;
				$switchA = 1;
			}
			if ($struc[$a] eq $symbolEnd && $switchB == 0){
				$start_block2 = $a;
				$switchB = 1;
			}
			if ($struc[$a] eq $symbolStart){
				$end_block1 = $a;
			} elsif ($struc[$a] eq $symbolEnd){
				$end_block2 = $a;
			} else {
				next;
			}
		}
	} else {
		print_error("Lenghts of strings are not the same in $code");
	}
	return ($start_block1, $end_block1, $start_block2, $end_block2);
}

sub analysis_stockholm_file {
	my $stoFile = shift;
	my $structure;
	my $str = Bio::AlignIO->new(-file => $stoFile, -format => 'stockholm');
	my $aln = $str->next_aln();
	my $consensusNT = $aln->consensus_string(1); # Get the consensus along the alignment, 1 means 1% of Identity, to take all nts 
	my $matchComplete = $aln->match_line(); # Symbols indicating strong, medium or low conservation
	open my $IN, "< $stoFile" or die;
	while (<$IN>){
		chomp;
		if($_ =~ /^\#\=GC SS_cons/){
			$_ = (split /\s+|\t/, $_)[2];
			$structure = $_; #Structure line
			last;
		}
	}											
	close $IN;
	return ($consensusNT, $matchComplete, $structure);
}

sub calculate_edit_distance_trees {
	my ($new, $old) = @_;
	my $tree1 = structure_parsing($old);
	my $tree2 = structure_parsing($new);
	my $editDistance = RNA::tree_edit_distance($tree1, $tree2);
	return $editDistance;
}

sub structure_parsing {
	my $line = shift;
	my $conversion = RNA::b2C($line); #Convert to Coarse-grained nomenclature
	my $tree = RNA::make_tree($conversion); #Convert to tree
	return $tree;
}

sub obtain_str_line {
	my $file =  shift;
	open my $IN, "< $file" or die "The Stockholm file $file is empty\n";
	my $line;
	while (<$IN>){
		chomp;
		if ($_ =~ /^\#\=GC SS_cons/){
			$line = (split /\s+|\t/, $_)[2]; #Get Structure line
			$line =~ s/\#//g;
			last;
		} else {
			next;
		}
	}
	return $line;
}

sub check_structure {
	my $file =  shift;
	open my $IN, "< $file" or print "The Stockholm file is empty\n";
	if (-z $file || !-s $file){
		return 0;
	}
	while (<$IN>){
		chomp;
		if ($_ =~ /^\#\=GC SS_cons/){
			if ($_ =~ m/\(|\)/){
				return 1;	
			} else {
				return 0;
			}	
		} else {
			next;
		}
	}
}

sub load_accepted {
	my $table = shift;
	my $db = process_table($table); 
	return $db;
}

sub load_discarted {
	my $table = shift;
	my $db = process_table($table);
	return $db;
}

sub process_table {
	my $table = shift;
	my %temp;
	open my $FILE, "< $table" or die;
	while (<$FILE>){
		chomp;
		my @table = split /\s+|\t/, $_;
		$temp{$table[0]} = $_;		
	}
	return \%temp;
}


sub compare_positions_miRNAs {
	my ($newLine, $oldLine) = @_;	
	my ($SOld, $SNew, $EOld, $ENew, $Sdiff, $Ediff);
	$SOld = (split /\s+|\t/, $oldLine)[6];
	$EOld = (split /\s+|\t/, $oldLine)[7];
	$SNew = (split /\s+|\t/, $newLine)[6];
	$ENew = (split /\s+|\t/, $newLine)[7];
	$Sdiff = ($SOld - $SNew) + 1;
	$Ediff = ($EOld - $ENew) + 1;
	my $distanceFinal = abs($ENew - $SNew) + 1;
	return ($Sdiff, $Ediff, $distanceFinal);
}

sub generate_output_files {
	my $shift = shift;
	open (my $POSITIVE, "<", $shift->accepted_file->stringify);
	open (my $POSITIVE_NO_STR, "<", $shift->accepted_noStr_file->stringify);
	open (my $DISCARDED, "<", $shift->discarded_file->stringify);
	open (my $INFO, "<", $shift->mature_description->stringify);
	my (%info, %positiveHigh, %positiveMed, %positiveAccepted, %negativeLow);
	#H1601739630	X	-	RF00711	1	90	45932312	45932381	17.0	0.0034	no	mir-449	90	10	Blast	ANCA,DARE,XELA,XETR	anca148299,anca148674,anca149012,dare322,dare719,xela37007,xetr2120,xetr2352	miRNA	1
	#H1601747417	>patr-RF00711-42 H1601747417 Pan troglogytes mir-449 stem-loop H1601747417 MIMAT0041140/B1600189418/H1601747417 MIMAT0041140/B1600189418/H1601747417-star 10 31 45 67 RF00711 RF00711_1	Accepted	RF00711_1
	#H1601739630	>patr-RF00711-35 H1601739630 Pan troglogytes mir-449 stem-loop NA RF00711 RF00711_0	Discarded	RF00711_0
	while (<$INFO>){
		chomp;
		my @all = split /\s+|\t/, $_;
		my $id = $all[0];
		my $fam = $all[5];
		my $rfamily = (split /\-/, $all[1])[1];
		my $start5p = $all[10];
		my $end5p = $all[11];
		my $start3p = $all[12];
		my $end3p = $all[13];
		my $classification = $all[-2];
		$info{$id}{$classification} = [ $fam, $rfamily, $start5p, $end5p, $start3p, $end3p ];
	}
	while (<$POSITIVE>){
		chomp;
		my @all = split /\s+|\t/, $_;
		my $id = $all[0];
		my $chr = $all[1];
		my $strand = $all[2];
		my $start = $all[6];
		my $end = $all[7];
		my $evalue = $all[9];
		my $idfam = $all[11];
		my $idrfam = $all[3];
		my $support = "High_confidence";
		$positiveHigh{$id} = [$chr, $strand, $start, $end, $idfam, $idrfam, $evalue, $support];
		$positiveAccepted{$id} = [$chr, $strand, $start, $end, $idfam, $idrfam, $evalue, $support];
	}
	while (<$POSITIVE_NO_STR>){
		chomp;
		my @all = split /\s+|\t/, $_;
		my $id = $all[0];
		my $chr = $all[1];
		my $strand = $all[2];
		my $start = $all[6];
		my $end = $all[7];
		my $evalue = $all[9];
		my $idfam = $all[11];
		my $idrfam = $all[3];
		my $support = "Medium_confidence";
		$positiveMed{$id} = [$chr, $strand, $start, $end, $idfam, $idrfam, $evalue, $support];
		$positiveAccepted{$id} = [$chr, $strand, $start, $end, $idfam, $idrfam, $evalue, $support];
	}
	while (<$DISCARDED>){
		chomp;
		my @all = split /\s+|\t/, $_;
		my $id = $all[0];
		my $chr = $all[1];
		my $strand = $all[2];
		my $start = $all[6];
		my $end = $all[7];
		my $evalue = $all[9];
		my $idfam = $all[11];
		my $idrfam = $all[3];
		my $support = "NO_confidence";
		$negativeLow{$id} = [$chr, $strand, $start, $end, $idfam, $idrfam, $evalue, $support];
	}
	#Generate Accepted (High + Medium) output
	generate_gff($shift->gff_ACCEPTED_file, \%positiveAccepted, \%info);
	generate_bed($shift->bed_ACCEPTED_file, \%positiveAccepted, \%info);
	#Genererate High Output
	#generate_gff($shift->gff_high_file, \%positiveHigh, \%info);
	#generate_bed($shift->bed_high_file, \%positiveHigh, \%info);
	#Genererate Medium Output
	#generate_gff($shift->gff_med_file, \%positiveMed, \%info);
	#generate_bed($shift->bed_med_file, \%positiveMed, \%info);
	#Genererate Low Output
	generate_gff($shift->gff_NO_file, \%negativeLow, \%info);
	generate_bed($shift->bed_NO_file, \%negativeLow, \%info);
	return;
}

sub generate_gff {
	my ($out, $positive, $info) = @_;
	open (my $GFF, ">", $out) or die "GFF file fails to be created\n";
	if (-e $out && -z $out){ #Print header
		print $GFF "\#\#gff-version 3\n";
	}
	foreach my $id (sort keys %{$positive}){
		# Check if exist in positive
		print $GFF "$$positive{$id}[0]\tmiRNAture\tpre_miRNA\t$$positive{$id}[2]\t$$positive{$id}[3]\t$$positive{$id}[-2]\t$$positive{$id}[1]\t\.\tID=$id;NAME=$$positive{$id}[4];Alias=$$positive{$id}[5];Support=$$positive{$id}[-1]\n"; 
		if (exists $$info{$id} && exists $$info{$id}{"Accepted"}){
			my $start5p = $$positive{$id}[2] + $$info{$id}{"Accepted"}[2]; 
			my $end5p = $$positive{$id}[2] + $$info{$id}{"Accepted"}[3];
			my $start3p = $$positive{$id}[2] + $$info{$id}{"Accepted"}[4];
			my $end3p = $$positive{$id}[2] + $$info{$id}{"Accepted"}[5];
			print $GFF "$$positive{$id}[0]\tmiRNAture\tmiRNA\t$start5p\t$end5p\t\.\t$$positive{$id}[1]\t\.\tID=${id}_5p;NAME=$$positive{$id}[4];Alias=$$positive{$id}[5];Parent=$id;Comments=$$positive{$id}[-1]\n"; 
			print $GFF "$$positive{$id}[0]\tmiRNAture\tmiRNA\t$start3p\t$end3p\t\.\t$$positive{$id}[1]\t\.\tID=${id}_3p;NAME=$$positive{$id}[4];Alias=$$positive{$id}[5];Parent=$id;Comments=$$positive{$id}[-1]\n";
		} else {
			; #Now the negative is considered, then no matures could be detected, only print precursors 
		}
	}
	close $GFF;
	return;
}

sub generate_bed {
	my ($out, $positive, $info) = @_;
	open (my $BED, ">", $out) or die "GFF file fails to be created\n";
	foreach my $id (sort keys %{$positive}){
		# Check if exist in positive
		# Change from 1-based to 0-based
		$$positive{$id}[2] = $$positive{$id}[2] - 1;
		$$positive{$id}[3] = $$positive{$id}[3] - 1;
		print $BED "$$positive{$id}[0]\t$$positive{$id}[2]\t$$positive{$id}[3]\t$$positive{$id}[4]\t0\t$$positive{$id}[1]\n"; 
		if (exists $$info{$id} && exists $$info{$id}{"Accepted"}){
			my $start5p = $$positive{$id}[2] + $$info{$id}{"Accepted"}[2]; 
			my $end5p = $$positive{$id}[2] + $$info{$id}{"Accepted"}[3];
			my $start3p = $$positive{$id}[2] + $$info{$id}{"Accepted"}[4];
			my $end3p = $$positive{$id}[2] + $$info{$id}{"Accepted"}[5];
			$start5p = $start5p - 1;
			$end5p = $end5p - 1;
			$start3p = $start3p - 1;
			$end3p = $end3p - 1;
			print $BED "$$positive{$id}[0]\t$start5p\t$end5p\t$$positive{$id}[4]_5p\t0\t$$positive{$id}[1]\n"; 
			print $BED "$$positive{$id}[0]\t$start3p\t$end3p\t$$positive{$id}[4]_3p\t0\t$$positive{$id}[1]\n";
		} else {
			;
		}
	}
	close $BED;
	return;
}

sub organise_mess {
	my $shift = shift;
	my $working_path = shift;
	create_folders($working_path, "GFF3");
	create_folders($working_path, "BED");
	create_folders($working_path, "Additional_Support");
    my $short_specie = $shift->tag_spe_query;
	opendir DH, $working_path;
	while (my $file = readdir DH) {
		if ($file =~ /\.txt$/ && $file =~ /\_$short_specie/){
			move("$working_path/$file", "$working_path/Additional_Support");
		} elsif ($file =~ /\.gff3/ && $file =~ /\_$short_specie/){
			move("$working_path/$file", "$working_path/GFF3");
		} elsif ($file =~ /\.bed/ && $file =~ /\_$short_specie/){
			move("$working_path/$file", "$working_path/BED");
		} else {
			;
		}
	}
	return;
}

no Moose;
__PACKAGE__->meta->make_immutable;
