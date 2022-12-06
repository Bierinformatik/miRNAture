package MiRNAture::_BlastSearch_ResolveBlastMergings;

use Moose::Role;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);

my (@new, @sorted, @temporal, @all, @new_coord);
my $percentage_limit_defined = 40; #Percentage of minimum coverage from subject to query.
my $F_OUT;
my $F_OUT2;

# Subs
sub resolveBlastMergings {
	#hmiRNA127	scaffold1431-size16679	6215	6166	-1	50	76.00	1	49	60
	my $shift = shift;
	my $input_path_db_file = shift;
	my $output_merged = shift;
	my $output_special = shift;
	my %data;
	open (my $FILE, "<", $input_path_db_file) or die;
	open ($F_OUT, ">", $output_merged) or die;
	open ($F_OUT2, ">", $output_special) or die;
	#Check if empty
	my $input = $input_path_db_file;
	if (-z $input || !-e $input){
		return;
	}
	#Read File line by line
	while (<$FILE>){
		chomp;
		my @fields = split /\s+|\t/, $_;
		if ($fields[4] == -1){
			$fields[4] = "-";
		} elsif ($fields[4] == 1) {
			$fields[4] = "+";
		} else {
			die "Seems that your columns are misplaced!\n";
		}
		push @{ $data{"$fields[0]\t$fields[1]\t$fields[4]"} }, [ @fields[2,3,7,8,9,0,1,4] ];
	}
	foreach my $key ( sort keys %data){
		my $val = $data{$key};
		my $num = scalar @$val;
		direct_flow($num, $val, $F_OUT);
		CHANGE: 
		next;
	}
	close $F_OUT; close $F_OUT2; close $FILE;  
	return;
}

sub direct_flow {
	my ($number, $values, $F_OUT) = @_;
	if ($number == 0){
		die "Exists a fatal error related with number of values and keys from hash of candidates!\n";
	} elsif ($number == 1 ){
		#Check if have the coverage!
		foreach my $i (${$values}[0]){
			my $test_cov = test_coverage($$i[0], $$i[1], $$i[4]);
			if ($test_cov == 1){ #This register reported a 40 \% or more of coverage respect to query
				#sce1 dvexscf18546 +	 3359	 3424	 8	 73	 82
				print $F_OUT "$$i[5]\t$$i[6]\t$$i[7]\t$$i[0]\t$$i[1]\t$$i[2]\t$$i[3]\t$$i[4]\n";
			}
		}
	} else { #More than 1
		@sorted = sort { $a->[0] <=> $b->[0] } @$values; #Organize by the first one coordinate
		@all = @sorted;
		analyse_sorted(\@sorted);
	}
}

sub analyse_sorted {
	my $all_sorted = shift;
	@temporal = ();
	my ($S1, $E1, $S2, $E2);
	FLAG:
	push @temporal, $$all_sorted[0]; #Pair comparisons.
	push @temporal, $$all_sorted[1];
	FLAG1:
	my $sense0 = $$all_sorted[0][7];
	my $sense1 = $$all_sorted[1][7];
	$S1 = $temporal[0][0];
	$E1 = $temporal[0][1];
	$S2 = $temporal[1][0];
	$E2 = $temporal[1][1];
	my $overlap = compare_coord($S1,$E1,$S2,$E2);
	if ($overlap == 1){ #Coordinates from candidates are overlapping
		@new_coord = fusion_all($temporal[0], $temporal[1], 1);
		@new = ();
		filter_candidate(\@new_coord, $temporal[0], $temporal[1]); #Fusion and individual queries	
	} else { #There is no-overlap, so print the first if exists x > 2 elements
		to_print($temporal[0], 1);
		@temporal = ();
		shift @all; #Get first one from all list
		shift @sorted;	#Get first one from sorted list	
		test_length(\@all, 1); #number of elements on temporal
	}
	return;
}

sub compare_coord {
	my ($s1,$e1,$s2,$e2) = @_;
	my $response;
	if ( ($s1 <= $s2 && $e1 >= $e2) || ($s1 <= $s2 && $s2 <= $e1 && $e1 <= $e2) || ($s1 >= $s2 && $e2 >= $s1 && $e2 <= $e1)){
		$response = 1;
	} else {
		$response = 0;
	}
	return $response;
}

sub test_coverage {
	#$$i[0], $$i[1] $$i[2], $$i[3], $$i[4]
	my ($startS, $endS, $lenQ) = @_;
	my $limit = sprintf("%.0f",(($percentage_limit_defined * $lenQ)/100));
	my $upper_limit = sprintf("%.0f",((125 * $lenQ)/100)); #defined upper limit, maximum length
	my $dist_q = abs($endS - $startS) + 1;
	if (($dist_q >= $limit) && ($dist_q <= $upper_limit)){
		return 1;
	} else {
		return 0;
	}
}

sub fusion_all {
	my @cand1 = $_[0];
	my @cand2 = $_[1];
	my $mode = $_[2];
	my $new_line;
	if ($mode == 1){ #Merge with greater values and names with comma
		for (my $a=0; $a <= 7; $a++){ # Walking througt the array
			if ($a == 0 || $a == 2){ #Select the lowest, is the best.
				push_best($cand1[0][$a], $cand2[0][$a], -1);
			} elsif ($a == 1 || $a == 3 || $a == 4){ #Greater value is the best
				push_best($cand1[0][$a], $cand2[0][$a], 1);
			} elsif ($a == 5 || $a == 6 || $a == 7){ #Text comparison, merge by comma if different.
				push_best($cand1[0][$a], $cand2[0][$a], 0);
			} else { 
				die "It contains more fields to merge, check the input format\n";
			}
		}	
	}
	return @new;
}

sub test_equal {
	my ($cand1, $cand2) = @_;
	for (my $i=0; $i<=7;$i++){
		if (looks_like_number($$cand1[$i])){
			unless ($$cand1[$i] == $$cand2[$i]){
				return 0;
			}
		} else {# String comparison 
			unless ($$cand1[$i] eq $$cand2[$i]){
				return 0;
			}
		}
	}
	return 1;
}

sub filter_candidate {
	my @data = $_[0];
	my $cand1 = $_[1];
	my $cand2 = $_[2];
	my $eqTest = test_equal($cand1, $cand2);
	if ($eqTest == 1){ # The candidates are equal, remove first one
		@temporal = ();
		shift @sorted; #Remove [0]
		shift @all; #Remove [0]
		@new_coord = ();
		test_length(\@all, 1);
	} else {
		# Candidates differ
		my $evaluation_candidate = test_coverage($data[0][0], $data[0][1], $data[0][4]);
		#Evaluate coverage 40% <= x <= 120% 
		if($evaluation_candidate == 1){ #Fusion is accomplished
			@temporal = ();
			shift @sorted; #Remove [0]
			shift @sorted; #Remove [1]
			shift @all; #Remove [0]
			shift @all; #Remove [1]
			# Insert to the array the new mergen candidate
			unshift @sorted, [ @new_coord ]; 
			unshift @all, [ @new_coord ];
			@new_coord = ();
			test_length(\@all, 1);
		} elsif ($evaluation_candidate == 0) {
			#Here the evaluation of the merged candidates is not the best,
			#In this case, evaluate independently the candidates
			# If useful, print the candidates to special evaluation
			my $best = compare_overlapping_excess($cand1, $cand2);
			# Select the best candidate and discard the other to special.
			if ($best == 0){ # Print to special the oposite number: if the best
				# is 0, print as special 1.
				to_print($cand2, 0); # Print to SPECIAL file
				@temporal = ();
				shift @sorted; #Remove [0]
				shift @sorted; #Remove [1]
				shift @all; #Remove [0]
				shift @all; #Remove [1]
				# Mantain cand1 on the evaluation
				unshift @all, [ @$cand1 ];
				unshift @sorted, [ @$cand1 ]; 
				@new_coord = ();
				test_length(\@all, 1);
			} elsif ($best == 1){
				to_print($cand1, 0); # Print to SPECIAL file
				@temporal = ();
				shift @sorted; #Remove [0]
				shift @sorted; #Remove [1]
				shift @all; #Remove [0]
				shift @all; #Remove [1]
				unshift @all, [ @$cand2 ];
				unshift @sorted, [ @$cand2 ]; 
				@new_coord = ();
				test_length(\@all, 1);
			} else { #The lengths are equal print both as special, should be a long region
				@temporal = ();
				@new_coord = (); ##
				shift @sorted; #Remove [0]
				shift @sorted; #Remove [1]
				shift @all; #Remove [0]
				shift @all; #Remove [1]
				test_length(\@all, 1);
			}
		} else {
			die "Something went wrong on evaluation of new candidate!\n";
		}
	}
}	

sub compare_overlapping_excess {
	# Select the best between both candidates, the coverage of the merge is not 
	# enough or excess the distance.
	my ($first, $second) = @_;
	my $best;
	my $leng1 = (abs($$first[2] - $$first[3]))/$$first[4];
	my $leng2 = (abs($$second[2] - $$second[3]))/$$second[4];
	if ($leng1 >= $leng2){ #Compare relative coverage, respect to query
		$best = 0;
	} elsif ($leng1 == $leng2){
		$best = 2;
	} else {
		$best = 1;
	}
	return $best;
}

sub to_print {
	my $complete = $_[0];
	my $out = $_[1];
	my $test_cov = test_coverage($$complete[0], $$complete[1], $$complete[4]);		
	if ($out == 1){ #Print as candidate
		if ($test_cov == 1){ #Pass coverage test, print.
			print $F_OUT "$$complete[5]\t$$complete[6]\t$$complete[7]\t$$complete[0]\t$$complete[1]\t$$complete[2]\t$$complete[3]\t$$complete[4]\n";
		} else { #Do not print!, do not have the coverage!
			;
		}
	} else {
		print $F_OUT2 "$$complete[5]\t$$complete[6]\t$$complete[7]\t$$complete[0]\t$$complete[1]\t$$complete[2]\t$$complete[3]\t$$complete[4]\n";
	}
}

sub test_length {
	my $array = $_[0]; # Test the length of the big array
	my $numb = $_[1]; #Mode
	my $nelements = scalar @$array;
	if ($nelements == 0){
		goto CHANGE;
	} elsif ($nelements == 1){
		foreach my $i (${$array}[0]){
			my $test_cov = test_coverage($$i[0], $$i[1], $$i[4]);
			if ($test_cov == 1){ #This register reported a 40 \% or more of coverage respect to query
				#sce1 dvexscf18546 +	 3359	 3424	 8	 73	 82
				print $F_OUT "$$i[5]\t$$i[6]\t$$i[7]\t$$i[0]\t$$i[1]\t$$i[2]\t$$i[3]\t$$i[4]\n";
			}
			@temporal = ();
		}
	} else {
		if ($numb == 1){ # Mode
			goto FLAG;
		} else {
			goto FLAG1;
		}
	}
}

sub push_best {
	my ($test1, $test2, $ind) = @_;
	my $best;
	if ($ind == -1){ #menor
		if ($test1 <= $test2){
			$best = 0;
		} else {
			$best = 1;
		}
	} elsif ($ind == 1){ #mayor
		if ($test1 <= $test2){
			$best = 1;
		} else {
			$best = 0;
		}
	} elsif ($ind == 0){ #texto
		if ($test1 eq $test2){
			$best = 2; #$test1;	
		} else {
			$best = 3; #"($test1,$test2)";
		}
	} else {
		die "Wrong definitions on merging rules!\n";
	}
	push_to_new($best, $test1, $test2);
}

sub push_to_new {
	my ($cand, $text1, $text2) = @_;
	if ($cand == 0){
		push @new, $text1;
	} elsif ($cand == 1){
		push @new, $text2;
	} elsif ($cand == 2){ #Text are the same
		push @new, $text1;
	} elsif ($cand == 3){ #Join text, those are not the same
		my $all = "$text1,$text2";
		push @new, $all;
	} else {
		die "Error creating the new merging array!\n";
	}
}

no Moose::Role;
1;
