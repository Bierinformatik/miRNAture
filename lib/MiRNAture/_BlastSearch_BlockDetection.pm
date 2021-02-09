package MiRNAture::_BlastSearch_BlockDetection;

use Moose::Role;
use Data::Dumper;

sub generate_blocks {
	my $shift = shift;
	my $input_path_location_file = shift;
	my $group_id = 0;
	my $sumGroup = 0;
	my %base;
	open my $IN, "< $input_path_location_file" or die; #Here is the location files by Str
	open my $F_OUT, "> ${input_path_location_file}.blocks.coord" or die; #Location file
	open my $F_OUT2, "> ${input_path_location_file}.database" or die;
	while (<$IN>){
		chomp;
		my @fields = split /\s+|\t/, $_;
		push @{ $base{"$fields[1]\t$fields[2]"} }, [ @fields[3,4,5,6,7,0,1,2] ];
	}
	foreach my $key ( sort keys %base){
		my $val = $base{$key};
		my $num = scalar @$val;
		$group_id = direct_flow($num, $val, $F_OUT, $F_OUT2, $group_id);
	}
	close $IN;
	close $F_OUT; close $F_OUT2;
	return;
}

sub direct_flow {
	my ($number, $values, $F_OUT, $F_OUT2, $group_id) = @_;
	if ($number == 0){
		die "Exists a fatal error related with number of values and keys from hash of candidates!\n";
	} elsif ($number == 1 ){
		#Check if have the coverage!
		$group_id++;
		foreach my $i (${$values}[0]){
			#sce1 dvexscf18546 +	 3359	 3424	 8	 73	 82
			print $F_OUT2 "${group_id}A\t$$i[6]\t$$i[5]\t$$i[7]\t$$i[0]\t$$i[1]\t$$i[2]\t$$i[3]\t$$i[4]\n"; # Database 
			print $F_OUT "${group_id}A\t$$i[6]\t$$i[7]\t$$i[0]\t$$i[1]\n"; # Retrieve seq
		}
	} else { #More than 1
		my @sorted = sort { $a->[0] <=> $b->[0] } @$values; #Organize by the first one coordinate
		my ($grouped, $database, $group_f) = analyse_sorted(\@sorted, $group_id);
		# Create the file with the blocks positions and fasta sequences
		# to be tested.
		foreach my $group (sort {$a <=> $b} keys %{$grouped}){
			foreach my $chr (keys %{$$grouped{$group}}){
				foreach my $strand (keys %{$$grouped{$group}{$chr}}){
					#>JH126831.1.+.11776.11879.pema9515754.39.75.84
					print $F_OUT "$group\t$chr\t$strand\t$$grouped{$group}{$chr}{$strand}{'start'}\t$$grouped{$group}{$chr}{$strand}{'end'}\n";
				}
			}
		}
		foreach my $grouped (sort keys %{$database}){
			my $all = $$database{$grouped};
			foreach my $candidate (@$all){
				my $new_line = "$$candidate[6]\t$$candidate[5]\t$$candidate[7]\t$$candidate[0]\t$$candidate[1]\t$$candidate[2]\t$$candidate[3]\t$$candidate[4]";
				print $F_OUT2 "$grouped\t$new_line\n";
			}
		}
		$group_id = $group_f;
	}
	return $group_id;
}

sub correct_coordenates {
	my ($coordinate, $mappingq, $lengthq, $mode) = @_;
	my $new;
	if ($mode == -1) { # This is the start coordinate S1 - Sq .
		$new = ($coordinate - $mappingq) + 1;
		if ($new < 1 ){
			$new = 1;
		}
	} elsif ($mode == 1){ # This is the end coordinate E2 + (Eq - Lenf)
		my $dist = ($lengthq - $mappingq) + 1;
		$new = $coordinate + $dist;
	} else {
		die "Error mode\n";
	}
	return $new;
}

sub analyse_sorted {
	my ($all_sorted, $group_id) = @_;
	my (%data, %database);
	$group_id++; #Not collide to the last group
	foreach (@$all_sorted){
		$$_[0] = correct_coordenates($$_[0], $$_[2], $$_[4], -1); # Start and startq, mode -1
		$$_[1] = correct_coordenates($$_[1], $$_[3], $$_[4], 1); # End and endq mode 1
		if (exists $data{$group_id} && exists $data{$group_id}{$$_[6]} && exists $data{$group_id}{$$_[6]}{$$_[7]}){
			my ($S1, $E1, $S2, $E2);
			$S1 = $data{$group_id}{$$_[6]}{$$_[7]}{'start'};
			$E1 = $data{$group_id}{$$_[6]}{$$_[7]}{'end'};
			$S2 = $$_[0];
			$E2 = $$_[1];
			my $overlap = compare_coord($S1,$E1,$S2,$E2); #Check if ranges overlaped
			if ($overlap == 1){ # If overlapping update main
				my $min = push_best($S1, $S2, -1);
				my $max = push_best($E1, $E2, 1);
				my $distance = abs($max - $min) + 1;
				push @{ $database{$group_id}}, $_;
				$data{$group_id}{$$_[6]}{$$_[7]}{'start'} = $min;
				$data{$group_id}{$$_[6]}{$$_[7]}{'end'} = $max;
				$data{$group_id}{$$_[6]}{$$_[7]}{'distance'} = $distance;
				$data{$group_id}{$$_[6]}{$$_[7]}{'n'} += 1;
			} else {
				$group_id += 1;
				$data{$group_id}{$$_[6]}{$$_[7]}{'start'} = $$_[0];
				$data{$group_id}{$$_[6]}{$$_[7]}{'end'} = $$_[1];
				my $distance = abs($$_[1] - $$_[0]) + 1;
				$data{$group_id}{$$_[6]}{$$_[7]}{'distance'} = $distance;
				$data{$group_id}{$$_[6]}{$$_[7]}{'n'} += 1;
				push @{ $database{$group_id}}, $_;
			}
		} else {
			$data{$group_id}{$$_[6]}{$$_[7]}{'start'} = $$_[0];
			$data{$group_id}{$$_[6]}{$$_[7]}{'end'} = $$_[1];
			my $distance = abs($$_[1] - $$_[0]) + 1;
			$data{$group_id}{$$_[6]}{$$_[7]}{'distance'} = $distance;
			$data{$group_id}{$$_[6]}{$$_[7]}{'n'} += 1;
			push @{ $database{$group_id}}, $_;
		}
	}
	return (\%data, \%database, $group_id);
}

sub fusion_all {
	#my (@cand1, @cand2, $mode) = @_;
	my @cand1 = $_[0];
	my @cand2 = $_[1];
	my $mode = $_[2];
	my $new_line;
	my @new;
	if ($mode == 1){ #Merge with greater values and names with comma
		for (my $a=0; $a <= 7; $a++){ # Walking throught the array
			if ($a == 0 || $a == 2){ #Select the lowest, is the best.
				my $best =  push_best($cand1[0][$a], $cand2[0][$a], -1);
				push @new, $best;
			} elsif ($a == 1 || $a == 3 || $a == 4){ #Greater value is the best
				my $best = push_best($cand1[0][$a], $cand2[0][$a], 1);
				push @new, $best;
			} elsif ($a == 5 || $a == 6 || $a == 7){ #Text comparison, merge by comma if different.
				my $best = push_best($cand1[0][$a], $cand2[0][$a], 0);
				push @new, $best;
			} else { 
				die "It contains more fields to merge, check the input format\n";
			}
		}	
	}
	return @new;
}

sub push_best {
	my ($test1, $test2, $ind) = @_;
	my $best;
	if ($ind == -1){ #menor
		if ($test1 <= $test2){
			$best = $test1;
		} else {
			$best = $test2;
		}
	} elsif ($ind == 1){ #mayor
		if ($test1 <= $test2){
			$best = $test2;
		} else {
			$best = $test1;
		}
	} elsif ($ind == 0){ #text
		if ($test1 eq $test2){
			$best = $test1; #$test1;	
		} else {
			$best = "$test1,$test2"; #"($test1,$test2)";
		}
	} else {
		die "Wrong definitions on merging rules!\n";
	}
	return $best;
}

sub compare_coord {
	my ($s1,$e1,$s2,$e2) = @_;
	my $response;
	if ( ($s1 <= $s2 && $e1 >= $e2) || ($s1 <= $s2 && $s2 <= $e1 && $e1 <= $e2) || ($s1 >= $s2 && $e2 >= $s1 && $e2 <= $e1) || ($s1 >= $s2 && $s1 <= $e2 && $e1 <= $e2) ){
		$response = 1;
	} else {
		$response = 0;
	}
	return $response;
}

no Moose::Role;
1;
