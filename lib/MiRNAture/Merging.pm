package MiRNAture::Merging;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(resolve_mergings direct_flow analyse_sorted compare_coord to_print test_length get_uniq obtain_best_candidate get_data_best filter_candidate fusion_all error push_best push_to_new generate_final_ncRNAs);
use Moose::Role;
use Data::Dumper;
use File::Basename;
use MiRNAture::ToolBox;
use List::Util 'first';
use Math::BigFloat;

my (%rows, %cm_len, %names_r);
my (@temporal, @new_coord, @new_coord_bad, @sorted, @new, @all);
my ($new_coord, $new_coord_bad);
my ($MODE, $tto);
my ($NAMES, $LEN, $OUT, $SPECIAL,$IN);

with 'MiRNAture::Evaluate';

sub load_database_query {
	my ($folder, $pattern, $specie) = @_;
	my %database_names;
	my $IN;
	my @all_database_files = check_folder_files($folder, $pattern); #Files generated with define_final_CMs on all strategies
	foreach my $filesdb (@all_database_files){
		my $str = $filesdb;
		$str =~ s/(${specie}\_)([0-9]+)(\.miRNA\.tab\.db\.location\.database)/$2/g;
		my $complete = "${folder}${filesdb}";
        if(!-e $complete || -z $complete){
           print_error("The database file $complete did not exists"); 
           die;
        }
		open $IN, "< $complete" or die "The database file is corrupted\n";
		while (<$IN>){
			chomp;
			# 0	JH126831.1	anca148202	+	11795	11880	30	75	87
			my $group = (split /\s+|\t/, $_)[0];
			my $query = (split /\s+|\t/, $_)[2];
			$database_names{"$group.$str"} .= "$query,";
		}
	}
	return \%database_names;
}

sub discover_query_sequences {
	my ($group, $str, $database) = @_;
	my $queries;
	if (exists $$database{"$group.$str"}){
		$queries = $$database{"$group.$str"};
		$queries =~ s/^(.*)(\,)$/$1/g;
		return $queries;
	} else {
		print_error("The database of queries has a wrong reference at: $group and $str");
	}
}

sub resolve_mergings {
	my ($specie, $dir, $MODE, $str) = @_;
	my $file;
	my $tto = 0;
	my $database_grouped_queries;
	if ($str =~ /^HMM$/){
		$file = "${dir}/all_RFAM_$specie.truetable.clean"; #Adjusted HMM coordinates
	} elsif ($str =~ /^INFERNAL$|^OTHER_CM$/){
		$file = "${dir}/all_RFAM_$specie.truetable"; #INFERNAL direct coordinates
	} elsif ($str =~ /^ALL$|^COMPLETE$/) {
		my $database_folder = "${dir}/../../";
		my $pattern_file_db = "miRNA\\.tab\\.db\\.location\\.database"; #All database files
		# Load all database files to get query references
		$database_grouped_queries = load_database_query($database_folder, $pattern_file_db, $specie);
		$file = "${dir}/all_RFAM_${specie}_${str}.truetable"; #Adjusted Blast ALL str coordinates
    } elsif ($str =~ /^Final$/) {
		$file = "${dir}/all_RFAM_${specie}_${str}.truetable"; #Final file concatenated
	} elsif ($str =~ /^Merging$/){
		$file = "${dir}/all_RFAM_${specie}_Final.all";
	} elsif ($str =~ /^\d+$/) {
		$file = "${dir}/all_RFAM_${specie}_${str}.truetable.clean"; #Adjusted Blast specific str coordinates
	} else {
		print_error("$str is not recognized");
	}
	if (!-e $file || -z $file) {
		print_result("No valid candidates were found by $str method");
		return;
	}
	open $IN, "< $file" or die;
	open $OUT, "> $file.joined.table";
	open $SPECIAL, "> $file.discarded.table";
	%rows = ();
	while (<$IN>){
		chomp;
		my @database = split /\s+|\t/, $_;
		if ($MODE == 2){
			push @{$rows{"$database[0] $database[2]"}}, [ @database[3,4,5,6,7,8,9,10,11,12,1,2,0] ];
		} elsif ($MODE == 1){
			push @{$rows{"$database[0] $database[2]"}}, [ @database[3,4,5,6,7,8,9,10,11,12,1,2,0] ];
			#push @{ $rows{"$database[0] $database[1]"} } , [ @database[3,4,5,6,7,8,9,10,11,12,2,1,0] ];
			#push @{ $rows{"$database[0] $database[2]"} } , [ @database[3,4,5,6,7,8,9,10,11,12,1,2,0] ]; #Considering strand

		} elsif ($MODE == 4){ #Blast ALL
			# JH126831_0	1_1	-_2	1_3	74_4	75661_5	75722_6	15.5_7	0.0021_8	no_9	mir-10_10	74_11	1_12
			$database[1] = discover_query_sequences($database[1], $database[-1], $database_grouped_queries);
			push @{$rows{"$database[0] $database[2]"}}, [ @database[3,4,5,6,7,8,9,10,11,12,1,2,0] ];
		} elsif ($MODE == 5){
			#scaffold101126-size1290_0 -_1       cin7,cin9_2       1_3       107_4     1006_5    1105_6    46.6_7    2.2e-09_8 no_9      RF00026_10 104_11     3_12
			push @{$rows{"$database[0] $database[1]"}}, [ @database[3,4,5,6,7,8,9,10,11,12,2,1,0] ];
		} elsif ($MODE == 6){
			#scaffold101126-size1290_0 -_1       cin7,cin9_2       1_3       107_4     1006_5    1105_6    46.6_7    2.2e-09_8 no_9      RF00026_10 104_11     3_12
			push @{$rows{"$database[0] $database[1]"}}, [ @database[3,4,5,6,7,8,9,10,11,12,2,1,0,14,15] ];  #At the end, add specie [14], and original name query [15]. 
		} elsif ($MODE == 3){
			if ($database[9] eq "+"){ #Not considering sense
				push @{$rows{"$database[0] $database[9]"}} , [ @database[5,6,7,8,14,15,10,2,-1],$tto,@database[2,9,0] ];
			} else {
				push @{$rows{"$database[0] $database[9]"}} , [ @database[5,6,8,7,14,15,10,2,-1],$tto,@database[2,9,0] ];
			}
		}
	}

	foreach my $key ( sort keys %rows){
		my $val = $rows{$key};
		my $num = scalar @$val;
		direct_flow($num, $val);
		CHANGE:	
		next;
	}
	close $OUT; close $SPECIAL; close $IN;  
	return;
}

sub direct_flow {
	my ($number, $values) = @_;
	if ($number == 0){
		error();
	} elsif ($number == 1 ){
		foreach my $i (${$values}[0]){
			my $fam2 = get_uniq($$i[10]);
			my $fam9 = get_uniq($$i[6]);
			my $fam10 = get_uniq($$i[7]);
			my $fam11 = get_uniq($$i[8]);
			my $fam12 = get_uniq($$i[9]); 
			print $OUT "$$i[12]\t$$i[11]\t$fam2\t$$i[0]\t$$i[1]\t$$i[2]\t$$i[3]\t$$i[4]\t$$i[5]\t$fam9\t$fam10\t$fam11\t$fam12\n";
		}
	} else {
		#my @sorted;
		@sorted = sort { $a->[2] <=> $b->[2] } @$values; #Organize by the first one coordinate
		@all = @sorted;
		analyse_sorted(\@sorted);
	}
}

sub analyse_sorted {
	my $all_sorted = shift;
	@temporal = ();
	my ($S1, $E1, $S2, $E2);
	FLAG:
	push @temporal, $$all_sorted[0]; #Tomar de a pares, para comparar.
	push @temporal, $$all_sorted[1];
	FLAG1:
	my $sense0 = $$all_sorted[0][11];
	my $sense1 = $$all_sorted[1][11];
	my $interq_diff = abs($$all_sorted[1][2] - $$all_sorted[0][3]);
	$S1 = $temporal[0][2];
	$E1 = $temporal[0][3];
	$S2 = $temporal[1][2];
	$E2 = $temporal[1][3];
	my $overlap = compare_coord($S1,$E1,$S2,$E2);
	if ($overlap == 1){ #Coordinates from candidates are overlapping
		my $CM1 = $temporal[0][4];
		my $CM2 = $temporal[1][4];
		my $eval = obtain_best_candidate($CM1, $CM2, 1); # 1 es el mayor
		if ($eval == 1){ #el mejor es el array 1
			($new_coord, $new_coord_bad) = get_data_best($eval, \@temporal, 0);
			@new = (); #Crear un empty arrays!!
		} elsif ($eval == 2){ #Son iguales los CM
			my $EVAL1 = $temporal[0][5];
			my $EVAL2 = $temporal[1][5];
			my $eval = obtain_best_candidate($EVAL1, $EVAL2, -1); # -1 es el menor
			if ($eval == 1){ #el mejor es el array 1
				($new_coord, $new_coord_bad) = get_data_best($eval, \@temporal, 0);	
				@new = (); #Crear un empty arrays!!
			} elsif ($eval == 2){ #Son iguales los evalue, mirar longitud
				my $long1 = abs($temporal[0][3]-$temporal[0][2]);
				my $long2 = abs($temporal[1][3]-$temporal[1][2]);
				my $eval = obtain_best_candidate($long1, $long2, 1); # 1 es el menor
				if ($eval == 1){ #el mejor es el array 1, por longitud
					($new_coord, $new_coord_bad) = get_data_best($eval, \@temporal, 0);
					@new = (); #Crear un empty arrays!!
				} elsif ($eval == 2){ #Son iguales, comprobar si son iguales los ACCs CMs
					my $rfamname1 = $temporal[0][7];
					my $rfamname2 =	$temporal[1][7];
					if ($rfamname1 eq $rfamname2){
						$new_coord = fusion_all($temporal[0], $temporal[1], 0);
						$new_coord_bad = $new_coord;
						@new = ();
					} else {
						$new_coord = fusion_all($temporal[0], $temporal[1], 1);
						$new_coord_bad = $new_coord;
						@new = ();
					}
				} elsif ($eval == 0){ #el mejor es el array 0
					($new_coord, $new_coord_bad) = get_data_best($eval, \@temporal, 1);
					@new = (); #Crear un empty arrays!!
				} else {
					# Here outside rules
					error();
				}
			} elsif ($eval == 0){ #el mejor es el array 0
				($new_coord, $new_coord_bad) = get_data_best($eval, \@temporal, 1);
				@new = (); #Crear un empty arrays!!
			} else { 
				error();
			}
		} elsif ($eval == 0){ #el mejor es el array 0
			($new_coord, $new_coord_bad) = get_data_best($eval, \@temporal, 1);
			@new = (); #Crear un empty arrays!!
		} else {
			error();
		}
		# With the build $new_coord array with the merged/selected
		# values, next evaluate it
		filter_candidate($new_coord, 1);	
	} else { #There is not overlap, so print the first if exists x > 2 elements
		filter_candidate($temporal[0], 2);
	}
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

sub to_print {
	my $complete = $_[0];
	my $out = $_[1];
	if ($out == 1){ #As an isolated candidate, just print
		my $fam2 = get_uniq($$complete[10]);
		my $fam9 = get_uniq($$complete[6]);
		my $fam10 = get_uniq($$complete[7]);
		my $fam11 = get_uniq($$complete[8]);
		my $fam12 = get_uniq($$complete[9]);
		print $OUT "$$complete[12]\t$$complete[11]\t$fam2\t$$complete[0]\t$$complete[1]\t$$complete[2]\t$$complete[3]\t$$complete[4]\t$$complete[5]\t$fam9\t$fam10\t$fam11\t$fam12\n";
	} else {
		my $fam2 = get_uniq($$complete[10]);
		my $fam9 = get_uniq($$complete[6]);
		my $fam10 = get_uniq($$complete[7]);
		my $fam11 = get_uniq($$complete[8]);
		my $fam12 = get_uniq($$complete[9]);
		print $SPECIAL "$$complete[12]\t$$complete[11]\t$fam2\t$$complete[0]\t$$complete[1]\t$$complete[2]\t$$complete[3]\t$$complete[4]\t$$complete[5]\t$fam9\t$fam10\t$fam11\t$fam12\tNoCoverage\n";
	}
}

sub test_length {
	my $array = $_[0]; # Test the length of the big array
	my $numb = $_[1];
	my $nelements = scalar @$array;
	if ($nelements == 0){
		goto CHANGE;	
	} elsif ($nelements == 1){
		my $diff = ((abs($$array[0][1] - $$array[0][0]) * 100)/$$array[0][8]);
		my $defined_limit = 70; #at least 70% del CM
		if($diff >= $defined_limit){ #Succesfull coverage
			foreach my $i (${$array}[0]){
				my $fam2 = get_uniq($$i[10]);
				my $fam9 = get_uniq($$i[6]);
				my $fam10 = get_uniq($$i[7]);
				my $fam11 = get_uniq($$i[8]);
				my $fam12 = get_uniq($$i[9]);
				print $OUT "$$i[12]\t$$i[11]\t$fam2\t$$i[0]\t$$i[1]\t$$i[2]\t$$i[3]\t$$i[4]\t$$i[5]\t$fam9\t$fam10\t$fam11\t$fam12\n";
				@temporal = ();
				#last; #END ALL
			}
		} else { #Here print those ones with insufficient coverate respect to CM length
			foreach my $i (${$array}[0]){
				my $fam2 = get_uniq($$i[10]);
				my $fam9 = get_uniq($$i[6]);
				my $fam10 = get_uniq($$i[7]);
				my $fam11 = get_uniq($$i[8]);
				my $fam12 = get_uniq($$i[9]);
				print $SPECIAL "$$i[12]\t$$i[11]\t$fam2\t$$i[0]\t$$i[1]\t$$i[2]\t$$i[3]\t$$i[4]\t$$i[5]\t$fam9\t$fam10\t$fam11\t$fam12\n";
				@temporal = ();
			}
			#END ALL with no final candidate
		}
	} else {
		if ($numb == 1){
			goto FLAG;
		} else {
			goto FLAG1;
		}
	}
}

sub get_uniq {
	my $str = $_[0];
	my @data = split /\,/, $str;
	my @unique = do { my %seen; grep { !$seen{$_}++ } @data };
	my $new_str = join ",", (sort @unique);
	return $new_str;
}

sub obtain_best_candidate {
	my ($param1, $param2, $dir) = @_; #Parametro1, parametro2, (best_major || best_menor)
	my $best;
	if($dir == 1){ #El mejor es el mayor
		if($param1 > $param2){
			$best = 0;
		} elsif ($param1 == $param2){
			$best = 2;
		} else {
			$best = 1;
		}
	} elsif ($dir == -1){ #El mejor es el menor
		if($param1 > $param2){
			$best = 1;
		} elsif ($param1 == $param2){
			$best = 2;
		} else {
			$best = 0;
		}
	} elsif ($dir == 0){ #Los mejores son los iguales
		if($param1 > $param2){
			$best = 4;
		} elsif ($param1 == $param2){
			$best = 2;
		} else {
			$best = 4;
		}
	} else {
		&error();
	}
	return $best;
}

sub get_data_best {
	my $nArray = $_[0];
	my $todo = $_[1];
	my $x_array = $_[2];
	my (@new_true, @new_bad);
	my $len = (scalar @{ $$todo[0] }) - 1;
	for (my $j = 0; $j <= $len; $j++){
		push @new_true, $$todo[$nArray][$j];
	}
	for (my $j = 0; $j <= $len; $j++){
		push @new_bad, $$todo[$x_array][$j];
	}
	return (\@new_true, \@new_bad);		
}

sub filter_candidate {
	my $data = $_[0];
	my $mode = $_[1];
	if ($$data[8] == 0){ #Here the length of the CM is 0, die and correct in CM scores file
		print_error("Seems that $$data[-1] has a DB size == 0, please correct it on the all_other_scores.txt");
	}
	my $diff = ((abs($$data[1] - $$data[0]) * 100)/$$data[8]);
	my $defined_limit = 70; #at least 70% del CM
	if ($mode == 2){ # Here only the first candidate is evaluated
		if($diff >= $defined_limit){
			to_print($temporal[0], 1);
			@temporal = ();
			shift @all;
			shift @sorted;		
			test_length(\@all, 1); #number of elements on temporal
		} else { # The fusion of two candidates is evaluated
			#Not printing.
			@temporal = ();
			shift @all;
			shift @sorted;		
			test_length(\@all, 1); #number of elements on temporal
		}
	} else {
		if($diff >= $defined_limit){
			@temporal = ();
			shift @sorted; #Remove [0]
			shift @sorted; #Remove [1]
			shift @all; #Remove [0]
			shift @all; #Remove [1]
			unshift @sorted, $new_coord; 
			unshift @all, $new_coord;
			$new_coord = ();
			$new_coord_bad = ();
			test_length(\@all, 1);
		} else {
			to_print($data, 0); # Print to SPECIAL file
			@temporal = ();
			$new_coord = (); ##
			shift @sorted; #Remove [0]
			shift @sorted; #Remove [1]
			unshift @sorted, $new_coord_bad;
			shift @all; #Remove [0]
			shift @all; #Remove [1]
			unshift @all, $new_coord_bad;
			$new_coord_bad = ();
			test_length(\@all, 1);
		}
	}
}

sub fusion_all {
	#my (@cand1, @cand2, $mode) = @_;
	my @cand1 = $_[0];
	my @cand2 = $_[1];
	my $mode = $_[2];
	my $new_all;
	my $new_line;
	my $len = (scalar @{$cand1[0]}) - 1;
	if($mode == 0){ #CM names are equal, so merge with greater values
		for (my $a=0;$a <= $len; $a++){
			if ($a == 0 || $a == 2 || $a == 5){ #El mejor es el menor
				push_best($cand1[0][$a], $cand2[0][$a], -1);
			} elsif ($a == 1 || $a == 3 || $a == 4 || $a == 8){ #El mejor es el mayor
				push_best($cand1[0][$a], $cand2[0][$a], 1);
			} elsif ($a == 6 || $a == 7 || $a == 9 || $a == 10 || $a == 11 || $a == 12){ #Text comparison, merge
				push_best($cand1[0][$a], $cand2[0][$a], 0);
			} else { #text comparison, merge by comma
				error();
			}
		}	
	} elsif ($mode == 1){ #CM names are diferent merge with greater values and names with comma
		for (my $a=0;$a <= $len; $a++){
			if ($a == 0 || $a == 2 || $a == 5){ #El mejor es el menor
				push_best($cand1[0][$a], $cand2[0][$a], -1);
			} elsif ($a == 1 || $a == 3 || $a == 4 || $a == 8){ #El mejor es el mayor
				push_best($cand1[0][$a], $cand2[0][$a], 1);
			} elsif ($a == 6 || $a == 7  || $a == 9 || $a == 10 || $a == 11 || $a == 12 || $a == 13){ #Text comparison, merge
				push_best($cand1[0][$a], $cand2[0][$a], 0);
			} else { #text comparison, merge by comma
				error();
			}
		}	
	}
	#Here @new was populated, then return it
	my @all_new = @new;
	return \@all_new;
}

sub error{
	print_error("There is something really wrong!");
	die;
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
		error();
	}
	#return $best;
	push_to_new($best, $test1, $test2);
	return;
}

sub push_to_new {
	# Here I fill the global @new with the fusion
	my ($cand, $text1, $text2) = @_;
	#my @temp;
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
		error();
	}
	return;
}

sub generate_final_ncRNAs {
	my ($blast_file, $hmm_file, $infernal_file, $other_file, $out_folder, $specie) = @_;
	my (@all_files, @all_filesT);
	push @all_filesT, $blast_file;
	push @all_filesT, $hmm_file;
	push @all_filesT, $infernal_file;
	push @all_filesT, $other_file;
	foreach my $file (@all_filesT){
		if (-z $file || !-e $file){
			print_result("$file is empty or is missing!");
		} else {
			push @all_files, $file;
		}
	}
	concatenate_true_cand($specie, $out_folder, \@all_files, "Final"); #Concantenate all true
	resolve_mergings($specie, $out_folder, "5", "Final");
	return;
}

no Moose::Role;
1;
