package MiRNAture::Evaluate;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(load_all_databases concatenate_true_cand cleancmsearch get_family_name format_final_table analyse_involved_species get_specie_name analyse_strategies get_str_name perform_detection_repeated_loci);

use Moose::Role; #Set of common tools that miRNAture uses
use Data::Dumper;
use File::Basename;
use List::Util 'first';
with 'MiRNAture::ToolBox';

my $DIR;
my (%names_r, %names_r_inverse, %bitscores, %lengs, %families_names, %len_other);

=head1 load_all_databases
    Title: load_all_databases
    Usage: load_all_databases(mode, path_data, return_data);
    Function: Index all needed databases by miRNAture. Due different running modes, 
	 this subroutine can load 3 different modes: C<Basic> load only the 
	 scores from RFAM file. C<Additional> are loaded in case new scores 
	 are needed to run the structural validation. C<Joined> loads both score's
	 files. 
    Returns: A set of hashes that saves the specific values for each RFAM family,
	 as: bitscore, length, accession numbers, and ncRNA family names. If 
	 mode is set as B<Additional>, only it is returned the index of the 
	 reported lengths.
=cut 

sub load_all_databases {
	my ($modeT, $datapath, $user) = @_;
	if ($modeT =~ /^Basic$/){ #Only include RFAM
		open my $BASIC, "< $datapath/all_RFAM_scores.txt" or die "Critical error: Does not have the file $datapath/all_RFAM_scores.txt, which provides all scores to miRNAture\n";
		#RF00006    34.00   96  Vault   misc_RNA
		while (<$BASIC>){
			chomp;
			my @fields = split /\s+|\t/, $_;
			$bitscores{$fields[0]} = $fields[1];
			$lengs{$fields[0]} = $fields[2];
			$names_r{$fields[0]} = $fields[3];
			$names_r_inverse{$fields[3]} = $fields[0];
			$families_names{$fields[0]} = $fields[-1];
		}
		close $BASIC;
		return \%bitscores, \%lengs, \%names_r, \%names_r_inverse, \%families_names;
	} elsif ($modeT =~ /^Additional$/){ #Include miRBase and User
		open my $OTHER, "< $datapath/all_other_scores.txt" or die "Critical error: Does not have the file $datapath/all_other_scores.txt, which provides all scores to miRNAture\n";
		#RF00006    34.00   96  Vault   misc_RNA
		while (<$OTHER>){
			chomp;
			my @fields = split /\s+|\t/, $_;
			$bitscores{$fields[0]} = $fields[1];
			$lengs{$fields[0]} = $fields[2];
			$names_r{$fields[0]} = $fields[3];
			$names_r_inverse{$fields[3]} = $fields[0];
			$families_names{$fields[0]} = $fields[-1];
		}
		close $OTHER;
		open my $USER, "< $user/all_user_scores.txt" or die "Critical error: Does not have the file $user/all_other_scores.txt, which provides all scores to miRNAture\n";
		#RF00006    34.00   96  Vault   misc_RNA
		while (<$USER>){
			chomp;
			my @fields = split /\s+|\t/, $_;
			$bitscores{$fields[0]} = $fields[1];
			$lengs{$fields[0]} = $fields[2];
			$names_r{$fields[0]} = $fields[3];
			$names_r_inverse{$fields[3]} = $fields[0];
			$families_names{$fields[0]} = $fields[-1];
		}
		close $USER;
		return \%bitscores, \%lengs, \%names_r, \%names_r_inverse, \%families_names;
	} elsif ($modeT =~ /^mirbase$/){ #Include miRBase and User
		open my $OTHER, "< $datapath/all_other_scores.txt" or die "Critical error: Does not have the file $datapath/all_other_scores.txt, which provides all scores to miRNAture\n";
		#RF00006    34.00   96  Vault   misc_RNA
		while (<$OTHER>){
			chomp;
			my @fields = split /\s+|\t/, $_;
			$bitscores{$fields[0]} = $fields[1];
			$lengs{$fields[0]} = $fields[2];
			$names_r{$fields[0]} = $fields[3];
			$names_r_inverse{$fields[3]} = $fields[0];
			$families_names{$fields[0]} = $fields[-1];
		}
		close $OTHER;
		return \%bitscores, \%lengs, \%names_r, \%names_r_inverse, \%families_names;
	} elsif ($modeT =~ /^Joined$/){ #When both modes are used other + infernal 
		#Clear previously defined hashes
		undef %names_r;
		undef %names_r_inverse;
		undef %bitscores;
		undef %lengs;
		undef %families_names;
		undef %len_other;
		#
		open my $BASIC, "< $datapath/all_RFAM_scores.txt" or die "Critical error: Does not have the file $datapath/all_RFAM_scores.txt, which provides all scores to miRNAture\n";
		#RF00006    34.00   96  Vault   misc_RNA
		while (<$BASIC>){
			chomp;
			my @fields = split /\s+|\t/, $_;
			$bitscores{$fields[0]} = $fields[1];
			$lengs{$fields[0]} = $fields[2];
			$names_r{$fields[0]} = $fields[3];
			$names_r_inverse{$fields[3]} = $fields[0];
			$families_names{$fields[0]} = $fields[-1];
		}
		close $BASIC;
		open my $OTHER, "< $datapath/all_other_scores.txt" or die "Critical error: Does not have the file $datapath/all_other_scores.txt, which provides all scores to miRNAture\n";
		#RF00006    34.00   96  Vault   misc_RNA
		while (<$OTHER>){
			chomp;
			my @fields = split /\s+|\t/, $_;
			$bitscores{$fields[0]} = $fields[1];
			$lengs{$fields[0]} = $fields[2];
			$names_r{$fields[0]} = $fields[3];
			$names_r_inverse{$fields[3]} = $fields[0];
			$families_names{$fields[0]} = $fields[-1];
		}
		close $OTHER;
		open my $USER, "< $user/all_user_scores.txt" or die "Critical error: Does not have the file $user/all_user_scores.txt, which provides all scores to miRNAture\n";
		#RF00006    34.00   96  Vault   misc_RNA
		while (<$USER>){
			chomp;
			my @fields = split /\s+|\t/, $_;
			$bitscores{$fields[0]} = $fields[1];
			$lengs{$fields[0]} = $fields[2];
			$names_r{$fields[0]} = $fields[3];
			$names_r_inverse{$fields[3]} = $fields[0];
			$families_names{$fields[0]} = $fields[-1];
		}
		close $USER;
		return \%bitscores, \%lengs, \%names_r, \%names_r_inverse, \%families_names;
	} elsif ($modeT =~ /^JoinedN$/){ #When both modes are used other + infernal 
		#Clear previously defined hashes
		undef %names_r;
		undef %names_r_inverse;
		undef %bitscores;
		undef %lengs;
		undef %families_names;
		undef %len_other;
		#
		open my $BASIC, "< $datapath/all_RFAM_scores.txt" or die "Critical error: Does not have the file $datapath/all_RFAM_scores.txt, which provides all scores to miRNAture\n";
		#RF00006    34.00   96  Vault   misc_RNA
		while (<$BASIC>){
			chomp;
			my @fields = split /\s+|\t/, $_;
			$bitscores{$fields[0]} = $fields[1];
			$lengs{$fields[0]} = $fields[2];
			$names_r{$fields[0]} = $fields[3];
			$names_r_inverse{$fields[3]} = $fields[0];
			$families_names{$fields[0]} = $fields[-1];
		}
		close $BASIC;
		open my $OTHER, "< $datapath/all_other_scores.txt" or die "Critical error: Does not have the file $datapath/all_other_scores.txt, which provides all scores to miRNAture\n";
		#RF00006    34.00   96  Vault   misc_RNA
		while (<$OTHER>){
			chomp;
			my @fields = split /\s+|\t/, $_;
			$bitscores{$fields[0]} = $fields[1];
			$lengs{$fields[0]} = $fields[2];
			$names_r{$fields[0]} = $fields[3];
			$names_r_inverse{$fields[3]} = $fields[0];
			$families_names{$fields[0]} = $fields[-1];
		}
		close $OTHER;
	       	return \%bitscores, \%lengs, \%names_r, \%names_r_inverse, \%families_names;
	} else {
		print_error("Covariance models are restricted to RFAM, miRBase and/or user created models. Selection of parameters failed");
	}	
}

sub concatenate_true_cand {
	my ($specie, $dir, $files_true, $str) = @_;
	my $output_file;
	if ($str !~ /^HMM$|^INFERNAL$|^Final$|^OTHER_CM$/){ #Str is not Infernal or HMM, only blast
		foreach my $input_file (@$files_true){
			system "cat $dir/$input_file >> $dir/all_RFAM_${specie}_${str}.truetable.temp";
		}
		$output_file = "$dir/all_RFAM_${specie}_${str}.truetable.temp";
	} elsif ($str =~ /^Final$/){
		foreach my $input_file (@$files_true){
			system "cat $input_file >> $dir/all_RFAM_${specie}_${str}.truetable.temp";
		}
		$output_file = "$dir/all_RFAM_${specie}_${str}.truetable.temp";
	} else {  # Here HMM, INFERNAL   
		foreach my $input_file (@$files_true){
			system "cat $dir/$input_file >> $dir/all_RFAM_$specie.truetable.temp";
		}
		$output_file = "$dir/all_RFAM_$specie.truetable.temp";
	}
	my $output_file_final = $output_file;
	$output_file_final =~ s/(.*)(\.temp)/$1/g;
	open my $GENERATED, "< $output_file";
	open my $OUT, "> $output_file_final";
	if (!-e $output_file || -z $output_file) {
		return;
	}
	while (<$GENERATED>){
		chomp;
		# JH126831.1.-.75502.75785.41	-	mir-10	RF00104	cm	1	74	159	220	+	no	1	0.26	0.2	15.5	0.0021	!	-	74
		if ($str =~ /^ALL$|^Final$/){
			print $OUT "$_\n";
		} else { #Here Blast, Infernal and HMM modes
			my @split = split /\s+|\t/, $_;
			my @other = split /\./, $split[0];
			print $OUT "$_\n";
		}
	}
	close $OUT;
	return;
}

sub cleancmsearch {
	my ($tabfile, $defined_nGA, $mode, $bitscores, $lengs, $names_r, $minBitscore) = @_; #Mode=1 with GA score, 2: no GA score
	$tabfile =~ s/(\.\.\/.*\/|\/.*\/|.*\/)(.*)/$1$2/g;
	my $path_data = $1;
	my $filename = $2;
	open my $IN, "< $tabfile" or die; #CM Output;
	open my $OUT1, "> ${path_data}/Final/${filename}.true.table" or die;
	open my $OUT2, "> ${path_data}/Final/${filename}.false.table" or die;
	my (@bits, @bitscore, @truetab, @nonaccepted, @accepted);
	my ($cm, $bitscore, $bitscore78,$max, $max70, $starth, $endh, $trunc);
	while (<$IN>){
		chomp;
		next if ($_ =~ m/^\#/);
		my @valuesT = split /\s+|\t/, $_;
		if ($valuesT[2] =~ /^RF[0-9]+/ || $valuesT[2] =~ /^\-$/ || $valuesT[3] =~ /^\-$/){ #Missing real name in column 2
			if ($valuesT[3] =~ /^\-$/){
				if ($valuesT[2] =~ /^\-$/){
					print_error("Output file from cmsearch has not information about RFAM family");
				} else { #Based on the name on column 2, infer column 3. Then search the correct RFAM name and finally replace it on column 2
					#Specific modifications RFAM
					$valuesT[3] = format_acc_name($valuesT[2]); #Based on column 3, infer column 4 as RF* name
					$valuesT[2] = $$names_r{$valuesT[3]}; 
				}
			} elsif ($valuesT[3] =~ /^RF[0-9]+$/){ #Based on Column 3, give the RFAM name
				$valuesT[2] = $$names_r{$valuesT[3]}; 
			} else {
				$valuesT[2] = $valuesT[3];
			}
		}
		if ($valuesT[5] < $valuesT[6]){
			$starth = $valuesT[5];
			$endh = $valuesT[6];
		} else {
			$starth = $valuesT[6];
			$endh = $valuesT[5];
		}
		my $diffh = abs($endh - $starth) + 1; #Coverage respect to the model
		my $evalue = $valuesT[15]; #E-value reported from Infernal for each record.
		my $bitsc = $valuesT[14]; #BitScore
		my $bitscoreC;
		if ($mode == 1){
			if ($$bitscores{$valuesT[3]} == 0){ #CMs did not reported bitscore threshold
				$bitscoreC = 1;
			} else {
				$bitscoreC = $valuesT[14]/$$bitscores{$valuesT[3]}; #nGA
			}
			if (exists $$lengs{$valuesT[2]}){
				$max = $$lengs{$valuesT[2]};
			} elsif (exists $$lengs{$valuesT[3]}){
				$max = $$lengs{$valuesT[3]};
			}
		} elsif ($mode == 2){
			$bitscoreC = 1; #
			if (exists $$lengs{$valuesT[2]}){
				$max = $$lengs{$valuesT[2]};
			} elsif (exists $$lengs{$valuesT[3]}){
				$max = $$lengs{$valuesT[3]};
			}
		} 
		if (!$max){
			print_error("Seems that the family: $_\n did not have a reference into the scores list, please add to scores files\n");
		}
		my $ln = join "\t", @valuesT;
		if ($evalue > 0.01){ #miRNAture
			#if ($evalue > 100){
			print $OUT2 "$ln\t$max\n";
		} else {
			if ($bitscoreC < $defined_nGA || $bitsc <= $minBitscore){ #Defined log2(N) <= x & nx >= nGA 
				print $OUT2 "$ln\t$max\n";
			} else {
				my $max70 = sprintf("%.1f", ($max * 0.7)); #miRNAture
				#my $max70 = sprintf("%.1f", ($max * 0)); # Default 
				if ($diffh < $max70){ #<<<<<<< Length
					print $OUT2 "$ln\t$max\n";
				} else {
					print $OUT1 "$ln\t$max\n";
				}
			}
		}
	}
	close $IN; close $OUT1; close $OUT2;
	return;
}

sub get_family_name {
	my ($acc, $families_names) = @_;
	my $name;
	if (exists $$families_names{$acc}){
		$name = $$families_names{$acc};
	} else {
		$name = "lncRNA";
		print_error("Your $acc has not description of the ncRNA family!");
	}
	return $name;
}

sub index_new_old_contig_names {
	my $index_genome_file = shift;
	open my $IN, "< $index_genome_file" or die "The genome index file was not created\n";
	my %index;
	while (<$IN>){
		my @all = split /\s+|\t/, $_;
		$index{$all[0]} = $all[1];
	}
	close $IN;
	return \%index;
}

sub format_final_table {
	my ($input_folder, $specie, $families_names, $names_r_inverse, $database_names_contigs) = @_;
	my $file_in = "$input_folder/all_RFAM_${specie}_Final.truetable.joined.table"; 
	my $file_out = "$input_folder/all_RFAM_${specie}_Final.ncRNAs_homology.txt.temp";
	my $index_names = index_new_old_contig_names($database_names_contigs); # tagnameNumb => contigName
	open my $INF, "< $file_in" or die "The file $file_in not exists\n";  #all_RFAM_Dive_Final.truetable.joined.table
	open my $OUTF, "> $file_out" or die;  
	while (<$INF>){
		chomp;
		#print "$_\n";
		my @split = split /\s+|\t/, $_;
		my ($ncRNA_name, $family_ncRNA);
		if ($split[10] !~ /,/){
			$ncRNA_name = $names_r_inverse{$split[10]};
			$family_ncRNA = $families_names{$ncRNA_name};
		} else { #MIPF0000204,MIPF0001319 Different families
			$ncRNA_name = $split[10];
			$family_ncRNA = "miRNA";
		}
		my $involved_species = analyse_involved_species($split[2], $split[-1]);
		my $current_strategies = analyse_strategies($split[-1]);
		#$split[0] =~ s/(.*)(\.\d+)/$1/g;
		if (exists $$index_names{$split[0]}){
			$split[0] = $$index_names{$split[0]};
		} else {
			print_error("The $split[0] reference was not found on the indexed genome\n");
		}
		my $newline = "$split[0]\t$split[1]\t$ncRNA_name\t$split[3]\t$split[4]\t$split[5]\t$split[6]\t$split[7]\t$split[8]\t$split[9]\t$split[10]\t$split[11]\t$split[-1]\t$current_strategies\t$involved_species\t$split[2]\t$family_ncRNA\n";
		print $OUTF "$newline";
	}
	close $INF; close $OUTF;
}

sub perform_detection_repeated_loci {
	my ($final_table, $mode, $threshold_repeat, $selection_number) = @_;
	open my $IN, "< $final_table" or die;
	$final_table =~ s/(.*)(\.temp)$/$1/g;
	open my $TEMP, "> $final_table" or die; #Here the *.db is created
	open my $POTENTIAL, "> ${final_table}-potential.txt" or die; 
	my %database;
	while(<$IN>){
		chomp;
		my @split = split /\s+|\t/, $_;
		if ($split[10] =~ /,/){
			my @fams = split /,/, $split[10];
			foreach my $fam (@fams){
				$split[2] = $fam; #Reeplace for the current one
				$split[10] = $fam; #Reeplace the detected group of families by the current one
				push @{$database{$fam}}, [ @split ];
			}
		} else {
			push @{$database{$split[10]}}, [ @split ]; #Group by family
		}
	}
	foreach my $mirna_acc (sort keys %database){
		my $candidates = $database{$mirna_acc};
		my $number = scalar @$candidates;
		if ($number >= $threshold_repeat){ #Select those families that reported high number loci
			my @sorted_bit = sort { $a->[7] <=> $b->[7] } @$candidates; #Organize by bitscore
			my $selected = select_best_score_candidates(\@sorted_bit, $selection_number); #Select the top X candidates
			my $potential = select_low_score_candidates(\@sorted_bit, $selection_number); #Select the remaining
			foreach my $ln (@$selected){
				print $TEMP "$ln\n";
			}
			foreach my $ln2 (@$potential){
				print $POTENTIAL "$ln2\n";
			}
		} else {
			foreach my $ln (@$candidates){
				my $ln_out;
				for (my $j = 0; $j <= 16; $j++){
					$ln_out .= "$$ln[$j]\t";
				}
				$ln_out =~ s/^(.*)\s+$/$1/g;
				print $TEMP "$ln_out\n";
			}
		}
	}
	close $IN;
	close $TEMP;
	close $POTENTIAL;
	return $final_table;
}

sub select_best_score_candidates {
	my ($sorted, $threshold) = @_;
	my @selected;
	for (my $i=1; $i<=$threshold; $i++){ #Last started on -1
		my $ln;
		for (my $j = 0; $j <= 16; $j++){
			$ln .= "$$sorted[-$i][$j]\t";
		}
		$ln =~ s/^(.*)\s+$/$1/g;
		push @selected, $ln;
	}
	return \@selected;
}

sub select_low_score_candidates {
	my ($sorted, $threshold) = @_;
	my @selected;
	my $num = scalar @$sorted;
	for (my $i=$threshold+1; $i<=$num; $i++){ #Last started on -1
		my $ln;
		for (my $j = 0; $j <= 16; $j++){
			$ln .= "$$sorted[-$i][$j]\t";
		}
		$ln =~ s/^(.*)\s+$/$1/g;
		push @selected, $ln;
	}
	return \@selected;
}

sub analyse_involved_species {
	my ($species, $str) = @_;
	my $set_species;
	if ($species =~ /,/){
		my @all = split /\,/, $species;
		my @temp;
		foreach my $cand (@all){
			my $spe = get_specie_name($cand, $str);
			push @temp, $spe;
		}	
		my @unique = do { my %seen; grep { !$seen{$_}++ } @temp };
		$set_species = join ",", @unique;
	} else {
		$set_species = get_specie_name($species, $str);	
	}
	return $set_species;
}

sub get_specie_name {
	my ($spe, $str) = @_;
	my $specie_out;
	if ($spe =~ /^\w+/ && $str =~ /^0$/){
		$specie_out = "Infernal";
	} elsif ($spe =~ /^\w{2,4}\d+/ && $str =~ /\d+/ && $str !~ /^0$/){
		#} elsif ($spe =~ /^anca|^brfl|^cael|^ciin|^cisa|^dare|^lach|^oidi|^pema|^sce|^xetr/){
		$specie_out = $spe;
		$specie_out =~ s/^(\w{2,4})(\d+|\_.*)/$1/g;
		$specie_out = uc($specie_out);
	} elsif ($spe =~ /^\w+/ && $str !~ /^0$/) {
		$specie_out = "HMM";
	} else {
		die "This combination of parameters are not allowed\n";
	}
	return $specie_out;
}

sub analyse_strategies {
	my $query = shift;
	my $set_str;
	if ($query =~ /,/){
		my @all = split /\,/, $query;
		my @temp = ();
		foreach my $cand (@all){
			my $spe = get_str_name($cand);
			push @temp, $spe;
		}	
		my @unique = do { my %seen; grep { !$seen{$_}++ } @temp };
		$set_str = join ",", @unique;
	} else {
		$set_str = get_str_name($query);	
	}
	return $set_str;
}

sub get_str_name {
	my $query_2 = shift;
	my $q_out;
	if ($query_2 =~ /^0$|^0D$/){
		$q_out = "Infernal";
	} elsif ($query_2 =~ /^HMM/){
		$q_out = "HMM";
	} else { #numbers from 1 to 10
		$q_out = "Blast";
	}
	return $q_out;
}

sub format_acc_name {
	my $name = shift;
	#print $name."\n";
	#model: RF\d.cm
	if ($name =~ m/^RF[0-9]/) {
		if ($name !~ m/RF[0-9]+\.cm/){
			$name =~ s/^(RF[0-9]+)(\.|\_|\-)(.*)$/$1/g;    
		} else {
			$name =~ s/^(RF[0-9]+)(\.cm)$/$1/g;
		}
		return $name; #obtain only the RF0* tag
	} else { #Return the name in column 2
		return $name;
	}
}

no Moose::Role;
1;
