package miRNAnchor::Validate;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(process_query_fasta_specific generate_index_specific_ids obtain_flip_sequences create_environment_detailed copy_rfam_files_group create_genome_file remove_flip_sequences_group process_query_fasta_group validate_mirnas_group print_results clean_redundancy_file create_environment_detailed_subset copy_rfam_files_group_subset create_genome_file_subset process_query_fasta_group_subset clean_redundancy_file_subset copy_best_element_subset);

use Moose::Role;
use Data::Dumper;
use File::Copy;
use lib 'Modules';
use lib 'Modules/JSON-MaybeXS-1.004000/lib';
use JSON;
use miRNAnchor::Check;

with 'miRNAnchor::Tools';

sub process_query_fasta_specific {
	my ($fasta_query_path, $family, $working_path, $tag) = @_;
	my $fasta_query = (check_folder_files($fasta_query_path, ".*\_$tag.$family\.*\.fasta"))[0]; 	
	#>aeae-RF00103-856 B1566465451 Aedes aegypti RF00103 stem-loop
	if (!-e $fasta_query_path."/".$fasta_query && -z $fasta_query_path."/".$fasta_query) {
		print_error("The input fasta file ".$fasta_query_path."/".$fasta_query." is empty or not exists!");
	}
	open (my $IN, "<", $fasta_query_path."/".$fasta_query) or die "The fasta file: $! is not valid!\n";
	##Modify this list if we want to subset the group of sequences
	#my $specific_ids = "/scr/k70san/bogota_unal/miRNAsTunicata/miRNAture/Code/RFAM_Processing/Code_Final/ValidateNoFlip/completeLIST";
	my (%candidates, %fasta);
	## Activate to list specific candidates
	#my $identification = generate_index_specific_ids($specific_ids);
	my $head;
	while (<$IN>){
		chomp;
		next if $_ =~ /^$/;
		if ($_ =~ /^\>/){
			$head = $_;
			my $id = (split /\s+|\t/, $_)[1];
			## Identify specific list of candidates
			#if (exists $$identification{$id}){
			#	$candidates{$id} = $_;
			#}
		} else {
			$_ =~ tr/Tt/Uu/;
			push @{$fasta{$head}}, $_;
		}
	}
	close $IN; 
	#return (\%candidates, \%fasta);
	return (\%fasta);
}

sub generate_index_specific_ids {
	my $list_file = shift;
	open my $LIST, "< $list_file" or die "Please provide a specific list to process\n";
	my %all;
	while (<$LIST>){
		chomp;
		#RFAM Specific_code
		my $code = (split /\s+|\t/, $_)[1];
		$all{$code} = 1;
	}
	return \%all;
}

sub obtain_flip_sequences {
	my $json_result = shift;
	my %flip_candidates;
	if (-e $json_result && !-z $json_result){
		open my $IN, "< $json_result" or die;
		while (<$IN>){
			chomp;
			$_ =~ s/^\"(.*)\"$/$1/g;
			$_ =~ s/\\//g;
			my $test = decode_json $_;
			my $fippledNumber = $$test{'Flipped(changed) precursors'}{Number};
			if ($fippledNumber > 0){
				my @candidates = $$test{'Flipped(changed) precursors'}{'IDs'};
				foreach my $elements (@{$candidates[0]}){
					$flip_candidates{$elements} = 1;
				}
			}	
		return \%flip_candidates;
		}
	} else {
		return \%flip_candidates;
	}
}

=head1 create_folder_environment_detailed
	Title: create_folder_environment_detailed
	Usage: create_folder_environment_detailed(OUT_FOLDER, )
	Function: Generates a complete file tree required to organize all input/output data from MIRFix.
	Returns: The complete subtree of folders, organized by accession number of the target miRNAs. 
		This example of distribution is noted in the following example:
		ACC_FAMILY_RFAM/
			|-CODE_SPECIFIC_FAMILY
				|-BaseFiles
				|-Families
				|-Output
		
=cut 

sub create_environment_detailed {
	#Test one to one each sequence against the RFAM corrected sequences
	my ($out_path, $family, $code, $address, $tag) = @_;
	create_folder("$out_path/$address");
	create_folder("$out_path/$address/$tag");
    create_folder("$out_path/$address/$tag/$family");
	create_folder("$out_path/$address/$tag/$family/$code");
	create_folder_environment("$out_path/$address/$tag","Validation", "Validation", "$family/$code", "Complete"); #Specific environment from mirfix;
	return;
}

sub create_environment_detailed_subset {
	#Test one to one each sequence against the RFAM corrected sequences
	my ($out_path, $family, $code, $address, $subset, $tag) = @_;
	create_folder("$out_path/$address");
    create_folder("$out_path/$address/$tag");
	create_folder("$out_path/$address/$tag/$family");
    create_folder("$out_path/$address/$tag/$family/$subset");
	create_folder("$out_path/$address/$tag/$family/$subset/$code");
	create_folder_environment("$out_path/$address/$tag","Subset", "Subset", "$family/$subset/$code", "Complete"); #Specific environment from mirfix;
	return;
}

=head1 copy_rfam_files_group
	Title: copy_rfam_files_group 
	Usage: copy_rfam_files_group(OUT_FOLDER, )
	Function: Copy all input files required by MIRfix.
    Returns: All required files into created folders		
=cut 

sub copy_rfam_files_group {
	#Here input files are copied and the names modified
	my ($working_path, $family, $old_path, $code, $address, $tag) = @_;
	my $new_path = "$working_path/$address/$tag/$family/$code";
	copy_files("$old_path/BaseFiles/${family}_list.txt","$new_path/BaseFiles"); #list
	copy_files("$old_path/${family}_mapping_updated.txt","$new_path/BaseFiles"); #mapping updated
	copy_files("$new_path/BaseFiles/${family}.mapping_original.txt","$new_path/BaseFiles"); #mapping exp
	copy_files("$old_path/${family}_mature_updated.fa","$new_path/BaseFiles"); #mature updated
	copy_files("$new_path/BaseFiles/${family}_mature_original.fa","$new_path/BaseFiles"); #mature exp
	copy_files("$old_path/Output/$family.out/${family}-Final.fasta","$new_path/Families"); #Updated fasta hairpin generated by MIRFix
	return;
}

sub copy_rfam_files_group_subset {
	#Here input files are copied and the names modified
	my ($working_path, $family, $old_path, $code, $address, $subset, $tag) = @_;
	my $new_path = "$working_path/$address/$tag/$family/$subset/$code";
    my $family2 = ${subset};
    my $number = $subset;
    $number =~ s/(RF.*)([0-9]+)/$2/g;
	copy_files("$old_path/BaseFiles/${family2}_list.txt","$new_path/BaseFiles"); #list
	copy_files("$old_path/${family2}_mapping_updated.txt","$new_path/BaseFiles"); #mapping updated
	copy_files("$new_path/BaseFiles/${family2}.mapping_original.txt","$new_path/BaseFiles"); #mapping exp
	copy_files("$old_path/${family2}_mature_updated.fa","$new_path/BaseFiles"); #mature updated
	copy_files("$new_path/BaseFiles/${family2}_mature_original.fa","$new_path/BaseFiles"); #mature exp
	copy_files("$old_path/Output_${number}/$family.out/${family}-Final.fasta","$new_path/Families"); #Updated fasta hairpin generated by MIRFix
    system "mv $new_path/Families/${family}.fa $new_path/Families/$family2.fa";
	return;
}

sub create_genome_file {
	my ($working_path, $family, $old_path, $code, $address, $genome_location, $tag) = @_;
	my $new_path = "$working_path/$address/$tag/$family/$code";
	open (my $OUT, ">", "$new_path/BaseFiles/${family}_genomes_list.txt") or die "Can't create mapping file\n";
	open (my $FASTA, "<", "$old_path/Output/$family.out/${family}-Final.fasta"); #Final Fasta file
	my @genome;
	open (my $GENOME, "<", $genome_location) or die "Genome file couldn't open\n";
	while(<$GENOME>){
		chomp;
		push @genome, $_;
	}
	close $GENOME;
	while (<$FASTA>){
		if ($_ =~ /^>/){
			my $specie_code = (split /\s+|\t/, $_)[2];
			create_genome_list_mirfix($new_path, "BaseFiles", $family, \@genome, $specie_code);		
		} else {
			next;
		}
	}
	close $FASTA;
	close $OUT;
	return;
}

sub create_genome_file_subset {
	my ($working_path, $family, $old_path, $code, $address, $genome_location, $subset, $tag) = @_;
	my $new_path = "$working_path/$address/$tag/$family/$subset/$code";
    my $family2 = $subset;
    my $num = $subset;
    $num =~ s/(RF.*)([0-9]+)/$2/g;
	open (my $OUT, ">", "$new_path/BaseFiles/${family2}_genomes_list.txt") or die "Can't create mapping file\n";
	open (my $FASTA, "<", "$old_path/Output_${num}/$family.out/${family}-Final.fasta") or die; #Final Fasta file
	my @genome;
	open (my $GENOME, "<", $genome_location) or die "Genome file couldn't open\n";
	while(<$GENOME>){
		chomp;
		push @genome, $_;
	}
	close $GENOME;
	while (<$FASTA>){
		if ($_ =~ /^>/){
			my $specie_code = (split /\s+|\t/, $_)[2];
			create_genome_list_mirfix($new_path, "BaseFiles", $family2, \@genome, $specie_code);		
		} else {
			next;
		}
	}
	close $FASTA;
	close $OUT;
	return;
}

sub remove_flip_sequences_group {
	#Remove possible flipped sequences from final hairpin fasta file
	my ($working_path, $rfamily, $pathbase, $flipped_seq, $code, $address) = @_;
	my $new_path = "$working_path/$address/$rfamily/$code";
	#Modify the final fasta file:
	my $finalFasta = "$new_path/Families/$rfamily.fa";
	my %finalFasta;
	if (-e $finalFasta && !-z $finalFasta ){
		open my $FASTA, "< $finalFasta" or die "Updated fasta file is corrupted!\n";
		open my $OUT, "> $finalFasta.temp" or die "Temporal file couldn't be created\n";
		my $head;
		while (<$FASTA>){
			chomp;
			if ($_ =~ /^\>/){
				$head = $_;
			} else {
				push @{ $finalFasta{$head} }, $_;
			} 
		}
		close $FASTA;
		foreach my $head (sort keys %finalFasta){
			my $code = (split /\s+|\t/, $head)[1];
			if (exists $$flipped_seq{$code}){
				next;
			} else {
				print $OUT "$head\n";
				my $seqs = $finalFasta{$head};
				foreach my $ln (@$seqs){
					print $OUT "$ln\n";
				}
			}	
		}	
		close $OUT;
		system "mv $finalFasta.temp $finalFasta";
	} else {
		return;
	}
	return;
}

sub process_query_fasta_group {
	my ($family, $working_path, $code, $header, $seqs, $address, $tag) = @_;
	my $new_path = "$working_path/$address/$tag/$family/$code/Families/$family.fa";
	open my $FASTA_BASE_RFAM, ">> $new_path" or die "Not possible to overwrite $!\n";
	if (exists $$seqs{$header}){
		print $FASTA_BASE_RFAM "$header\n";
		my $all = $$seqs{$header};
		foreach my $line (@$all){
			print $FASTA_BASE_RFAM "$line";
		}
		print $FASTA_BASE_RFAM "\n";
	}
	close $FASTA_BASE_RFAM;
	return;
}

sub process_query_fasta_group_subset {
	my ($family, $working_path, $code, $header, $seqs, $address, $subset,$tag) = @_;
	my $new_path = "$working_path/$address/$tag/$family/$subset/$code/Families/$subset.fa";
	open my $FASTA_BASE_RFAM, ">> $new_path" or die "Not possible to overwrite $!\n";
	if (exists $$seqs{$header}){
		print $FASTA_BASE_RFAM "$header\n";
		my $all = $$seqs{$header};
		foreach my $line (@$all){
			print $FASTA_BASE_RFAM "$line";
		}
		print $FASTA_BASE_RFAM "\n";
	}
	close $FASTA_BASE_RFAM;
	return;
}

sub validate_mirnas_group {
	my ($working_path, $rfam_family, $candidates, $code, $mode, $subset) = @_;
	my $final_file; 
    if ($mode eq "Subset"){
        $final_file = "$working_path/$rfam_family/$subset/$code/Output/$subset.out/$subset-FinalCoor.txt";
    } elsif ($mode eq "Normal"){
        $final_file = "$working_path/$rfam_family/$code/Output/$rfam_family.out/$rfam_family-FinalCoor.txt";
    }
	my $sanity = check_if_file_exists($final_file);
	my (%db, %fasta_info);
	my $result;
	if ($sanity == 1){
		open my $FINAL, "< $final_file" or die "Look up for the final file\n"; 
		#Retrieve info from outfile mirfix
		while (<$FINAL>){
			chomp;
			my $id = (split /\s+|\t/, $_)[0];
			$db{$id} = $_;	
		}
		# Retrieve info from original fasta
		foreach my $id_fasta (sort keys %{$candidates}){
			my $code_id = (split /\s+|\t/, $id_fasta)[1];
			$fasta_info{$code_id} = $id_fasta;
		}
		if (exists $db{$code} && $fasta_info{$code} && $db{$code} !~ m/\s-1\s/){
			$result = "Accepted\t$code\t$fasta_info{$code}\t$db{$code}\t$rfam_family\t$subset";	
		} else {
			#Here discart those ones that reported -1 numbers
			$result = "Discarded\t$code\t$fasta_info{$code}\tNA\t$rfam_family\t$subset"; 
		}
	} else {
        	$result = "NA";
        	return $result;
	}
	return $result;
}

sub copy_sto {
	my ($old, $new) = @_;
	copy($old, $new); #RF000*_X -> RF000*
	return;
}

sub rename_files {
    my ($working_path, $family, $best) = @_;
    #Copied Folder
    rename("$working_path/Output/${family}_${best}.out", "$working_path/Output/${family}.out");
    #Final File
    rename("$working_path/Output/${family}.out/${family}_${best}-FinalCoor.txt", "$working_path/Output/${family}.out/${family}-FinalCoor.txt");
    #Sto File
    rename("$working_path/Output/${family}.out/${family}_${best}.stk", "$working_path/Output/${family}.out/${family}.stk");
    #Final Fasta
    rename("$working_path/Output/${family}.out/${family}_${best}-Final.fasta", "$working_path/Output/${family}.out/${family}-Final.fasta");
    #Final Align Dialign
    rename("$working_path/Output/${family}.out/${family}_${best}-Final.fa", "$working_path/Output/${family}.out/${family}-Final.fa");
    #Json file
    rename("$working_path/Output/${family}.out/${family}_${best}.json", "$working_path/Output/${family}.out/${family}.json");
    #Summary file
    rename("$working_path/Output/${family}.out/${family}_${best}-summ.txt", "$working_path/Output/${family}.out/${family}-summ.txt");
    my $corrected_sto = "$working_path/Output/${family}.out/${family}_${best}Corrected.stk";
    my $corrected_fasta = "$working_path/Output/${family}.out/${family}_${best}-FinalCorrected.fasta";
    if (-e $corrected_sto && !-z $corrected_sto){
        #STO Corrected
        rename("$working_path/Output/${family}.out/${family}_${best}Corrected.stk", "$working_path/Output/${family}.out/${family}Corrected.stk");
    }
    if (-e $corrected_fasta && !-z $corrected_fasta){
        #STO Corrected
        rename("$working_path/Output/${family}.out/${family}_${best}FinalCorrected.fasta", "$working_path/Output/${family}.out/${family}FinalCorrected.fasta");
    }
    return;
}

sub get_best_folding {
    my ($sto0, $sto1, $current_dir) = @_;
    my $best;
    my $test1 = validate_secondary_structure_alignment($sto0, $current_dir);
    my $test2 = validate_secondary_structure_alignment($sto1, $current_dir);
    if ($test1 =~ m/^Valid/ && $test2 =~ m/^Valid/){ #TT
        $best = 0;
    } elsif ($test1 =~ m/^Valid/ && $test2 =~ m/^No_valid|^NA/){ #TF
        $best = 0;
    } elsif ($test1 =~ m/^No_valid|^NA/ && $test2 =~ m/^Valid/){ #FT
        $best = 1;
    } elsif ($test1 =~ m/^No_valid|^NA/ && $test2 =~ m/^No_valid|^NA/){ #FF
        $best = 0; #Here send 0, it would be in Classify.pm classified as No_valid Str
    } else {
        print_error("Not possible result at folding evaluation, check Validate.pm:get_best_folding");
    }
    return $best;
}

sub copy_best_element_subset {
    my ($data, $output_folder, $rfam_defined_models, $current_dir) = @_;
    foreach my $class (sort keys %{ $data }){
        if (exists $$data{$class}){    
            my $num = keys %{ $$data{$class} };
            if ($num == 2){ #Here evaluate the best folding
                my $code = $class; #(split /\s+|\t/, $$data{$class}{$model})[1];
                my $sto0 = "$output_folder/$rfam_defined_models/${rfam_defined_models}_0/$code/Output/${rfam_defined_models}_0.out/${rfam_defined_models}_0.stk";
                my $sto1 = "$output_folder/$rfam_defined_models/${rfam_defined_models}_1/$code/Output/${rfam_defined_models}_1.out/${rfam_defined_models}_1.stk";
                my $best = get_best_folding($sto0, $sto1, $current_dir);
                #if ($best eq "NA"){ #Both were wrong alignments > discard
                #    return;
                #}
                my $old = "$output_folder/$rfam_defined_models/${rfam_defined_models}_${best}/$code";
                my $new = "$output_folder/$rfam_defined_models/$code";
                copy_folders($old, $new);
                #Change name sto file
                rename_files($new, $rfam_defined_models, $best);
            } elsif ($num == 1){ #Guess who is the accepted one
		        my $code = $class; #(split /\s+|\t/, $$data{$class}{$model})[1];
                my $best = (sort(keys %{ $$data{$class} }))[0]; #Print the unique one #$$data{$class};#(split /\s+|\t/, $$data{$class}{$model})[-1]; 
                $best =~ s/(RF.*)([0-9]+)/$2/g; #Here detect which group was the best
                my $old = "$output_folder/$rfam_defined_models/${rfam_defined_models}_${best}/$code";
                my $new = "$output_folder/$rfam_defined_models/$code";
                copy_folders($old, $new);
                rename_files($new, $rfam_defined_models, $best);
            } elsif ($num == 0){ #Neither 0 or 1 validated, return empty
                print_error("The number of accepted candidates is wrong, please refer to Validate.pm:copy_best_element_subset");
            }
        } else {
            #The reference did exists: discard
            return;
        }
    }
    return;
}

sub print_results {
	my ($results, $working_path, $rfam_family) = @_;
	my $out_path = "$working_path/$rfam_family/$rfam_family.out";
	open my $OUT, ">> $out_path" or die;
	foreach my $classification (sort keys %{$results}){
	    my $specific = $$results{$classification}{$rfam_family}; #Only focused unique RFAM
        # Iterate over the codes included in this family
        foreach my $code (sort keys %{ $specific }){
			foreach my $subset (sort keys %{ $$specific{$code}}){
                print $OUT "$code\t$$specific{$code}{$subset}\t$classification\t$subset\n";
            }
        }
	}
	return;
}

sub clean_redundancy_file {
	my ($working_path, $family, $code, $address, $tag) = @_;
	my $file = "$working_path/$address/$tag/$family/$code/BaseFiles/${family}_genomes_list.txt";
	open my $IN, "< $file" or die "Your file $_ is not available to be cleaned\n";
	open my $OUT, "> $file.temp" or die "The output file couldn't be created\n";
	my %lines;
	while (<$IN>){
		chomp;
		next if ($_ =~ /^$/); #Avoid empty lines
		$lines{$_} = 1;
	}
	foreach my $line (sort keys %lines){
		print $OUT "$line\n";
	}
	close $OUT;
	close $IN;
	system "mv $file.temp $file";
	return;
}

sub clean_redundancy_file_subset {
	my ($working_path, $family, $code, $address, $subset, $tag) = @_;
	my $file = "$working_path/$address/$tag/$family/$subset/$code/BaseFiles/${subset}_genomes_list.txt";
	open my $IN, "< $file" or die "Your file $_ is not available to be cleaned\n";
	open my $OUT, "> $file.temp" or die "The output file couldn't be created\n";
	my %lines;
	while (<$IN>){
		chomp;
		next if ($_ =~ /^$/); #Avoid empty lines
		$lines{$_} = 1;
	}
	foreach my $line (sort keys %lines){
		print $OUT "$line\n";
	}
	close $OUT;
	close $IN;
	system "mv $file.temp $file";
	return;
}

no Moose::Role;
1;