package MiRNAture::BlastPrepareQueries;

use Moose::Role;
use Bio::SeqIO;
with 'MiRNAture::ToolBox'; #Import set of subroutines

=head1 detect_blast_queries
	Title: detect_blast_queries
	Usage: NAME-CONSTRUCTOR->detect_blast_queries
	Function: Normalize the fasta format from query(ies) from ncRNAs to 
		be searched on the subject genome.
	Returns: Creates new query files, with designed tag for each query
		specie. All the results are located on the output folder
		assigned to store the data to be blasted. 
=cut 

sub detect_blast_queries {
	my $queriesFolder = shift; 
	my @fasta_files = check_folder_files($queriesFolder, "\.fa\$\|\.fasta\$\|\.fna\$"); #Look for fasta files in reported folder
	my $numFiles = scalar @fasta_files;
	if ($numFiles == 0){
		print_error("Have not found fasta files on the declared query(ies) folder");
	}
	my $count = 1;
	my $speciesTag = generate_tags($queriesFolder); #Return hash based on queries_description.txt file
	foreach my $fastaF (@fasta_files){ #Iterate along all the defined query files
		next if ($fastaF =~ /\.new\.fasta/); #Skip database created files
		my ($tag, $sum); #Specie tag
		if (exists $$speciesTag{$fastaF}){
			$tag = $$speciesTag{$fastaF};
			$sum = generate_map_file_and_sizes("$queriesFolder/$fastaF", $tag, $count); #combine all to create new headers, based on mapping and generated tag for specie
			if ($sum =~ /^NA$/){ #Avoid to create map and sizes if exists
				print_process("$fastaF has been processed and indexed");
				next;
			} else {
				$count += 1; #$sum; #Avoid duplicated mapped headers
			}
		} 
		#else {
		#	print_warning("Sequence tag for query sequence: $fastaF is not available.");
		#}
	}
	if (-s "$queriesFolder/queries_description.txt"){
		open my $METADATA, "< $queriesFolder/queries_description.txt" or die "Metadata file is corrupted\n";
		while (<$METADATA>){
			chomp;
			check_fields_metadata_fasta($_);
		}
	} else {
		create_metadata_file_fasta($queriesFolder, \@fasta_files);
	}
	return $speciesTag;
}

=head1 generate_tags
	Title: generate_tags 
	Usage: generate_tags(queries_folder_path)
	Function: Check the existence of the file <queries_description.txt> 
		which would have the information about the origin of the 
		query files, based in the following format:
		<file.fasta name> <type_ncRNA_molecule> <Specie_name: (Genera specie)>
	Returns: It returns a hash reference for the files and the created tag for each
		specie.
=cut 

sub generate_tags {
	my $queriesFolder = shift;
	my %tags;
	if (-s "$queriesFolder/queries_description.txt"){
		open my $METADATA, "< $queriesFolder/queries_description.txt" or die;
		#snRNA.fasta	snRNA	Ciona intestinalis
		while (<$METADATA>){
			chomp;
			my @all = split /\t|\s+/, $_;
			my $file = $all[0]; #snRNA.fasta
			my $specie = "$all[-2] $all[-1]"; #Ciona intestinalis
			my $tag = create_tag_specie($specie);
			$tags{$file} = $tag;	
		}
	} else {
		$tags{"Unknown"} = "unkn";
	}
	return \%tags;
}

=head1 create_tag_specie
	Title: create_tag_specie
	Usage: create_tag_specie(specie name)
	Function: Create a tag based on the Genera and specie names from 
		the reported specie. This name must be the scientific name
		written as: Genera specie (i.e Homo sapiens).
	Returns: The created tag based on the first two letters of the 
		input name. (i.e Homo sapiens -> hosa)
=cut 

sub create_tag_specie {
	my $specie = shift;
	my $tag = lc($specie);
	$tag =~ s/^([a-z]{2})[a-z]+\s+([a-z]{2})[a-z]+$/$1$2/g; #Take first 2 letters both, genera and specie.
	return $tag;
}

=head1 create_metadata_file_fasta 
	Title: create_metadata_file_fasta
	Usage: create_metadata_file_fasta(specie name)
	Function: In cases where the queries_description.txt file has not
		been created, the program will generate a new one, based on the
		recovered fasta files from the indicated query(ies) folder.
	Returns: queries_description.txt file with the follow structure:
		<Name_fasta_file> <Unknown> <Unknown specie>
=cut 
sub create_metadata_file_fasta {
	my ($queriesFolder, $filesFasta) = @_;
	open my $METADATA, "> $queriesFolder/queries_description.txt" or die;
	foreach my $files (@$filesFasta){
		print $METADATA "$files\tmiRNA\tUnknown specie\n";	
	}
	close $METADATA;
	return;
}

=head1 generate_map_file_and_sizes
	Title: generate_map_file_and_sizes
	Usage: generate_map_file_and_sizes(fasta_file, tag_generated_hash, count)
	Function: Normalize the format of the provided query(ies) fasta files. 
		It creates a name database (*.map) file, and a new fasta with 
		new headers (*.new.fasta).
	Returns: Will generate a new fasta file with shortened headers and it
		corresponding mapping database. At the same time, returns the
		count variable to avoid repeat the headers.
=cut 

sub generate_map_file_and_sizes {
	my ($file, $short_lab, $count) = @_;
	my (%DB, %DB_Real,%sizesDB);
	if ((-e "$file.$short_lab.map" && !-z "$file.$short_lab.map") && (-e "$file.$short_lab.new.fasta" && -z "$file.$short_lab.new.fasta")){
		return "NA";
	}
	my $files2 = Bio::SeqIO->new(-file => $file, -format => "fasta");
	my $filename = $file;
	$filename =~ s/(.*)(\.fasta|\.fa|\.fna)/$1/g;
	while (my $files = $files2->next_seq){
		my $def = $count;
		$DB{$def} = $files->seq;
		my $id = $files->display_id;
		my $len = $files->length;
		$sizesDB{$def} = $len;
		$DB_Real{$def} = $files->id." ".$files->desc;
		$count++;	
	}
	open my $OUTF, "> $filename.$short_lab.new.fasta";
	foreach my $assign (keys %DB){
		print $OUTF ">".${short_lab}.$assign."\n".$DB{$assign}."\n";
	}
	open my $OUT, "> $file.$short_lab.map" or die;
	foreach my $assign (keys %DB_Real){
		print $OUT ${short_lab}.$assign."\t".$DB_Real{$assign}."\n";
	}
	open my $OUT2, "> $filename.$short_lab.new.fasta.len" or die;
	foreach my $label (keys %sizesDB){
		print $OUT2 ${short_lab}.$label."\t".$sizesDB{$label}."\n";
	}
	close $OUTF;
	close $OUT;
	close $OUT2;
	return $count;
}

=head1 check_fields_metadata_fasta {
	Title: check_fields_metadata_fasta
	Usage: check_fields_metadata_fasta(CURRENT_FILE_LINE)
	Function: Sanity check of the metadata file <queries_description.txt>.
		Evaluates line by line, each field, based on defined rules.
	Returns: Returns exception and die in case fields are not well formatted.
=cut 

sub check_fields_metadata_fasta {
	my $line = shift;
	my @all = split /\t|\s+/, $line;
	my $scientific_name = "$all[-2] $all[-1]";
	my ($check1, $check2, $check3);
	if ($all[0] =~ m/\.fa|\.fasta|\.fna/){
		$check1 = 1;
	} else {
		die "The line $_ seems to be bad formated: $all[0] should be a fasta file: <name.f(a|sta|na)>\n";
	}	
	if ($all[1] =~ m/antisense|antitoxin|Cis-reg|CRISPR|Intron|lncRNA|miRNA|misc_RNA|ribozyme|rRNA|snRNA|sRNA|tRNA|Unknown/i){ #Restricted only to ncRNA types on RFAM or Unknown.
		$check2 = 1;
	} else {
		die "The line $_ seems to be bad formated: $all[1] should be a recognized ncRNA type by RFAM: \n";
	}
	if ($scientific_name =~ m/^([A-Z]{1}[a-z]+)\s+([a-z]+)$|^Unknown specie$/){ #It has to be as scientific name species
		$check3 = 1;
	} else {
		die "The line $_ seems to be bad formated: $scientific_name should be a scientific name of the query specie as: Didemnum vexillum (Genera + specie)\n";
	}
	return;
}

no Moose::Role;
1;
