package MiRNAnchor::Main;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(identify_RFAM_families check_existence_folder_output delete_temporal_folders get_genome_validation_list infer_list_genomes recognize_families_homology);

use Moose;
use MooseX::Types::Path::Class;
use Data::Dumper;
use File::Copy;

has 'fasta_sequences' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1
);

has 'mirfix_path' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1
);

has 'output_folder' => (
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

has 'parallel_running_confirmation' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'current_dir' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1
);

has 'subject_genome_path' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1
);

has 'tag_spe_query' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

with 'MiRNAnchor::Tools';

sub recognize_families_homology {
	my $shift = shift;
	my $db_models_relation = shift;
	my $specie = $shift->tag_spe_query;
	my $targetfolder = $shift->output_folder->stringify;
	my $table_db_target = "$targetfolder/miRNA_prediction/Final_Candidates/all_RFAM_${specie}_final.ncRNAs_homology.txt.db";
	my $target_fasta = $shift->fasta_sequences->stringify;
	modify_table($table_db_target, $db_models_relation);
	modify_fasta($target_fasta, $specie, $db_models_relation);
	return;
}

sub identify_RFAM_families {
	my $shift = shift;
	my %rfam_families_list;	
	my @files = check_folder_files_miranchor($shift->fasta_sequences->stringify, "\_".$shift->tag_spe_query."\.\.\*\.fa|\_".$shift->tag_spe_query."\.\.\*"."\.fasta"); #Get all fasta files
	if (scalar @files == 0){
		die "Does not exists fasta candidates in the declared folder\n";
	}
	foreach my $file_fasta (@files)	{
		if ($file_fasta =~ m/^all/){ #all_RFAM_Hosa.MIPF0000725.fasta
			my $code_rfam = (split /\./, $file_fasta)[1];
			$rfam_families_list{$code_rfam} = 1;	
		} else {
			next;
		}
	}
	return \%rfam_families_list;
}	

sub check_existence_folder_output {
	my $shift = shift;
	if (!-d $shift->output_folder->stringify){
		print_process("Creating the Output folder for you");
		create_folder($shift->output_folder->stringify);
	}
	return;
}

sub delete_temporal_folders {
	my $shift = shift;
	my $tempdir = shift;
	my $folder = $shift->current_dir->stringify."/".$tempdir;
	system "rm -r $tempdir";
	return;	
}

sub get_genome_validation_list {
	my $shift = shift;
	my $genomes_folder = $shift->precalculated_path->stringify."/Genomes_miRBase";
	my $genomes_list;
	# Check the existence of genomes folder
	if (-d $genomes_folder){
		$genomes_list = infer_list_genomes($genomes_folder);	
	} else {
		print_error("The genomes folder Genomes_miRBase did not exists! please copy accordingly inside the Pre-calculated folder ".$shift->precalculated_path->stringify);
	}
    return $genomes_list;
}

sub infer_list_genomes {
	my $path_genomes = shift;
	my @final_fasta_files;
	opendir(DIR, $path_genomes) or die "can't opendir $path_genomes: $!";
	my @folders = grep { (/\w+/) } readdir(DIR);
	foreach my $folder (@folders){
		if (-d $path_genomes."/".$folder){
			opendir(DIR2, $path_genomes."/".$folder) or die;
			my @files_specific = grep { (/\.fa$|\.fasta$|\.gz$/ ) && -f "$path_genomes/$folder/$_" } readdir(DIR2);
			foreach my $file (@files_specific){
				$file = "$path_genomes/$folder/$file";
				push @final_fasta_files, "$file";	
			}
			close DIR2;
		}
	}
	close DIR;
	return \@final_fasta_files;
}

sub modify_table {
	my $table = shift;
	my $db =  shift;
	open my $TAB, "< $table" or die;
	open my $NEWF, "> $table.temp" or die; 
	while (<$TAB>){
		chomp;
		my @all = split /\s+|\t/, $_;
		$all[3] = reeplace_model($all[3], $db);
		my $new = join "\t", @all;
		print $NEWF "$new\n";
	}
	close $NEWF;
	move("$table.temp", $table);
	return;
}

sub check_folder_files_index {
	my ($dir, $prefix) = @_;
	opendir(DIR, $dir);
	my @files = grep(/$prefix$/, readdir(DIR));
	closedir(DIR);
	return \@files;
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

sub move_fasta {
	my ($path, $fasta) = @_;
	create_folders($path, "Temp");
	foreach my $file (@$fasta){
		move("$path/$file", "$path/Temp/");
	}
	return;
}

sub modify_fasta {
	my ($path, $specie, $db) = @_;
	#Move original fasta to the Temp/ folder
	my $fasta_query = check_folder_files_index($path, ".*\_".$specie."\.*\.fasta");
	move_fasta($path, $fasta_query);
	#Index again on the new folder
	my $fasta_moved_files = check_folder_files_index("$path/Temp", ".*\_".$specie."\.*\.fasta");
	foreach my $file (@$fasta_moved_files){
		my $model = (split /\./, $file)[1]; #all_RFAM_Hosa.F-let-7-1.fasta
		my $new_name = $model;
		if (exists $$db{$model}){ #Assuming 1:1 on Rfam models
			$new_name = $$db{$model};
		}
		#Create new file with new model concatenated in the fasta folder
		open my $OUT, ">> $path/all_RFAM_$specie.$new_name.fasta" or die;
		open my $IN, "< $path/Temp/$file" or die;
		while (<$IN>){
			chomp;
			print $OUT "$_\n";
		}
		close $IN;
		close $OUT;
	}
	return;
}

sub reeplace_model {
	my $modelold = shift;
	my $db = shift;
	my $modelnew = $modelold;
	# Takes the last value provided in the files
	if (exists $$db{$modelold}){
		$modelnew = $$db{$modelold};
	}
	return $modelnew;
}

no Moose;
__PACKAGE__->meta->make_immutable;
