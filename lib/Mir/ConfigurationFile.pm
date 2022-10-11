package Mir::ConfigurationFile;

use Moose;
use MooseX::Types::Path::Class;
use YAML::Tiny;
use Bio::SeqIO;
use Data::Dumper;
#use File::Share ':all';

use lib "lib/MiRNAture";

with 'MiRNAture::ToolBox';

has 'stage' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'specie_name' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'specie_tag' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'specie_genome' => (
	is => 'ro',
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1,
);


has 'specie_genome_new' => (
	is => 'rw',
	isa => 'Path::Class::File',
	coerce => 1,
);

has 'specie_genome_new_length' => (
	is => 'rw',
	isa => 'Path::Class::File',
	coerce => 1,
);

has 'specie_genome_new_database' => (
	is => 'rw',
	isa => 'Path::Class::File',
	coerce => 1,
);

has 'data_path' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1,
);

has 'mode' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'blast_strategy' => (
	is => 'ro',
	isa => 'ArrayRef[Str]',
	lazy => 1, #Only called when attribute is accessed
	default => sub { #Test if previously defined mode
		my $shift =  shift;
		my $mode = $shift->mode or die "Mode didn't specified\n";
		if ($mode eq "blast"){
			return $shift->blast_strategy;
		} else {
			print_error("Declaring strategies only makes sense if you are running blast mode");
		}
	},
);

has 'blast_queries_path' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	lazy => 1, 
	default => sub { 
		my $shift =  shift;
		my $mode = $shift->mode or die "Mode didn't specified\n";
		if ($mode eq "blast"){
			return $shift->blast_queries_path->stringify;
		} else {
			print_error("Declaring the location of queries to blast only makes sense if you are running blast mode");
		}
	},
);

has 'output_folder' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1,
);

has 'current_folder' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	required => 1,
);

has 'mirfix_path' => (
    is => 'ro',
    isa => 'Maybe[Str|Path::Class::File]',
    required => 1,
    lazy => 0,
    default => sub {
        my $shift = shift;
        return "";
        },
    trigger => \&_path_set,
);

has 'config_file' => (
	is => 'rw',
	isa => 'Path::Class::File',
	coerce => 1
);

has 'parallel' => (
	is => 'ro',
	isa => 'Int',
	required => 1
);

has 'parallel_linux' => (
	is => 'ro',
	isa => 'Int',
	required => 1
);

has 'debug_mode' => (
	is => 'ro',
	isa => 'Int',
	required => 0
);

has 'user_folder' => (
    is => 'ro',
    isa => 'Maybe[Str|Path::Class::File]',
    required => 0,
    lazy => 1,
    default => sub {
        my $shift = shift;
        return;
        },
    trigger => \&_folder_set,
);

has 'model_list' => (
	is => 'ro',
	isa => 'Maybe[Str|Path::Class::File]',
	required => 1,
	lazy => 1,
	default => sub {
		my $shift = shift;
		if ($shift->mode->stringify eq "rfam"){
			#TODO: Change those references
			return $shift->data_path->stringify."/Data/RFAM_14-4/rfam_models_list.txt";
		} elsif ($shift->mode->stringify eq "mirbase"){
			return $shift->data_path->stringify."/Data/Mirbase/mirbase_models_list.txt";
		}
	},
	trigger => \&_size_set,
);

has 'repetition_rules' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
	lazy => 1,
    default => "default,200,100"
);

has 'nbitscore_cut' => (
	is => 'ro',
	isa => 'Num',
	required => 1
	#lazy => 1
	#default => 0.32 # By default this nBit= 0.32, cutoff bitscore
);

sub _folder_set {
    my ($self) = @_;
    if (!$self->user_folder){
	    return;
    }
    if (length $self->mirfix_path > 0){
        $self->{user_folder} = $self->user_folder;
    } else {
        ; 
    }
    return;
}

sub _path_set {
    my ($self) = @_;
    if (!$self->mirfix_path){
	return;
    }
    if (length $self->mirfix_path > 0){
        $self->{mirfix_path} = $self->mirfix_path;
    } else {
        ; 
    }
    return;
}

sub _size_set {
    my ($self) = @_;
    if (!$self->model_list){
		if ($self->mode =~ /^rfam$/){
			#TODO: Change those references
			$self->{model_list} = $self->data_path."/Data/RFAM_14-4/rfam_models_list.txt";
		} elsif ($self->mode =~ /^mirbase$/){
			$self->{model_list} = $self->data_path."/Data/Mirbase/mirbase_models_list.txt";
		} elsif ($self->mode =~ /,/){
			$self->{model_list} = $self->data_path."/Data/concatenated_models_list.txt";
		} elsif ($self->mode =~ /^user$/){
			$self->{model_list} = $self->user_folder."/user_models_list.txt";
		} elsif ($self->mode =~ /^final$/){
			$self->{model_list} = $self->data_path."/Data/RFAM_14-4/rfam_models_list.txt";
		}
        return;
    }
    if (length $self->model_list > 0){
        $self->{model_list} = $self->model_list;
    } 
    return;
}

# Create list of models provided for the user
sub infer_list_models_user {
	my $self = shift;
	my $path_scores = shift;
	if (!-e $path_scores || -z $path_scores){
		print_error("The user score file: $path_scores is incorrect");
	}
	my $outfile = $self->user_folder."/user_models_list.txt";
	if (!-e $outfile || -z $outfile){
		system "cut -f 1 $path_scores | sort | uniq >  $outfile";
	}
	return;
}

sub models_list_hash {
	my $self = shift;
	my %reference;
	my $mirbase = $self->data_path."/Data/Mirbase/mirbase_models_list.txt"; 
	my $rfam = $self->data_path."/Data/RFAM_14-4/rfam_models_list.txt";
	my $user = $self->user_folder."/user_models_list.txt";
	push @{$reference{'mirbase'}}, $mirbase;
	push @{$reference{'rfam'}}, $rfam;
	push @{$reference{'user'}}, $user; 
	push @{$reference{'hmm'}}, $rfam;	
	push @{$reference{'hmm'}}, $mirbase;
	push @{$reference{'blast'}}, $rfam;	
	push @{$reference{'blast'}}, $mirbase;
	push @{$reference{'final'}}, "NA";
	if (-e $user && !-z $user){
		push @{$reference{'hmm'}}, $user;
		push @{$reference{'blast'}}, $user;
	}
	return \%reference;
}


#Concatenate models in case both rfam and mirbase modes are indicated
sub concatenate_models {
	my $self =  shift;
	my $references = shift;
	my $mode = $self->mode;
	my @modes = split /,/, $mode;
	my $num = scalar @modes;
	# Dinamically create/update file based on indicated modes
	my @final_paths;
	my $out = $self->data_path."/Data/concatenated_models_list.txt";
	if (-e $out){
		system "rm $out";
		system "touch $out";
	}
	for (my $i = 0; $i <= $num - 1; $i++) { #Iterate modes
		if (exists $$references{$modes[$i]}){
			my $path = $$references{$modes[$i]}; # Array
			my $num2 = scalar @$path;
			for (my $j = 0; $j <= $num2 - 1; $j++) {
				push @final_paths, $$path[$j];
			}
		}
	}
	# Remove redundants
	my @unique_paths = do { my %seen; grep { !$seen{$_}++} @final_paths }; 
	for (my $l = 0; $l <= $#unique_paths; $l++) {
		unless ($unique_paths[$l] eq "NA"){
			system "cat $unique_paths[$l] >> $out";
		}
	}
	return;
}

sub generate_copy_genome {
	my $shift = shift;
	if (!-e $shift->output_folder){
        ##print_error("Output folder did not exist. Please create it to run miRNAture.");
        print_process("Output folder not detected. Creating...");
        create_folders($shift->output_folder, "");
	}
	my $data_folder = $shift->output_folder->stringify."/TemporalFiles"; #Genome folder in Temporal Folder
    create_folders($shift->output_folder->stringify, "TemporalFiles");
	create_folders($data_folder, $shift->specie_name);
	my $new_folder = $data_folder."/".$shift->specie_name;
	my $new_genome_file = $data_folder."/".$shift->specie_name."/".$shift->specie_genome->basename.".new.fa"; # Name of the new genome
	my $new_dababase_names = $data_folder."/".$shift->specie_name."/".$shift->specie_genome->basename.".map"; # Names DB
	my $new_genome_len = $new_genome_file.".len"; #Len Chr
	if (-e $new_genome_file && !-z $new_genome_file){ #Test if exists or process
		print_process("miRNAture detected a previous genome reference for the ".$shift->specie_name." genome");
		#generate_length_file_genome($new_genome_file, $shift->specie_tag, $data_folder."/".$shift->specie_name,$shift->specie_genome->basename.".new.fa");
	} else {
		print_process("Generating a carbon-copy of your ".$shift->specie_tag." genome");
		create_map_genome($shift->specie_genome, $shift->specie_tag, $data_folder."/".$shift->specie_name, $shift->specie_genome->basename);
		generate_length_file_genome($new_genome_file, $shift->specie_tag, $data_folder."/".$shift->specie_name, $shift->specie_genome->basename.".new.fa");
		#copy_files($shift->specie_genome, $new_folder);
	}
	$shift->specie_genome_new($new_genome_file);	
	$shift->specie_genome_new_database($new_dababase_names);
	$shift->specie_genome_new_length($new_genome_len);
	return;
}

sub generate_length_file_genome {
	my ($in, $short, $data_path, $name) = @_;
	open my $NEWLEN, "> $data_path/$name.len" or die;
	my $file = Bio::SeqIO->new(-file => $in, -format => "fasta");
	while (my $seq = $file->next_seq){ ## selects one sequence at a time
		my $id = $seq->display_id;
		my $len = $seq->length;
		print $NEWLEN "$id\t$len\n";
	}
	close $NEWLEN;
	return; 
}

sub create_map_genome {
	my ($in, $short, $data_path, $name) = @_;
	my %NR;
	my %NRori;
	my $file = Bio::SeqIO->new(-file => $in, -format => "fasta");
	my $num = 1;
	while (my $fasta = $file->next_seq) {
		my $def = $num;
		$NR{$def} = $fasta->seq;
		$NRori{$num} = $fasta->id ." " .$fasta->desc;
		$num++;
	}
	# Generate new fasta with labels
	open my $NEWFA, "> $data_path/$name.new.fa" or die;
	foreach my $defs(keys %NR) {
		my $seq = $NR{$defs};
		$seq =~ s/(.{60})/$1\n/g; #include jump each 60 char
		print $NEWFA ">".$short.$defs."\n".$seq."\n";
	}
	close $NEWFA;
	# Generate database file 
	open my $OUTMAP, "> $data_path/$name.map" or die;
	foreach my $defs(keys %NRori) {
		print $OUTMAP $short.$defs."\t".$NRori{$defs}."\n";
	}
	close $OUTMAP;
	return;
}

sub test_name {
	my $name = shift;
	$name =~ s/(.*)(\.yaml)/$1/g;
	my $count = 1;
	NAME:
	my $test_name = "$name-$count.yaml";
	if (-e $test_name){
		$count++;
		goto NAME;        
	} else {
		return $test_name;
	}    
}

#sub read_data_share {
#	my $data_location = dist_dir('Bio-miRNAture');
#	#print "$data_location\n";
#	return $data_location;
#}

sub create_config_file {
	my $shift = shift;
	my $yamlNew = YAML::Tiny->new();
	my $final_file = $shift->output_folder->stringify."/miRNAture_configuration_".$shift->specie_tag.".yaml"; #Allow multiple experiments at the same time
	if (-e $final_file && !-z $final_file){
		print_process("A previous configuration file was detected, then miRNAture will create an alternative one.");
		my $final_file_alternative = test_name($final_file);
		$yamlNew->write($final_file_alternative);
		$shift->config_file($final_file_alternative);
	} else {
		$yamlNew->write($final_file);
		$shift->config_file($final_file);
	}
	return;
}

sub write_header {
	my $shift = shift;
	my $yaml_file = YAML::Tiny->read($shift->config_file);
	my ($makeblastdb, $blastn, $nhmmer, $cmsearch, $cmcalibrate, $cmbuild, $clustalo, $RNAalifold, $MIRfix, $mirnatureHomology, $mirnanchor) = obtain_paths_programs();
	$yaml_file->[0]->{"miRNAture"}{"Author"} = 'Cristian A. Velandia-Huerto, Joerg Fallmann, Peter F. Stadler';
	$yaml_file->[0]->{"miRNAture"}{"Version"} = 'v.1.1';
	$yaml_file->[0]->{"miRNAture"}{"Date"} = 'Sept 21 2022';
	#$yaml_file->[0]->{"miRNAture"}{"Date"} = 'Mon Mar 16 19:02:25 CET 2020';
	$yaml_file->[1]->{"Data_user"}{"User"} = `whoami | tr -d '\n'`;
	$yaml_file->[1]->{"Data_user"}{"Hostname"} = `hostname | tr -d '\n'`;
	$yaml_file->[1]->{"Data_user"}{"Running_date"} = `date | tr -d '\n'`;
	$yaml_file->[2]->{"Program_locations"}{"miRNAture"} = $mirnatureHomology;
	$yaml_file->[2]->{"Program_locations"}{"miRNAnchor"} = $mirnanchor;
	#$yaml_file->[2]->{"Program_locations"}{"miRNAlien"} = $shift->current_folder->stringify."/src/miRNAlien.pl";
	$yaml_file->[2]->{"Program_locations"}{"makeblastdb"} = $makeblastdb;
	$yaml_file->[2]->{"Program_locations"}{"blastn"} = $blastn;
	$yaml_file->[2]->{"Program_locations"}{"nhmmer"} = $nhmmer;
	$yaml_file->[2]->{"Program_locations"}{"cmsearch"} = $cmsearch;
	$yaml_file->[2]->{"Program_locations"}{"cmcalibrate"} = $cmcalibrate;
	$yaml_file->[2]->{"Program_locations"}{"cmbuild"} = $cmbuild;
	$yaml_file->[2]->{"Program_locations"}{"clustalo"} = $clustalo;
	$yaml_file->[2]->{"Program_locations"}{"RNAalifold"} = $RNAalifold;
	if (length $shift->mirfix_path > 0){
		$yaml_file->[2]->{"Program_locations"}{"MIRfix"} = $shift->mirfix_path;
	} else {
		$yaml_file->[2]->{"Program_locations"}{"MIRfix"} = $MIRfix;
	}
	$yaml_file->write($shift->config_file);
	return;
}

sub obtain_paths_programs {
	my $makeblastdb = collect_data("makeblastdb");
	my $blastn = collect_data("blastn");
	my $nhmmer = collect_data("nhmmer");
	my $cmsearch = collect_data("cmsearch");
	my $cmcalibrate = collect_data("cmcalibrate");
	my $cmbuild = collect_data("cmbuild");
	my $clustalo = collect_data("clustalo");
	my $RNAalifold = collect_data("RNAalifold");
	my $mirfix = collect_data("MIRfix.py");
    ###my $mirnatureHomology = collect_data("miRNAture.pl");
    ###my $mirnanchor = collect_data("miRNAnchor.pl");
    ##TEMPORAL
    my $mirnatureHomology = "/homes/biertank/cristian/Projects/miRNAture_v1/script/miRNAture.pl";
    my $mirnanchor = "/homes/biertank/cristian/Projects/miRNAture_v1/script/miRNAnchor.pl";
	##TODO
    return ($makeblastdb, $blastn, $nhmmer, $cmsearch, $cmcalibrate, $cmbuild, $clustalo, $RNAalifold, $mirfix, $mirnatureHomology, $mirnanchor);
}

sub collect_data {
	my $program = shift;
	my $path = `which $program 2> /dev/null | tr -d '\n'`;
	if ($path){
		return $path;
	} else {
		print_error("The required program: $program has not been found, please install it");
	}
}

sub write_config_file {
	my $shift = shift;
	my $data_path = shift;
	my $yaml = YAML::Tiny->read($shift->config_file);
	$yaml->[3]->{Specie_data}{"Tag"} = $shift->specie_tag;
	$yaml->[3]->{Specie_data}{"Name"} = $shift->specie_name;
	$yaml->[3]->{Specie_data}{"Old_Genome"} = $shift->specie_genome->stringify;
	$yaml->[3]->{Specie_data}{"Genome"} = $shift->specie_genome_new->stringify;
	$yaml->[3]->{Specie_data}{"Genome_length"} = $shift->specie_genome_new_length->stringify;
	$yaml->[3]->{Specie_data}{"Database_names_genome"} = $shift->specie_genome_new_database->stringify;
	$yaml->[3]->{Default_folders}{"Current_dir"} = $shift->current_folder->stringify;
	$yaml->[3]->{Default_folders}{"Output_folder"} = $shift->output_folder->stringify;
	#$yaml->[3]->{Default_folders}{"Temp_folder"} = $shift->output_folder->stringify."/Temp";
	
	$yaml->[3]->{Default_folders}{Pre_calculated_validation_data} = "$data_path/Data/Validation_mature_data";
	$yaml->[3]->{Default_folders}{"Data_folder"} = "$data_path/Data";
	$yaml->[3]->{Default_folders}{"CM_folder"} = "$data_path/Data/RFAM_14-4/CMs"; #Default Models
	$yaml->[3]->{Default_folders}{"HMM_folder"} = "$data_path/Data/RFAM_14-4/HMMs"; #Modified Lach
	$yaml->[3]->{Default_folders}{"Other_CM_folder"} = "$data_path/Data/Mirbase/CMs";
	$yaml->[3]->{Default_folders}{"Other_HMM_folder"} = "$data_path/Data/Mirbase/HMMs";

	$yaml->[3]->{Default_folders}{"User_folder"} = $shift->user_folder;
	if (length $shift->user_folder > 0) {
		$yaml->[3]->{Default_folders}{"User_CM_folder"} = $shift->user_folder."/CMs";
		$yaml->[3]->{Default_folders}{"User_HMM_folder"} = $shift->user_folder."/HMMs";
	} else { # Case when the user folder is not defined
		$yaml->[3]->{Default_folders}{"User_CM_folder"} = "NA";
		$yaml->[3]->{Default_folders}{"User_HMM_folder"} = "NA";
	}

	$yaml->[3]->{Default_folders}{"Basic_files_miRNAture"} = "$data_path/Data/Basic_files";

	$yaml->[3]->{Default_folders}{"List_cm_miRNAs"} = $shift->model_list; 
	$yaml->[3]->{Default_folders}{"Blast_queries"} = $shift->blast_queries_path->stringify;
	$yaml->[3]->{Homology_options}{"Mode"} = $shift->mode;
	$yaml->[3]->{Homology_options}{"Parallel"} = $shift->parallel;
	$yaml->[3]->{Homology_options}{"Parallel_linux"} = $shift->parallel_linux;
	if ($shift->blast_strategy){
		$yaml->[3]->{Homology_options}{"Blast_strategies"} = $shift->blast_strategy;
	}
	if (length $shift->repetition_rules == 0){
		$yaml->[3]->{Homology_options}{Repetition_threshold} = "default,200,100";
	} else {
		$yaml->[3]->{Homology_options}{Repetition_threshold} = $shift->repetition_rules;
	}
	if (length $shift->nbitscore_cut == 0){
		$yaml->[3]->{Homology_options}{Threshold_bitscore} = "0.32";
	} else {
		$yaml->[3]->{Homology_options}{Threshold_bitscore} = $shift->nbitscore_cut;
	}
	$yaml->write($shift->config_file);
	return;
}

sub read_final_file {
	my $shift = shift;
	my $final = YAML::Tiny->read($shift->config_file);
	return $final;
}

sub read_last_file {
	my $shift = shift;
	my $final_file = $shift->output_folder->stringify."/miRNAture_configuration_".$shift->specie_tag.".yaml";
    my $last_name = get_last_name($final_file);
	my $final = YAML::Tiny->read($last_name);
	return $final;
}

sub get_last_name {
	my $name = shift;
    my $nameOriginal = $name;
	$name =~ s/(.*)(\.yaml)/$1/g;
	my $count = 1;
	NAME:
	my $test_name = "$name-$count.yaml";
    my $last_name;
	if (-e $test_name && !-z $test_name){
		$last_name = $test_name;
		$count++;
		goto NAME;        
	} else {
		if (length $last_name > 0){
			return $last_name;
		} else {
			return $nameOriginal;
		}
	}    
}
sub start {
	print "\n";
	print ".__________________________________________________________________________.\n";
	print "|                                                                          |\n";
	print "|________________________________________________.........---------........|\n";
	print "             .--. .   .    .    .                ._________________________.\n";
       	print "          o  |   )|\\  |   / \\  _|_               |                         |\n";
       	print ".--.--.   .  |--' | \\ |  /___\\  |  .  . .--..-.  |........---------........|\n";
       	print "|  |  |   |  |  \\ |  \\| /     \\ |  |  | |  (.-'  ._______.         ._______.\n";
       	print "'  '  `--' `-'   `'   ''       ``-'`--`-'   `--' |.......|         |.......|\n";
	print "Computational detection of microRNA candidates\n";
    	#print "v.1.0 Mar 16, 2020\n";
        #      "v.1.0 Feb 1, 2021\n";
	print "v.1.1 Sept 21th, 2022\n";
	print "Cristian A. Velandia-Huerto, Jöerg Fallmann, Peter F. Stadler\n";
	print "Bioinformatics Leipzig\n";
	print "University of Leipzig\n";
	print "\n";
	return;
}

no Moose;
__PACKAGE__->meta->make_immutable;
