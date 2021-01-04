package Mir::ConfigurationFile;

use Moose;
use MooseX::Types::Path::Class;
use YAML::Tiny;
use Data::Dumper;
use Bio::SeqIO;
use lib "lib/miRNAture";

with 'miRNAture::ToolBox';

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
		if ($mode eq "BLAST"){
			return $shift->blast_strategy;
		} else {
			print_error("Declaring strategies only makes sense if you are running BLAST mode");
		}
	},
);

has 'blast_queries_path' => (
	is => 'ro',
	isa => 'Path::Class::Dir',
	coerce => 1,
	lazy => 1, #Only called when attribute is accessed
	default => sub { #Test if previously defined mode
		my $shift =  shift;
		my $mode = $shift->mode or die "Mode didn't specified\n";
		if ($mode eq "BLAST"){
			return $shift->blast_queries_path->stringify;
		} else {
			print_error("Declaring the location of queries to blast only makes sense if you are running BLAST mode");
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
	isa => 'Path::Class::File',
	coerce => 1,
	required => 1,
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

has 'model_list' => (
    is => 'ro',
    isa => 'Maybe[Str|Path::Class::File]',
    required => 1,
    lazy => 1,
    default => sub {
        my $shift = shift;
        return $shift->current_folder->stringify."/Data/miRNA_RFAM14-2_ACC_metazoa_no_coelacanth.lista";
        },
    trigger => \&_size_set,
);

has 'repetition_rules' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
    default => 'default,200,100',
);

sub _size_set {
    my ( $self, $size, $old_size ) = @_;
    if (!$self->model_list){
        $self->{model_list} = $self->current_folder->stringify."/Data/miRNA_RFAM14-2_ACC_metazoa_no_coelacanth.lista";
        return;
    }
    if (length $self->model_list > 0){
        $self->{model_list} = $self->model_list;
    } else {
        $self->{model_list} = $self->current_folder->stringify."/Data/miRNA_RFAM14-2_ACC_metazoa_no_coelacanth.lista";
    }
    return;
}

sub generate_copy_genome {
	my $shift = shift;
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

sub create_config_file {
	my $shift = shift;
	my $yamlNew = YAML::Tiny->new();
	my $final_file = $shift->current_folder->stringify."/miRNAture_configuration_".$shift->specie_tag.".yaml"; #Allow multiple experiments at the same time
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
	my ($makeblastdb, $blastn, $nhmmer, $cmsearch, $cmcalibrate, $cmbuild, $clustalo, $RNAalifold) = obtain_paths_programs();
	$yaml_file->[0]->{"miRNAture"}{"Author"} = 'Cristian A. Velandia-Huerto';
	$yaml_file->[0]->{"miRNAture"}{"Version"} = 'v.1.0';
	$yaml_file->[0]->{"miRNAture"}{"Date"} = 'Mon Mar 16 19:02:25 CET 2020';
	$yaml_file->[1]->{"Data_user"}{"User"} = `whoami | tr -d '\n'`;
	$yaml_file->[1]->{"Data_user"}{"Hostname"} = `hostname | tr -d '\n'`;
	$yaml_file->[1]->{"Data_user"}{"Running_date"} = `date | tr -d '\n'`;
	$yaml_file->[2]->{"Program_locations"}{"miRNAture"} = $shift->current_folder->stringify."/src/miRNAture.pl";
	$yaml_file->[2]->{"Program_locations"}{"miRNAnchor"} = $shift->current_folder->stringify."/src/miRNAnchor.pl";
	#$yaml_file->[2]->{"Program_locations"}{"miRNAlien"} = $shift->current_folder->stringify."/src/miRNAlien.pl";
	$yaml_file->[2]->{"Program_locations"}{"makeblastdb"} = $makeblastdb;
	$yaml_file->[2]->{"Program_locations"}{"blastn"} = $blastn;
	$yaml_file->[2]->{"Program_locations"}{"nhmmer"} = $nhmmer;
	$yaml_file->[2]->{"Program_locations"}{"cmsearch"} = $cmsearch;
	$yaml_file->[2]->{"Program_locations"}{"cmcalibrate"} = $cmcalibrate;
	$yaml_file->[2]->{"Program_locations"}{"cmbuild"} = $cmbuild;
	$yaml_file->[2]->{"Program_locations"}{"clustalo"} = $clustalo;
	$yaml_file->[2]->{"Program_locations"}{"RNAalifold"} = $RNAalifold;
	$yaml_file->[2]->{"Program_locations"}{"MIRfix"} = $shift->mirfix_path->stringify;
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
	return ($makeblastdb, $blastn, $nhmmer, $cmsearch, $cmcalibrate, $cmbuild, $clustalo, $RNAalifold);
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
    $yaml->[3]->{Default_folders}{Pre_calculated_validation_data} = $shift->current_folder->stringify."/Data/ValidationMature"; # Let-7 Experiment
    #$yaml->[3]->{Default_folders}{Pre_calculated_validation_data} = $shift->current_folder->stringify."/Data/ValidationMature/HumanValidation"; #Hosa Experiment
    #$yaml->[3]->{Default_folders}{Pre_calculated_validation_data} = $shift->current_folder->stringify."/Data/ValidationMature/ReAnnotation"; #Re_annotation Experiment
	$yaml->[3]->{Default_folders}{"Data_folder"} = $shift->current_folder->stringify."/Data";
    #$yaml->[3]->{Default_folders}{"RFAM_folder"} = $shift->current_folder->stringify."/Data/RFAM_14-2";
    #$yaml->[3]->{Default_folders}{"CM_folder"} = $shift->current_folder->stringify."/Data/RFAM_14-2/CM/miRNANoLatimeria"; #Modified Lach
    $yaml->[3]->{Default_folders}{"CM_folder"} = $shift->current_folder->stringify."/Data/RFAM_14-2/CM/All"; #Default Models
    #$yaml->[3]->{Default_folders}{"CM_folder"} = "/scr/k70san/bogota_unal/miRNAturePaper/Code/ReAnnotationFilter/Data/SelectedWeirdFamilies/CM"; #Re_annotation exp. No seq
    $yaml->[3]->{Default_folders}{"CM_folder"} = "/scr/k70san/bogota_unal/miRNAturePaper/Code/Let7validationHuman/Data/CM/Calibrated"; #Let7
    #$yaml->[3]->{Default_folders}{"Other_CM_folder"} = $shift->current_folder->stringify."/Data/RFAM_14-2/CM/Metazoa/miRNAs";
    $yaml->[3]->{Default_folders}{"Other_CM_folder"} = "/scr/k70san/bogota_unal/miRNAturePaper/Code/Let7validationHuman/Data/CM/Calibrated"; #Let7
    #$yaml->[3]->{Default_folders}{"Other_CM_folder"} = "/scr/k70san/bogota_unal/miRNAturePaper/Code/ReAnnotationFilter/Data/SelectedWeirdFamilies/CM"; #OtherFAmilies
    #$yaml->[3]->{Default_folders}{"Other_CM_folder"} = "/scr/k70san/bogota_unal/miRNAturePaper/Code/ReAnnotationFilter/Data/SelectedWeirdFamilies/CM"; #ReAnnotation exp. No seq
    #$yaml->[3]->{Default_folders}{"Other_CM_folder"} = "/scr/k70san/bogota_unal/miRNAturePaper/Code/AnnotationHomoSapiens/Data/miRBase22/NoHuman/CM"; #OtherFAmilies without human
    #$yaml->[3]->{Default_folders}{"HMM_folder"} = $shift->current_folder->stringify."/Data/RFAM_14-2/HMM/miRNANoLatimeria"; #Modified Lach
    #$yaml->[3]->{Default_folders}{"HMM_folder"} = $shift->current_folder->stringify."/Data/RFAM_14-2/HMM/All"; #Modified Lach
    #$yaml->[3]->{Default_folders}{"HMM_folder"} = "/scr/k70san/bogota_unal/miRNAturePaper/Code/ReAnnotationFilter/Data/SelectedWeirdFamilies/HMM"; #Reannotation Exp. NoSeq
    $yaml->[3]->{Default_folders}{"HMM_folder"} = "/scr/k70san/bogota_unal/miRNAturePaper/Code/AnnotationHomoSapiens/Data/miRBase22/NoHuman/HMM"; #withouthuman
    $yaml->[3]->{Default_folders}{"Basic_files_miRNAture"} = $shift->current_folder->stringify."/Data/Basic_files";
	$yaml->[3]->{Default_folders}{"List_cm_miRNAs"} = $shift->model_list; 
	$yaml->[3]->{Default_folders}{"Blast_queries"} = $shift->blast_queries_path->stringify;
	$yaml->[3]->{Homology_options}{"Mode"} = $shift->mode;
	$yaml->[3]->{Homology_options}{"Parallel"} = $shift->parallel;
	if ($shift->blast_strategy){
		$yaml->[3]->{Homology_options}{"Blast_strategies"} = $shift->blast_strategy;
	}
    $yaml->[3]->{Homology_options}{Repetition_threshold} = $shift->repetition_rules;
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
    my $final_file = $shift->current_folder->stringify."/miRNAture_configuration_".$shift->specie_tag.".yaml"; #Allow multiple experiments at the same time
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
	if (-e $test_name){
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
	my $file = "lib/Mir/.banner2";
	open my $IN, "< $file" or die;
	print "\n";
	while (<$IN>){
		chomp;
		print $_."\n";
	}
	close $IN;
	print "Computational detection of microRNA candidates\n";
	print "v.1.0 Mar 16, 2020\n";
	print "Cristian A. Velandia-Huerto\n";
	print "Bioinformatics Leipzig\n";
	print "University of Leipzig\n";
	print "\n";
	return;
}

no Moose;
__PACKAGE__->meta->make_immutable;
