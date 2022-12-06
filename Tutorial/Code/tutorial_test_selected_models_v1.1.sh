#!/bin/bash

current=$( pwd )
species_tag="Lach"
species_genome="$current/../Data/latimeria_chalumnae_genome.fa"
species_name="Latimeria_chalumnae"

workdir="$current/../Results"
mode="blast,rfam,mirbase,hmm,final"
strategy="5,6,ALL"
blast_queries_folder="$current/../Data/QueriesToTest"
#user_models="$current/User_Test_Data"
data_precalculated_folder="$current/Dataset_mirnature_tutorial"
sublistcms="${current}/list_miRNAs_to_search.txt"

# Run miRNAture complete
miRNAture -stage complete -sublist ${sublistcms} -nbitscore_cut 0.32 \
          -dataF ${data_precalculated_folder} -speG ${species_genome} -speN ${species_name} \
          -speT ${species_tag} -w ${workdir} -m ${mode} -pe 1 -str ${strategy} -blastq ${blast_queries_folder} \
          -rep relax,150,100
