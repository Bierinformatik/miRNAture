#!/bin/bash

current=$( pwd )
specie_tag="Lach"
specie_genome="$current/../Data/latimeria_chalumnae_genome.fa"
specie_name="Latimeria_chalumnae"

workdir="$current/../Results"
mode="blast,rfam,mirbase,hmm,final"
strategy="5,6,ALL"
blast_queries_folder="$current/../Data/QueriesToTest"
#user_models="$current/User_Test_Data"
data_precalculated_folder="$current/Dataset_mirnature_tutorial"
sublistcms="${current}/list_miRNAs_to_search.txt"

# Run miRNAture complete
miRNAture -stage complete -sublist ${sublistcms} -nbitscore_cut 0.32 \
          -dataF ${data_precalculated_folder} -speG ${specie_genome} -speN ${specie_name} \
          -speT ${specie_tag} -w ${workdir} -m ${mode} -pe 1 -str ${strategy} -blastq ${blast_queries_folder} \
          -rep relax,150,100
