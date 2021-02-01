#!/bin/bash

current=$( pwd )
specie_tag="Latch"
specie_genome="$current/../Data/latimeria_chalumnae_genome.fa"
specie_name="Latimeria_chalumnnae"

echo "$specie_tag $specie_genome $specie_name"

workdir="$current/../Results"
mirfix_path="/homes/biertank/cristian/Projects/MIRfix/scripts/MIRfix.py"
mode="Blast,HMM,Infernal,Final"
#strategy="1,2,3,4,9,10,ALL"
strategy="9,10,ALL"
blastQueriesFolder="$current/../Data/QueriesToTest"

### Copy the correct list of CM, temporarly Basic_files/all_cm_scores_rfam14-2_metazoa_no_coelacanth.txt

cd $current/../../../Code/
#cp Data/Basic_files/all_cm_scores_rfam14-2_metazoa_no_coelacanth.txt Data/Basic_files/all_RFAM_scores.txt
###

# Run only homology-searches
#./miRNAture -stage homology -speG $specie_genome -speN $specie_name -speT $specie_tag -w $workdir -mfx $mirfix_path -m $mode -pe 0 -str $strategy -blastq $blastQueriesFolder
# Run detection matures
#./miRNAture -stage validation -speG $specie_genome -speN $specie_name -speT $specie_tag -w $workdir -mfx $mirfix_path -m $mode -pe 0 -str $strategy -blastq $blastQueriesFolder
# Run the complete analysis
#./miRNAture -stage evaluation -speG $specie_genome -speN $specie_name -speT $specie_tag -w $workdir -mfx $mirfix_path -m $mode -pe 0 -str $strategy -blastq $blastQueriesFolder
# Create summarise report
#./miRNAture -stage summarise -speG $specie_genome -speN $specie_name -speT $specie_tag -w $workdir -mfx $mirfix_path -m $mode -pe 0 -str $strategy -blastq $blastQueriesFolder

./miRNAture -stage complete -speG $specie_genome -speN $specie_name -speT $specie_tag -w $workdir -mfx $mirfix_path -m $mode -pe 0 -str $strategy -blastq $blastQueriesFolder