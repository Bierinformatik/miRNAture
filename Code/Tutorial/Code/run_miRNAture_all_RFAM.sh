#!/bin/bash

current=$( pwd )
specie_tag="Lach"
specie_genome="$current/../Data/latimeria_chalumnae_genome.fa"
specie_name="Latimeria_chalumnae"

echo "$specie_tag $specie_genome $specie_name"

workdir="$current/../Results"
mode="Blast,HMM,Infernal,Other_CM,Final"
strategy="5,6,ALL"
blastQueriesFolder="$current/../Data/QueriesToTest"

cd $current/../../../Code/

# Run only homology-searches
#./miRNAture -stage homology -sublist $current/list_miRNAs_to_search.txt -speG $specie_genome -speN $specie_name -speT $specie_tag -w $workdir -m $mode -pe 0 -str $strategy -blastq $blastQueriesFolder -rep default,150,100
# Run detection matures
#./miRNAture -stage validation -speG $specie_genome -speN $specie_name -speT $specie_tag -w $workdir -m $mode -pe 0
# Run the complete analysis
#./miRNAture -stage evaluation -speG $specie_genome -speN $specie_name -speT $specie_tag -w $workdir -m $mode -pe 0
# Create summarise report
#./miRNAture -stage summarise -speG $specie_genome -speN $specie_name -speT $specie_tag -w $workdir -m $mode -pe 0

# Run miRNAture complete
./miRNAture -stage complete -sublist $current/list_miRNAs_to_search.txt -speG $specie_genome -speN $specie_name -speT $specie_tag -w $workdir -m $mode -pe 0 -str $strategy -blastq $blastQueriesFolder -rep default,150,100
