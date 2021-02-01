#!/bin/bash

final_results="../Results/miRNA_prediction/Final_Candidates/all_RFAM_Latch_Final.ncRNAs_homology.txt"
original_annotation="../Data/annotated_miRNAs_latch.txt"

if [[ -f $final_results ]]; then
    echo "Contig TotalReported Total TBlast THMM TInfernal"
    while read contig
    do
        total=$( grep -c "$contig" $final_results)
        totalANN=$( grep "$contig" $original_annotation| grep -c "\smiRNA\s")
        totalBlast=$( grep "$contig" $final_results | grep -c "Blast" )
        totalHMM=$( grep "$contig" $final_results | grep -c "HMM" )
        totalInfernal=$( grep "$contig" $final_results | grep -c "Infernal" )
        echo "$contig $totalANN $total $totalBlast $totalHMM $totalInfernal"
        # Here generate all the files to plot on GenomeBrowser
    done < ./list_contigs_test.txt
else
    echo "Did you generated the results from coelacanth? please use miRNAture"
    exit
fi

awk '{print $1".1\t"$6"\t"$7"\t"$11"\t"$8"\t"$2}' $final_results >> predicted_miRNAture.bed
