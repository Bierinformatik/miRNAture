#!/bin/bash
/home/bioinf/anaconda3/envs/mirnature/bin/blastn -db /home/bioinf/Desktop/miRNAture/Tutorial/Code/../Results/TemporalFiles/Latimeria_chalumnae/latimeria_chalumnae_genome.fa.new.fa -query /home/bioinf/Desktop/miRNAture/Tutorial/Code/../Data/QueriesToTest/Xenopus_laevis.new.fasta -num_threads 4 -outfmt 6 -out /home/bioinf/Desktop/miRNAture/Tutorial/Code/../Results/miRNA_prediction/Blast/Lach/Lach_6.miRNA.xela.tab