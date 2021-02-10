#!/bin/bash

####
# miRNAture
# Cristian A. Velandia-Huerto, Joerg Fallmann, Peter F. Stadler
# Feb 1 2021
# v.1.0
####
# Session data:
# Hostname: unknownd89c6795a673
# Date:Wed Feb 10 13:57:37 CET 2021
# User:bioinf
# Program: miRNAture
####
#BLAST searches
/home/bioinf/perl5/bin/miRNAture.pl  -l /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/RFAM_14-4/rfam_models_list.txt -m BLAST -str 5,6,ALL -blstq /home/bioinf/Desktop/miRNAture/Tutorial/Code/../Data/QueriesToTest -pe 0 -spe Lach -n_spe Latimeria_chalumnae -w /home/bioinf/Desktop/miRNAture/Tutorial/Code/../Results/miRNA_prediction -data /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data -cmp /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/RFAM_14-4/CMs /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/Other_CM /CMs -hmmp /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/RFAM_14-4/HMMs /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/Other_HMM /HMMs -rep default,150,100
#HMM searches
/home/bioinf/perl5/bin/miRNAture.pl -l /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/RFAM_14-4/rfam_models_list.txt -m HMM -pe 0 -spe Lach -n_spe Latimeria_chalumnae -w /home/bioinf/Desktop/miRNAture/Tutorial/Code/../Results/miRNA_prediction -data /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data -cmp /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/RFAM_14-4/CMs /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/Other_CM /CMs -hmmp /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/RFAM_14-4/HMMs /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/Other_HMM /HMMs -rep default,150,100
#Infernal searches
/home/bioinf/perl5/bin/miRNAture.pl -l /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/RFAM_14-4/rfam_models_list.txt -m INFERNAL -pe 0 -spe Lach -n_spe Latimeria_chalumnae -w /home/bioinf/Desktop/miRNAture/Tutorial/Code/../Results/miRNA_prediction -data /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data -cmp /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/RFAM_14-4/CMs -hmmp /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/RFAM_14-4/HMMs -rep default,150,100
#Other searches
/home/bioinf/perl5/bin/miRNAture.pl -l /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/RFAM_14-4/rfam_models_list.txt -m OTHER_CM -pe 0 -nmodels /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/Other_CM /CMs -spe Lach -n_spe Latimeria_chalumnae -w /home/bioinf/Desktop/miRNAture/Tutorial/Code/../Results/miRNA_prediction -data /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data -cmp /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/Other_CM /CMs -hmmp /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/Other_HMM /HMMs -rep default,150,100
#Final searches
/home/bioinf/perl5/bin/miRNAture.pl -l /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data/RFAM_14-4/rfam_models_list.txt -m Final -pe 0 -spe Lach -n_spe Latimeria_chalumnae -w /home/bioinf/Desktop/miRNAture/Tutorial/Code/../Results/miRNA_prediction -data /home/bioinf/perl5/lib/perl5/auto/share/dist/Bio-miRNAture/Data -rep default,150,100
