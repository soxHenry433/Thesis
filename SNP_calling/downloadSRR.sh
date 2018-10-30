#!/bin/bash
#This program is to download SRA file from NCBI according to SraRunTable.txt

SRR=$(awk -F'\t' '/oler/{print $10}' SraRunTable.txt)

#Download using ascp recommended by NCBI
for i in $SRR; do
    ascp $i;
done;

#scp SSR to the server 
cd /mnt/g/SRA
SRA=`ls *.sra`
tgt="lab_311@140.112.74.31:~/lab_311/Sox/DATA/SRA/"
scp $SRA $tgt
