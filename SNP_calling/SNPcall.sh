#!/bin/bash

if [ -z ${1:+x} ];then 
    echo "No chromosome specified"
    exit 1
fi

TD="/mnt/c/Users/HP/Documents/TD"
cd $TD

#Creating Bowtie2 index for the reference
Index_bt="$TD/Bowtie2_index/TO1000_chr"
if [ ! -f "$Index_bt.dict" ] || [ ! -f "$Index_bt.fa.fai" ];
then
    if [ ! -f $Index_bt".fa" ];then 
        exit 1
    fi
    echo "Creating index..."
    samtools faidx $Index_bt".fa"
    picard="java -Xmx16g -jar /lop/picard.jar"
    $picard CreateSequenceDictionary R= $Index_bt".fa" O= $Index_bt".dict"
fi

if [ ! -d $1 ]
then
    echo "No BAM folder"
    exit 1
fi
cd $1 

if [ ! -e bamlist ];then 
    echo "No bamlist, creating one..."
    ls *bam > bamlist || exit 1
fi

#Calling SNP by bcftools
bcf="/mnt/e/ubuntu/bcftools/bcftools"
time samtools mpileup -b bamlist -uf $Index_bt".fa" -E -t DP | \
$bcf call --variants-only \
 --skip-variants indels \
 --multiallelic-caller \
 --threads 8 \
 -O v \
 -o ../GT_$1.vcf || exit 1

echo "Start to compress"
bgzip ../GT_$1.vcf 
tabix ../GT_$1.vcf.gz

