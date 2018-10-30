#!/bin/bash
#This program maps the reads to the reference
#From fastq to bam 

SRR=$1 
USE_TEMP_DIR=$2

#Export to the disk with faster IO (such as SSD)
SD="/mnt/g"  #Main store dir
TD=$SD 
TD="/mnt/c/Users/HP/Documents/TD"
if [ $USE_TEMP_DIR = T ]; then
    TD=$3 #Temp dir
    if [ -z ${3+x} ]; then
       TD="/mnt/c/Users/HP/Documents/TD/" 
    fi
fi


#Dump the fq 
SRAtool="/mnt/e/ubuntu/SRAtool/bin"
fq1="${TD}/FQ/${SRR}_1.fastq"
fq2="${TD}/FQ/${SRR}_2.fastq"
if [ ! -f $fq2 ] || [ ! -f $fq1 ];then
    $SRAtool/fastq-dump ${SD}/SRA/$SRR".sra" --split-3 -O ${TD}/FQ/
fi
echo "Dump Finished"
#ll=`head $fq1 -n 1 | awk -F' ' '{split($3,A,"=");print A[2]-10}'`
#if [ -z $ll ]; then
#    exits 1
#fi


#Mapping by bowtie2
bowtie="/mnt/e/ubuntu/bowtie2-2.3.3.1/bowtie2"
Index_bt="/mnt/g/Bowtie2_index/TO1000_chr"
sam="${TD}/BAM/${SRR}.sam"
$bowtie -p 8 --trim5 4 --trim3 8 -x $Index_bt -1 $fq1 -2 $fq2 -S $sam || exit


#Picard and samtools filter and indexing
picard="java -Xmx16g -jar /lop/picard.jar"
bam="${TD}/BAM/${SRR}.bam"
bam2="${TD}/BAM/${SRR}_2.bam"
MM="${TD}/BAM/metrics_$SRR.txt"
samtools view -q 20 -b -F 0x904 $sam > $bam2 || exit
$picard SortSam  INPUT=$bam2  OUTPUT=$bam SORT_ORDER=coordinate
$picard BuildBamIndex INPUT=$bam || exit

#Move back to the original folder
if [ $USE_TEMP_DIR = T ]; then
    mv $bam ${SD}/BAM/${SRR}.bam
    mv ${TD}/BAM/${SRR}.bai ${SD}/BAM/${SRR}.bai
    rm $fq1 $fq2 $sam $bam2
fi

