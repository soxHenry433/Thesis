#!/bin/bash
#This script call the phased SNP for cauliflower and broccoli from CH dataset


GATK="java -jar -Xmx16g /lop/GATK.jar "
picard="java -Xmx16g -jar /lop/picard.jar"
Index_bt="/mnt/c/Users/HP/Documents/TD/Bowtie2_index/TO1000_chr"
bam_dir="/mnt/g/BAM/BAM"
Var="/lop/SRR/GT/Var.table"

TD="/mnt/c/Users/HP/Documents/TD"
SD="/mnt/g/HaploCaller"

#First, arrange bam into format that can be processed by gatk
#Second, call these bam files by GATK haplotype caller
call_bam(){
    accession=$1
    cd $TD
    cp $bam_dir/${accession}.bam .
    $picard SortSam \
        INPUT=${accession}.bam \
        OUTPUT=${accession}2.bam \
        SORT_ORDER=coordinate
    $picard AddOrReplaceReadGroups \
        INPUT=${accession}2.bam \
        OUTPUT=${accession}.bam \
        CREATE_INDEX=1 \
        RGID=$accession \
        RGLB=library_$accession \
        RGPL=illumina \
        RGPU=$accession \
        RGSM=$accession 
    $picard BuildBamIndex INPUT=${accession}.bam
    rm ${accession}2.bam
    $GATK \
        -R $Index_bt".fa" \
        -T HaplotypeCaller \
        -I ${accession}.bam \
        -stand_call_conf 30 \
        --emitRefConfidence GVCF \
        -o ${accession}.g.vcf 
    mv ${accession}.g.vcf $SD/gVCF/${accession}.g.vcf
    mv ${accession}.g.vcf.idx $SD/gVCF/${accession}.g.vcf.idx
    mv ${accession}.bam $SD/BAM/${accession}.bam
    mv ${accession}.bai $SD/BAM/${accession}.bai
    echo "${accession} complete"
}

#Phased SNP calling
for i in $(ls $bam_dir/SRR*.bam | cut -d'.' -f1);do
    call_bam $id
done 

#create variable All_sample, which is used in GATK  CombineGVCFs 
All_sample=""
for i in $(ls SRR*.g.vcf | cut -d'.' -f1);do
    vcf_name=$i'.g.vcf'
    All_sample="$All_sample --variant $vcf_name "
done

#Combine all gvcf of different sample into one gvcf
#The process is carried out by chromosome
for i in {1..9};do
    $GATK \
        -T CombineGVCFs \
        -R $Index_bt".fa" \
        $All_sample \
        --intervals C${i} \
        -o $TD/C${i}.g.vcf
    mv $TD/C${i}.g.vcf /mnt/g/HaploCaller/C${i}.g.vcf
    mv $TD/C${i}.g.vcf.idx /mnt/g/HaploCaller/C${i}.g.vcf.idx
done

#For gVCF for rach chromosome
#First, filter gVCFs by bcftools and convert them into VCF
#Sec, Convert the file into tped format
#Third, devide tped into cauliflwoer and broccoli
for i in {1..9};do
    bcftools view  \
        -i 'Info/AN >= 81 & QUAL>=40' \
        -m2 -M2 -v snps \
        -O v \
        -o C${i}_f.vcf \
        $SD/VCF/C${i}_SRR.vcf 

    mawk 'BEGIN{OFS="\t"}
            !/#/{
            gsub($9,"GT",$9)
            for(i=10;i<=NF;i++){
                match($i, /[01]\|[01]/);
                if(RSTART==0){
                    $i=substr($i, 1, 3)
                }else{
                    $i=substr($i, RSTART, RLENGTH)
                }
            }
        }{  print $0 }' C${i}_f.vcf > C${i}_f_p.vcf

    rm C${i}_f.vcf 
    
    for j in 0 1;do
        VV=$(grep -n " $j" Popfile.txt | cut -d':' -f1 )
        COL=$(sed 's/ /,/g' <(echo $VV))
        grep -v "#" C${i}_f_p.vcf | \
        cut -f10-52 | \
        cut -f$COL | \
        awk '{
            gsub(/\/|\|/," ")
            gsub(/0/,"2");
            gsub(/\./,"0");
            gsub(/\t/," ")
            print $0 }' > C${i}_$j.tped
    done
    grep "^[^#]" C${i}_f_p.vcf  | awk '{print $1"_"$2,$1,$2,$4,$5}' > C${i}.MAP

    mv C${i}_f_p.vcf /lop/CB_VCF/GT/C${i}_f_p.vcf
    mv C${i}.MAP /lop/CB_VCF/rehh/C${i}.MAP
    mv C${i}_*.tped /lop/CB_VCF/rehh/

done
















