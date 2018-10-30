#!/bin/bash
#This script construct phylogenetic tree with SNPphylo

#Filter data 
vcf="/mnt/e/ubuntu/SRR/GT/GT_SRR.vcf.gz"
plink --vcf $vcf --maf 0.1 --geno 0.5 --recode vcf --out Tmp --allow-extra-chr
sed -i 's/^C//' Tmp.vcf 

#Plot tree 
SNP_tree="/mnt/e/ubuntu/SRR/SNPhylo/snphylo.sh"
$SNP_tree -v Tmp.vcf -M 0.5 -b 100 -P DATA/GT_bs -m 0.1 -r


