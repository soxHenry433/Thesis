#!/bin/bash
#This script does the following things:
#Fetch statistics in each region
#Find candidate genes in each candidate regions
#Plot the statistics and genotypes of these regions

vcf="/mnt/e/ubuntu/CB_VCF/Het/Red_result/GT_f_rep.vcf.gz"
Var="/mnt/e/ubuntu/CB_VCF/GT/Popfile.txt"
gff="/mnt/e/ubuntu/SRR/GT/Cabbage.gff"

Fst_file="/mnt/e/ubuntu/CB_VCF/Fst/PAIR_DIST_1_0.txt"
Pi_folder="/mnt/e/ubuntu/CB_VCF/TASSEL/Result/"
CB_rehh="/mnt/e/ubuntu/CB_VCF/rehh/CB_result.txt"
hapflk="/mnt/e/ubuntu/CB_VCF/hapflk/GT.hapflk_sc"
hapflk_tw="/mnt/e/ubuntu/All_GT/hapflk/CB/hapflk.hapflk_sc"
TW_Pai_B="/mnt/e/ubuntu/All_GT/TASSEL/Broccoli_pop_15_result"
TW_Pai_C="/mnt/e/ubuntu/All_GT/TASSEL/Cauliflower_pop_15_result"

awk '$3~/1/{print $1,$1 > "Tmp1"}
    $3~/0/{print $1,$1 > "Tmp0"}' $Var

do_it(){
    set -e
    VV="Tmp_$1"
    bcftools view $vcf \
        -O v \
        -r C$2:$3-$4 \
        -o $VV.vcf  
    awk -v Chr="C$2" -v Pos1="$3" -v Pos2="$4" \
        '$1==Chr && $3=="gene"{if($4 >= Pos1 && $5 <= Pos2){print}}' $gff > $VV.gff    
    awk -v Pos1="$3" -v Pos2="$4"  \
        'NR==1{
            next
        }{
            sub("Result/","",$1)
            P=($4+$5)/2
            if(P >= Pos1 && P <= Pos2){
                print $1,P,$12}}' ${Pi_folder}Broccoli_C${2}.result > $VV.Pi_B
    awk -v Pos1="$3" -v Pos2="$4"  \
        'NR==1{
            next
        }{
            sub("Result/","",$1)
            P=($4+$5)/2
            if(P >= Pos1 && P <= Pos2){
                print $1,P,$12}}' ${Pi_folder}Cauliflower_C${2}.result > $VV.Pi_C
    awk -v Chr="C$2" -v Pos1="$3" -v Pos2="$4" \
        '$1==Chr{if($2 >= Pos1 && $2 <= Pos2){
            print $1,$2,$5,$6,$8,$9}}' $CB_rehh > $VV.rehh
    awk -v Chr="$2" -v Pos1="$3" -v Pos2="$4" \
        '$1==Chr{if($2 >= Pos1 && $2 <= Pos2){
            print $1,$2,$8}}' $Fst_file > $VV.fst
    awk  -v Chr="C$2" -v Pos1="$3" -v Pos2="$4" -F' ' \
        '$2==Chr{if($3 >= Pos1 && $3 <= Pos2){
        print $2,$3,$5,$6}}' $hapflk > $VV.hapflk
    awk  -v Chr="C$2" -v Pos1="$3" -v Pos2="$4" -F' ' \
        '$2==Chr{if($3 >= Pos1 && $3 <= Pos2){
        print $2,$3,$5,$6}}' $hapflk_tw > ${VV}_tw.hapflk
    awk  -v Chr="$2" -v Pos1="$3" -v Pos2="$4" \
        '$3==Chr{
            P=($4+$5)/2
            if(P >= Pos1 && P <= Pos2){
                print $1,P,$12}}' $TW_Pai_B > ${VV}_tw.Pi_B
    awk  -v Chr="$2" -v Pos1="$3" -v Pos2="$4" \
        '$3==Chr{
            P=($4+$5)/2
            if(P >= Pos1 && P <= Pos2){
                print $1,P,$12}}' $TW_Pai_C > ${VV}_tw.Pi_C

     #Plot gene stats
    /lop/CB_VCF/Cand/Get_genestat.R $VV
    #Plot genotypes
    /lop/CB_VCF/GT/Plot_GT.sh -c "C$2" -p $3"-"$4 -o $1 -r 0 -d 1 
    }

    #../Plot_SNP.R ${VV}_0.tped ${VV}_1.tped Tmp_$1.MAP $1_01
    #../Plot_SNP.R ${VV}_1.tped ${VV}_0.tped Tmp_$1.MAP $1_10
    #rm Temp_*



export vcf gff CB_rehh Fst_file Pi_folder hapflk hapflk_tw TW_Pai_B TW_Pai_C
export -f do_it 


parallel --colsep '\t' do_it :::: $1 

#Plot_SNP.R Plot statistics and fetch detailed information for all genes 
../Plot_SNP.R $1 

rm Tmp_**




