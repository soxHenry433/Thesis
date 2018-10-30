#!/bin/bash
#Calculate hapFLK

GT="/mnt/e/ubuntu/CB_VCF/Het/Red_result/GT_f_rep.vcf.gz"
Var="/mnt/e/ubuntu/CB_VCF/GT/Popfile.txt"
hk="/usr/local/bin/hapflk"

#Convert vcf into ped
#Assign subpopulation to each accession
#Calculate hapFLK by chromosomes
for i in {1..9};do
    Chr="C$i"
    bcftools view $GT | awk -v CC=$Chr '/#/{print}$1==CC{print}' > Tmp.vcf
    plink --vcf Tmp.vcf --recode --out GT --allow-extra-chr
    awk -F' ' 'BEGIN{OFS="\t"}{$1="";sub(" ","");$3=-9;$4=-9;$5=-9;print $0}' GT.ped > Mytmp.ped 
    paste -d' ' <(cut -d' ' -f3 $Var) Mytmp.ped > GT.ped && rm  Mytmp.ped
    awk 'BEGIN{OFS="\t"}{print $1,$1"_"$4,-9,$4}' GT.map > myGT.map
    mv myGT.map  GT.map
    $hk --file GT -K 10 --nfit=10 --ncpu 4 
    mv hapflk.hapflk ${Chr}.hapflk
    mv hapflk.flk ${Chr}.flk
    mv hapflk.frq ${Chr}.frq
done

#Conbine results of 9 chromosome into one
cat $(ls C*.flk) > GT.flk
rm  C*.flk
rg "^C" --no-line-number GT.flk > myGT.flk
mv myGT.flk GT.flk
cat $(ls C*.hapflk) > GT.hapflk
rm  C*.hapflk
rg "^C" --no-line-number GT.hapflk > myGT.hapflk
mv myGT.hapflk GT.hapflk
cat $(ls C*.frq) > GT.frq
rm  C*.frq

#Convert hapFLK into p value
#The script is provided by the author of hapFLK
sed -i '1i rs chr pos hapflk' GT.hapflk 
/mnt/e/ubuntu/hapflk-1.4/utils/scaling_chi2_hapflk.py GT.hapflk 6 2

#Test for best k parameter that is used in hapFLK
fp="/mnt/e/ubuntu/fastPHASE"
plink --vcf $GT --recode fastphase --out GT --allow-extra-chr
sed '/P/d' GT.chr-C1.recode.phase.inp > GT_C1.inp
rm GT.chr-C*.recode.phase.inp
[ ! -d Test_k ] &&　mkdir Test_k
for i in 3 5 10 15 20 25 30 35 40;do
    $fp -K$i -T10 -C25 -K10 -H-1 -Z -oC1.k$i GT_C1.inp
done

/mnt/e/ubuntu/hapflk-1.4/utils/scaling_chi2_hapflk.py hapflk.hapflk 3 2


