#!/usr/bin
#This is the script to calculate LD by plink

Var="/mnt/e/ubuntu/SRR/GT/Var.table"

#Calcutate LD for each sub-population
for i in Cauliflower Broccoli Cabbage Kohlrabi;do
    awk -v type=$i -F'\t' 'BEGIN{OFS=" "}$2~type{print $1,$1,$3}' $Var > Popfile
    plink --bfile ../GT/GT \
        --keep Popfile \
        --geno 0.1 \
        --maf 0.05 \
        --ld-window-kb 100 \
        --ld-window 10000 \
        --r2 \
        --out $i \
        --allow-extra-chr
    for j in {1..9};do 
        doit $j $i".ld"
    done
    rm $i".ld"
done

#Calcutate LD for all 119 accessions
plink --bfile ../GT/GT \
    --geno 0.1 \
    --maf 0.05 \
    --r2 \
    --ld-window-kb 100 \
    --ld-window 100000 \
    --out All \
    --ld-window-r2 0 \
    --allow-extra-chr

#Calculate the mean of LD every 100 kb
doit(){
    awk -v Chr="C$1" '$1~Chr{D=$5-$2;print D,$7}' $2 > Tmp
    sort -n -t' ' -k1 Tmp > Tmp_s
    awk 'BEGIN{
            step=50;
            ss=0;
            Sum=0;
            LL=0
        }
        {
            while($1>=ss && $1<(ss+50)){
                Sum=Sum+$2;
                LL=LL+1;
                next
            }
            if(LL>0){
                Mean=Sum/LL;
                print ss,Mean,LL
            }
            ss=ss+50;
            Sum=0;
            LL=0
        }' Tmp_s > ${2%.ld}_C"$1"_50windows
    rm Tmp Tmp_s
}
export -f doit

for i in {1..9};do 
    doit $i 
done
