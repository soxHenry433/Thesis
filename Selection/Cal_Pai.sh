#!/bin/bash
#Calculate genome-wide nucleotide diversity

tassel="/mnt/e/ubuntu/Tassel/run_pipeline.pl -Xmx16g"
vcf="/mnt/e/ubuntu/CB_VCF/Het/Red_result/GT_f_rep.vcf.gz"
Var="/mnt/e/ubuntu/CB_VCF/GT/Popfile.txt"

#First filter and devide VCF into 9 chromosome
[ ! -d VCF_chr ] && mkdir VCF_chr
if [ -z "$(ls -A VCF_chr)" ]; then
    cd VCF_chr
    awk '{print "C"$1,$1}' <(seq 1 9) > Chr.txt
    makegoodvcf(){
        bcftools view $vcf \
            -r "C$1" \
            --min-af 0.01 \
            -O v \
            -o Tmp_ftr_$1.vcf
        bcftools annotate --rename-chrs Chr.txt -O v -o Tmp_$1.vcf Tmp_ftr_$1.vcf 
        rm Tmp_ftr_$1.vcf
    }
    export vcf
    export -f makegoodvcf
    parallel makegoodvcf ::: {1..9}
    cd ..
fi

WinSize=2250
[ ! -d Result ] && mkdir Result

#Calculate diversity 
Cal_div(){
    for i in 0 1; do
    # 0 is broccoli 
    # 1 is cauliflower
        for k in {1..9}; do
            #Filter VCF into broccoli or cualiflower
            bcftools view \
                -S <(awk -v GROUP="$i" '$3==GROUP{print $1}' $Var) \
                -O v -o Result/Tmp.vcf \
                VCF_chr/Tmp_$k.vcf

            case "$i" in
            "1")
                GROUP="Cauliflower";
                ;;
            "0")
                GROUP="Broccoli";
                ;;
            *)
                echo "Invalid"
                ;;
            esac

            #Calculate diversity by TASSEL
            $tassel -vcf Result/Tmp.vcf \
                -diversity \
                -diversitySlidingWin \
                -diversitySlidingWinStep 1 \
                -diversitySlidingWinSize $WinSize \
                -export Result/${GROUP}_C${k}_ \
                -exportType Table
            
            #Revise the resulting file in order to easier data processing 
            awk 'BEGIN{OFS="\t"}
            NR==1{print "Chr",$0}
            FNR>1{
            gsub("[A-Za-z]+_|_[0-9]*1.txt","",FILENAME);
            gsub(",","",$0);
            print FILENAME,$0}' Result/${GROUP}_C${k}_1.txt > Result/${GROUP}_C${k}.result
            
            rm  Result/*.txt
        done
    rm Result/Tmp.vcf
    done
}

Cal_div


