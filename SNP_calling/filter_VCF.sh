#!/bin/bash
#This file filter the VCF file

ll=$#
a=0
while [ $a -lt $ll ]; do
    case "$1" in
            -*)
                Arg="$Arg $1 $2";
                shift;
                shift;;
    esac
    a=$((a+1))
done

echo "Arguments: $Arg"

if [ -z ${1:+x} ];then
    echo "No file specified"
    exit 1
else
    echo "Files: $@"
fi

if [ "${1: -2}" != "gz" ];then
    echo "No compress file specified"
    exit 1
fi

if [ ! -f $1".tbi" ];then
    echo "No index file exists"
    exit 1
fi

bcf="/mnt/e/ubuntu/bcftools/bcftools"
cd /mnt/g/VCF

Filter_vcf(){
    $bcf view --types snps \
        -i 'Info/AN >= 250' \
        --min-alleles 2 \
        --min-af 0.01:minor \
        -O v \
        -o ${1%.vcf.gz}_filter.vcf \
        $thread \
        $Arg \
        $1 || exit 1
    sed -i '/##contig=<ID=Scaffold/d' ${1%.vcf.gz}_filter.vcf
    echo "${1%.vcf.gz}_filter.vcf"
}

if [ $# -eq 1 ];
then
    thread="--threads 1"
    Filter_vcf $1 
else
    thread="--threads 8"
    export thread Arg bcf
    export -f Filter_vcf
    parallel Filter_vcf ::: "$@" 
fi

