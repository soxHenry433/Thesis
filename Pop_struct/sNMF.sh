#!/bin/bash
#This script infers population structure by sNMF in R

#Convert the numeric format into .geno format
awk '
    {gsub("NA","9")
    gsub(1,2)
    gsub(0.5,1)
    gsub("\t","")
    print $0}' /mnt/e/ubuntu/SRR/GT/GT.num > GT.geno  

#perform sNMF with R package LEA
#Capture enropy output
R --vanilla --slave <<R_SCRIPT
library(LEA)
project.snmf <- snmf(input.file="GT.geno",
                     K=2:15,
     project= "continue",
     repetitions = 20,
     CPU = 8,
     alpha = 10, tolerance = 0.00001, entropy = TRUE,
     percentage = 0.05, I = 10000,
     iterations = 200, ploidy = 2, seed = -1
     )

AA=load.snmfProject("GT.snmfProject")
capture.output(for (i in 2:15) print(cross.entropy(AA,
                K=i)), file="entropy.output.txt")
R_SCRIPT

