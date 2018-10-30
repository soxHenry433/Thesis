#!/usr/bin/env Rscript
#Calculate XP-EHH 

args = commandArgs(trailingOnly=TRUE)
setwd("C:/Users/user/Desktop/rehh/")
library(rehh)
library(reshape2)

#Calculate XP-EHH by chromosomes
for(i in paste0("C",1:9)){
  #Parse argument
  args <- paste0(i,c(
    "_0.tped",
    "_1.tped",
    ".MAP",
    ".rehh"
  ))
  
  #read in file
  cat("\nmin_perc_geno.hap=50\nmin_perc_geno.snp=95\n")
  hap1 <- data2haplohh(hap_file=args[1],
                    map_file=args[3],
                    haplotype.in.columns=TRUE,
                    min_perc_geno.hap=50,
                    min_perc_geno.snp=95,
                    recode.allele=FALSE)
  
  hap2<-data2haplohh(hap_file=args[2],
                    map_file=args[3],
                    haplotype.in.columns=TRUE,
                    min_perc_geno.hap=50,
                    min_perc_geno.snp=95,
                    recode.allele=FALSE)
  
  #Scan for iHS
  hh1 <- scan_hh(hap1,threads=7)
  hh2 <- scan_hh(hap2,threads=7)
  nn <- intersect(rownames(hh1),rownames(hh2))
  hh1 <- hh1[nn,]
  hh2 <- hh2[nn,]
  
  #Calculate xp-ehh
  rsb <- ies2rsb(hh1,hh2,"Broc","Caul")
  xpehh <- ies2xpehh(hh1,hh2,"Broc","Caul")
  DA_rsb <- merge(rsb,xpehh,by=c("CHR","POSITION"),all=TRUE)
  DA_hh <- merge(hh1,hh2,by=c("CHR","POSITION"),all=TRUE)
  DA <- merge(DA_rsb,DA_hh,by=c("CHR","POSITION"),all=TRUE)
  
  write.table(DA,args[4],col.names = TRUE,row.names = FALSE,sep="\t",quote=FALSE)
}
