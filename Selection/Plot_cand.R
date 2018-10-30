#!/usr/bin/env Rscript
#This script  is used in Find_cand_gene.sh
#plot statistics for candidate regions
#fetch detailed information for genes in the candidate regions

#setwd("E:/ubuntu/CB_VCF/Cand/Pi_union/")
VV <- commandArgs(TRUE)

source("/mnt/e/ubuntu/CB_VCF/Cand/Plot_SNP_fun.R")
#args=commandArgs(TRUE)
#load_in1()

Raw=read.table(VV,header=FALSE,stringsAsFactors = FALSE)
for(SNP in Raw[,1]){

  try(DA_rsb <- plot_XPEHH(SNP))
  
  try(Fst <- plot_Fst(SNP))
  
  try(Get_gene(SNP,GO_Ortho))
  
  try(plot_Pai(SNP))
  
  try(plot_hapflk(SNP))
}


