#!/usr/bin/env Rscript
#Clusters hundreds of SNP into regions

#Read in results of selection scan and extract significant SNP
library(data.table)
Raw <- as.data.frame(fread("file:///E:/ubuntu/CB_VCF/TASSEL/Pi_pv.txt",header=TRUE),StringAsfactors=FALSE)
thred <- -1*log10(0.01)
sig_II <- which(Raw$log_pv > thred)
snp=Raw$SNP[sig_II]

Raw <- as.data.frame(fread("file:///E:/ubuntu/CB_VCF/hapflk/GT.hapflk_sc",header=TRUE),stringsAsFactors = FALSE)
hflk <- which(0.01 > Raw$pvalue)
snp=Raw$rs[hflk]
Raw$SNP <- Raw$rs
Raw$log_pv <- -1*log10(Raw$pvalue)

#Make all SNP into table
SNP_info <- data.frame(SNP=snp,
                       chr=substr(snp,2,2),
                       pos=as.numeric(substring(snp,4)))
Info_s <- split(SNP_info,SNP_info[,2])
Info_s <- lapply(Info_s,function(x){
  II=order(x[,3],decreasing=FALSE)
  x[II,]
})

#Cluster SNP with win size 250k into regions
win_size = 250000
AAA <- matrix(,nrow=0,ncol=4)
for(i in 1:length(Info_s)){
  da <- Info_s[[i]]
  AA <- c()
  while(nrow(da) > 0){
    II <- which(abs(da[1,3] - da[,3]) < win_size) 
    dd <- abs(da[II[length(II)],3] - da[-II,3])
    while(any(dd < win_size)){
      II2 <- which(dd < win_size)
      II <- c(II,length(II) + II2)
      dd <- abs(da[II[length(II)],3] - da[-II,3])
    }
    gr <- da[II,]
    da <- da[-II,]
    ll = round(range(gr[,3])[1] - win_size)
    uu = round(range(gr[,3])[2] + win_size)
    if(ll < 0 ) ll <- 0
    if(uu < 0 ) uu <- 0
    ii <- which(Raw$SNP %in% gr$SNP) 
    ii <- which.max(Raw$log_pv[ii])
    AA <- rbind(AA,c(as.vector(gr$SNP[ii]),as.vector(gr[1,2]),ll,uu))
  }
  AAA <- rbind(AAA,AA)
}

#Output
write.table(AAA,file="Range.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

