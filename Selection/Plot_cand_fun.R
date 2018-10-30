#The functions used in Plot_cand.R 
library(pegas)
library(rehh)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(data.table)

fread_packed <- function(file_dir){
  as.data.frame(fread(file_dir,header=FALSE,sep="\t",skip=1),stringsAsFactors=FALSE)
}

GO_Ortho <- fread_packed("/mnt/e/ubuntu/CB_VCF/Gene_data/Cabbage_go_ortho.assoc")

#cat("Platform: ",.Platform$OS.type)

#Function for Plot_SNP.R

plot_XPEHH <- function(SNP){
  DA<- read.table(paste0("Tmp_",SNP,".rehh"),header=FALSE)
  names(DA) <- c("CHR","POSITION","iES1","iES2","Pv_Rsb","Pv_XPEHH")
  DA_rsb <- melt(DA[c("CHR","POSITION","Pv_XPEHH")],id.vars = c("CHR","POSITION"))
  DA_hh <- melt(DA[c("CHR","POSITION","iES1","iES2")],id.vars = c("CHR","POSITION"))
  
  adj_num <- (max(DA_rsb$value,na.rm = TRUE)*0.9)
  
  rr <- range(DA_hh$value,na.rm=TRUE)
  DA_hh$iES_adj <- ((DA_hh$value - rr[1])/diff(rr))*adj_num
  
  p_rsb <- ggplot() + geom_point(data=DA_rsb,aes(x=POSITION,y=value,color=variable)) + 
    geom_step(data=DA_hh,aes(x=POSITION,y=iES_adj,color=variable)) +
    labs(x="Position",y="-log10(p)") + 
    theme(
      axis.text.x.top = element_text(vjust=0.5, hjust=0),
      axis.line = element_line(color="black"),
      panel.background=element_blank(), 
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      axis.title.y = element_text(size=9),
      axis.title.x = element_text(size=10),
      axis.text.x=element_text( colour="grey50"), 
      axis.text.y=element_text( colour="grey50"), 
      axis.ticks.x = element_line(colour=NA) 
    ) + 
    scale_y_continuous(expand =c(0.02,0),sec.axis=sec_axis(~./adj_num*diff(rr) + rr[1],name = "iES")) + 
    scale_color_brewer(palette = "Set1")
    
  ggsave(paste0(SNP,"_rsb.png"),p_rsb,units = "cm",width=25,height = 12)
}


plot_Pai <- function(SNP){
  filn <- paste0("Tmp_",SNP,".Pi_",c("C","B"))
  Pi_C <- read.table(filn[1],header=FALSE,sep=" ")
  Pi_B <- read.table(filn[2],header=FALSE,sep=" ")
  filn <- paste0("Tmp_",SNP,"_tw.Pi_",c("C","B"))
  Pi_C_tw <- read.table(filn[1],header=FALSE,sep=" ")
  Pi_B_tw <- read.table(filn[2],header=FALSE,sep=" ")
  
  DA <- data.frame(Pos = c(Pi_B$V2,Pi_C$V2,Pi_B_tw$V2,Pi_C_tw$V2)/1e6,
                   Pai = c(Pi_B$V3,Pi_C$V3,Pi_B_tw$V3,Pi_C_tw$V3),
                   Variety = rep(c("Broccoli","Cauliflower","Broccoli","Cauliflower"),c(nrow(Pi_B),nrow(Pi_C),nrow(Pi_B_tw),nrow(Pi_C_tw))),
                   Population = rep(c("Cheng et al","Cheng et al","Taiwan","Taiwan"),c(nrow(Pi_B),nrow(Pi_C),nrow(Pi_B_tw),nrow(Pi_C_tw))))
  DA$Variety <- relevel(DA$Variety,"Cauliflower")
  p <- ggplot(DA) + geom_point(aes(x=Pos,y=Pai,color=Variety,shape=Population)) +  
    labs(title=SNP,y="Nucleotide diversity") +
    theme_classic() +
    theme(axis.title = element_text(size=15),
          legend.title = element_text(size=15),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size=11),
          legend.text=element_text(size=11),
          plot.title = element_text(size=16,hjust=0.5))
  
  ggsave(paste0(SNP,"_pi.png"),p,units = "cm",width=25,height = 12)
  
}


plot_Fst <- function(SNP){

  DA_fst <- read.table(paste0("Tmp_",SNP,".fst"),header=FALSE)
  colnames(DA_fst) <- c("Chr","Position","Fst")
  
  filn <- paste0("Tmp_",SNP,".Pi_",c("C","B"))
  Pi_C <- read.table(filn[1],header=FALSE,sep=" ")
  Pi_B <- read.table(filn[2],header=FALSE,sep=" ")
  
  DA <- merge(Pi_C[,1:3],Pi_B[,1:3],by=c("V1","V2"))
  DA_pi <- data.frame(Chr = DA$V1,
                    Position = DA$V2,
                    Pi_CB = DA[,3]/DA[,4],
                    Pi_BC = DA[,4]/DA[,3])
  pos <- as.numeric(substring(SNP,4))/1e6
  
  DA <- merge(DA_fst,DA_pi,by=c("Chr","Position"),all=TRUE)
  rr <- range(c(DA$Pi_CB,DA$Pi_BC),na.rm=TRUE)
  DA$Pi_CB_adj = (DA$Pi_CB - rr[1])/diff(rr)
  DA$Pi_BC_adj = (DA$Pi_BC - rr[1])/diff(rr)
  DA <- melt(DA,id.var=c("Position"))
  DA <- subset(DA,DA$variable %in% c("Fst","Pi_CB_adj","Pi_BC_adj"))
  DA$value <- as.numeric(DA$value)
  DA$variable <- as.factor(as.vector(DA$variable))
  DA$Position <- DA$Position/1e6
  p_fst <- ggplot(DA) + geom_point(aes(x=Position,y=value,color=variable)) + 
    theme_bw() + 
    geom_smooth(aes(x=Position,y=value),alpha=0.5,se=FALSE,method="loess",color="red") + 
    scale_y_continuous(sec.axis=sec_axis(~.*diff(rr) + rr[1],name = expression(pi~" ratio"))) + 
    labs(title=SNP,y="Fst",color="",x="Physical position (Mbp)") +
    theme(axis.title = element_text(size=15),
          legend.title = element_text(size=15),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size=11),
          legend.text=element_text(size=11),
          plot.title = element_text(size=16,hjust=0.5)) +
    geom_vline(aes(xintercept=pos),color="red",alpha=0.5,linetype="dotted") + 
    scale_color_discrete(labels=c(expression("FST"),expression(pi~" ratio (B/C)"),expression(pi~" ratio (C/B)"))) +
    theme(legend.text.align = 0)
  
  ggsave(paste0(SNP,"_fst.png"),p_fst,units = "cm",width=25,height = 12)
  
}

plot_hapflk <- function(SNP){
  filn1 <- paste0("Tmp_",SNP,".hapflk")
  filn2 <- paste0("Tmp_",SNP,"_tw.hapflk")
  Raw_ncbi <- read.table(filn1,header=FALSE,sep=" ")
  Raw_tw <- read.table(filn2,header=FALSE,sep=" ")
  
  DA <- data.frame(Pos=c(Raw_ncbi$V2,Raw_tw$V2)/1e6,
                   Hapflk = c(Raw_ncbi$V3,Raw_tw$V3),
                   Pvalue = -1*log10(c(Raw_ncbi$V4,Raw_tw$V4)),
                   Population =c(rep("Cheng et al",nrow(Raw_ncbi)),rep("Taiwan",nrow(Raw_tw))))
  pos <- as.numeric(substring(SNP,4))
  print(SNP)
  p_hap <- ggplot(DA) + geom_point(aes(x=Pos,y=Pvalue,color=Population)) + 
    theme_bw() + 
    labs(title=SNP) +
    theme(axis.title = element_text(size=15),
          legend.title = element_text(size=15),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size=11),
          legend.text=element_text(size=11),
          plot.title = element_text(size=16,hjust=0.5)) +
    #geom_vline(aes(xintercept=pos),color="red",alpha=0.5,linetype="dotted") +
    labs(x="Physical position (Mbp)", y="-log10(P-value)") 
    
    
  ggsave(paste0(SNP,"_hapflk.png"),p_hap,units = "cm",width=24,height = 10)
  
}


Get_gene <- function(SNP,GO_Ortho){
  filn <- paste0("Tmp_",SNP,".gff")
  Raw <- readLines(filn)
  s_Gff <- strsplit(Raw,"\t")
  Gff <- do.call(rbind,lapply(s_Gff,`[`,1:9))
  Gff <- subset(Gff,Gff[,3] == "gene")
  reg <- regexpr("(?<=ID=gene:)Bo[1-9]g[0-9]+",Gff[,9],perl=TRUE)
  Gene_list <- regmatches(Gff[,9],reg)
  Gene_list <- paste0(Gene_list,".1")
  
  go <- GO_Ortho[which(GO_Ortho$V2 %in% Gene_list),c(2,5,8,9)]
  reg <- gregexpr("(?<=AGI_LocusCode:)AT[1-9]G[0-9]+",go$V8,perl=TRUE)
  At_orth <- regmatches(go$V8,reg)
  At_orth <- lapply(At_orth,paste,collapse=',')
  go$At_orth <- do.call('c',At_orth)
  go <- go[,c(1,2,5)]
  write.table(go,paste0(SNP,"_ortho.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
}





