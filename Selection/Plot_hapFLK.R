#!/usr/bin/Rscript
#Plot hapFLK

library(data.table)

#Read data
wd <- ifelse(.Platform$OS.type == "windows",
             "E:\\ubuntu\\CB_VCF\\hapflk",
             ".")
cat("Platform: ",.Platform$OS.type,"\n")
setwd(wd)
FLK <- as.data.frame(fread("GT.flk"),StringAsfactor=FALSE)
HapFLK <- as.data.frame(fread("GT.hapflk_sc",header=TRUE),stringsAsFactors = FALSE)

names(FLK) <- c("snp","chr","pos","p0","flk","pvalue")
names(HapFLK) <- c("snp","chr","pos","hapflk","hapflk_scaled","pvalue")

Last_base <- aggregate(pos~chr ,data=FLK,max)
cum_base <- c(0,cumsum(Last_base[,2])[-9])
tt <- table(FLK$chr)

DA <- data.frame(Chr = as.factor(FLK$chr),
                 Pos = FLK$pos + rep(cum_base,as.numeric(tt)),
                 FLK = FLK$flk,
                 hapFLK = -1*log10(HapFLK$pvalue))
Midbase <- aggregate(Pos~Chr,data=DA,median)

library(ggplot2)
library(cowplot)

p1 <- ggplot(DA) + geom_point(aes(x=Pos,y=FLK,color=Chr),alpha=0.3,size=0.8) + 
  labs(y="FLK") + 
  scale_y_reverse(expand =c(0.02,0)) + 
  theme(
    axis.text.x.top = element_text(vjust=0.5, hjust=0),
    axis.line = element_line(color="black"),
    panel.background=element_blank(), 
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    axis.title.y = element_text(size=9),
    axis.title.x = element_text(size=10),
    axis.text.x=element_text(size=6, colour="grey50"), 
    axis.text.y=element_text(size=6, colour="grey50"), 
    axis.ticks.x = element_line(colour=NA)  
  ) + scale_x_continuous(position = "top",expand =c(0,0),labels=paste0("C0",1:9),breaks = Midbase[,2]) +
  scale_color_brewer(palette="Paired") + 
  guides(color = FALSE) #+  scale_colour_manual(values = rep(c("blue4", "orange3")))
p2 <- ggplot(DA) + geom_point(aes(x=Pos,y=hapFLK,color=Chr),alpha=0.3,size=0.8) + 
  labs(y="-log10(P-value)", x = "Physical position") + 
  scale_y_continuous(expand =c(0.02,0)) + 
  theme(
    axis.line = element_line(color="black"),
    panel.background=element_blank(), 
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    axis.title.y = element_text(size=10),
    axis.title.x = element_text(size=12),
    axis.text.x=element_text(size=8, colour="grey50"), 
    axis.text.y=element_text(size=8, colour="grey50"), 
    axis.ticks.x = element_line(colour=NA)
  ) + scale_x_continuous(expand =c(0,0),labels=paste0("C0",1:9),breaks = Midbase[,2]) +
  scale_color_brewer(palette="Paired") + 
  geom_hline(aes(yintercept=-1*log10(0.005)),color="red",linetype="dotted",alpha=0.7) +
  guides(color = FALSE)# +  scale_colour_manual(values = rep(c("blue4", "orange3")))

aa <- ggdraw() + 
  draw_plot(p1,x=0 ,y=0.5 ,width = 1, height = 0.5) + 
  draw_plot(p2,x=0 ,y=0 ,width = 1, height = 0.5) + 
  draw_plot_label(label = c("a", "b"), size = 10,
                  y = c(1, 0.55), x = c(0,0))


ggsave("Manhattan_FLK.png",aa,width=9.5,height=3,limitsize = FALSE)
ggsave("Manhattan_hapFLK.png",p2,width=9.5,height=2.3,limitsize = FALSE)

thre<- round(nrow(HapFLK)*0.001)
hh <- order(HapFLK$hapflk,decreasing=TRUE)[1:thre]

AA <- HapFLK[hh,]


write.table(AA,"Sig_hapflk.snp",quote=FALSE,col.names = TRUE,row.names=FALSE)



