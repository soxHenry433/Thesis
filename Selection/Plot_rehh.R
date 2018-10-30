#!/usr/bin/env Rscript
#For plot genomewide manhattan plot

lf <- grep("\\.rehh",list.files(),value = TRUE)

library(data.table)

l <- lapply(lf, fread, sep="\t",header=TRUE)
dt <- rbindlist( l )
rm(l)

DA <- as.data.frame(dt[,c(1,2,3,5,10,15)])
names(DA) <- c("CHR","POSITION","Rsb","XPEHH","iES_Broc","iES_Caul")

Last_base <- aggregate(POSITION~CHR,data=DA,max)
cum_base <- c(0,cumsum(Last_base[,2])[-9])
tt <- table(DA$CHR)

DA$pos <- DA$POSITION + rep(cum_base,as.numeric(tt))
Midbase <- aggregate(pos~CHR,data=DA,median)

DA$CHR <- as.factor(DA$CHR)

adj_Pv <- function(x){
  stopifnot(length(x) > 1)
  ff <- ecdf(x)
  -1*log10(1-2*abs(ff(x)-0.5))
}

DA$P_Rsb = adj_Pv(DA$Rsb)
DA$P_XPEHH = adj_Pv(DA$XPEHH)

write.table(DA,"CB_result.txt", quote=FALSE, row.names=FALSE,sep="\t")

library(ggplot2)
library(cowplot)

p1 <- ggplot(DA) + geom_point(aes(x=pos,y=P_Rsb,color=CHR),alpha=0.3,size=0.8) + 
  labs(x="Rsb",y="-log10(p)") + 
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
  guides(color = FALSE) #+  scale_colour_manual(values = rep(c("blue4", "orange3")))
p2 <- ggplot(DA) + geom_point(aes(x=pos,y=P_XPEHH,color=CHR),alpha=0.3,size=0.8) + 
  labs(x="XPEHH",y="-log10(p)") + 
  scale_y_continuous(expand =c(0.02,0)) + 
  theme(
    axis.line = element_line(color="black"),
    panel.background=element_blank(), 
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size=10),
    axis.text.x=element_text(size=9, colour="grey50"), 
    axis.text.y=element_text(size=9, colour="grey50"), 
    axis.ticks.x = element_line(colour=NA)
  ) + scale_x_continuous(expand =c(0,0),labels=paste0("C0",1:9),breaks = Midbase[,2]) +
  scale_color_brewer(palette="Paired") + 
  guides(color = FALSE)# +  scale_colour_manual(values = rep(c("blue4", "orange3")))

aa <- ggdraw() + 
  draw_plot(p1,x=0 ,y=0.5 ,width = 1, height = 0.5) + 
  draw_plot(p2,x=0 ,y=0 ,width = 1, height = 0.5) + 
  draw_plot_label(label = c("a", "b"), size = 10,
                  y = c(1, 0.55), x = c(0,0))


ggsave("Manhattan_XPEHH.png",p2,width=9.5,height=2.3,limitsize = FALSE)


