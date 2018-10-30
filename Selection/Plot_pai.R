#!/usr/bin/Rscript
#Plot ratio of nucleotide diversity

wd <- ifelse(.Platform$OS.type == "windows",
             "E:\\ubuntu\\CB_VCF\\TASSEL\\Result",
             "/mnt/e/ubuntu/CB_VCF/TASSEL/Result")
cat("Platform: ",.Platform$OS.type)
setwd(wd)

library(ggplot2)
library(reshape2)
library(data.table)

#Read and merge data
lf <- list.files(pattern='C[1-9].result',recursive=TRUE)
type <- sub('/.+',"",lf)
l <- lapply(lf, fread, sep="\t",header=TRUE)
dt <- rbindlist( l, idcol=TRUE)
rm(l)

#Create data frame
DA <- data.table(
  True_Position = dt[.id <=9, round((StartChrPosition + EndChrPosition)/2)],
  Broc_pi = dt[.id <=9, PiPerBP],
  Caul_dt = dt[.id >=10, PiPerBP],
  Chr = factor(dt[.id <=9, Chromosome])
)
DA[,Pi_ratio := log10(Caul_dt/Broc_pi)]
rm(dt)

#Calculate P value for each SNP
u_Pi_ratio <- (DA$Pi_ratio - median(DA$Pi_ratio,na.rm = TRUE))/sd(DA$Pi_ratio,na.rm = TRUE)
ff <- ecdf(u_Pi_ratio)
DA$p_v <- -1*log10(1-2*abs(ff(u_Pi_ratio)-0.5))



#Plot ratio of pi
DA_stat = DA[,.(Max_pos = max(True_Position),
                Mid_pos = median(True_Position),
                N = .N),by = Chr]
cum_base <- c(0,cumsum(DA_stat[["Max_pos"]])[-9])
DA[, Position := True_Position + rep(cum_base, DA_stat[["N"]])] 
DA_stat[,Mid_pos := Mid_pos + cum_base]

thred = -1*log10(0.05)
p1 <- ggplot(DA) + geom_point(aes(x=Position,y=p_v,color=Chr),alpha=0.4,size=0.9) + 
  labs(x="Chr",y=expression(pi)) + 
  theme(
    axis.line = element_line(color="black"),
    panel.background=element_blank(), 
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    axis.title.y = element_text(size=9),
    axis.title.x = element_text(size=10),
    axis.text.x=element_text(size=6, colour="grey50"), 
    axis.text.y=element_text(size=6, colour="grey50"), 
    axis.ticks.x = element_line(colour=NA) 
  ) + scale_x_continuous(expand =c(0,0),labels=paste0("C0",1:9),breaks = DA_stat[["Mid_pos"]]) + 
  guides(color = FALSE) + 
  scale_color_brewer(palette="Paired") + 
  geom_hline(aes(yintercept=thred),color="red",alpha=0.7,linetype="dotted")
ggsave(paste0("Manhattan_Pi_ratio.png"),p1,width=9.5,height=3,limitsize = FALSE)


fwrite(x = DA[,c("Chr","True_Position","Position","p_v")],
       file = "Pi_ratio_CB2250.txt",
       sep = "\t",
       eol = "\n",
       row.names = FALSE,
       col.names = TRUE,
       quote = FALSE)


#Box plot per parameter
dt2 <- dt[c("Type","PiPerBP","ThetaPerBP","TajimaD")]
dt2 <- melt(dt2)
library(ggthemes)
p <- ggplot(dt2) + geom_boxplot(aes(x=Type,y=value,fill=Type)) +
  theme_few() + scale_colour_few() + 
  guides(fill=FALSE) + 
  theme(axis.text.x = element_text(size=9),
        strip.text = element_text(size=13)) +
  facet_wrap(~variable,nrow=3,scales="free_y") +
  labs(x="Chr",y="value")

png("Boxplot.png",res=400,width=3600,height = 2600)
p
dev.off()




