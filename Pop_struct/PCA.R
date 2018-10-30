#/usr/bin/env Rscript
#This is the script to perform PCA 

setwd("E:/ubuntu/SRR/PCA")
library(ggplot2)
library(data.table)

#Read data
GT <- as.data.frame(fread("E:/ubuntu/SRR/GT/GT_num_imp.txt"))
SNPval <- as.matrix(GT)
rm(GT)

#Scale data
#Not use scale function in the consideration of memory
for(i in 1:nrow(SNPval)){
  x <- SNPval[i,]
  x <- x[!is.na(x)]
  ss <- sqrt(sum(x^2)/max(1, length(x) - 1L))
  SNPval[i,]<- (x - mean(x))/ss
}

#Round and sample SNP for memory concern
SNPval <- round(SNPval,2)
to_be_rm <- sample(nrow(SNPval),7e5)
SNPval <- SNPval[-to_be_rm,]
rm(to_be_rm)

#Perform PCA with SVD function
svx <- svd(t(SNPval))
svx.pca <- t(SNPval) %*% svx$v
Varaince <- svx$d^2/sum(svx$d^2)*100

#Plot PCs X variance
DA <- data.frame(Index=1:length(svx$d),
                Varaince=Varaince,
                Cum_Varaince=cumsum(Varaince))
DA$Cum_v_adj <- (DA$Cum_Varaince - min(DA$Cum_Varaince))/diff(range(DA$Cum_Varaince)) * diff(range(Varaince)) * 0.9
library(ggplot2)
p2 <- ggplot(DA) + geom_point(aes(x=Index,y=Varaince),shape=1) + 
  geom_line(aes(x=Index,y=Varaince),color="gray50") + 
  theme_classic() + 
  labs(x="Index",y="Variance explained (%)") +
  theme(axis.title = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(vjust = 0.5),
        axis.text = element_text(size=12))

png("Var.png",res=400,height = 1500, width=3000)
p2
dev.off()

#Read the variety table and arrange
Var=read.delim("E:\\ubuntu\\SRR\\GT\\Var.table",header=TRUE)
Var <- Var[-1,]
"Cabbage" -> Var$ecotype_s[grep("Cabbage",Var$ecotype_s)]
"Kale" -> Var$ecotype_s[grep("Kale",Var$ecotype_s)]

#Plot PC1 X PC2
DATA=data.frame(svx.pca[,1:3],Var=Var$ecotype_s)
names(DATA)[1:3] = paste0("PC",1:3)
p <- ggplot(DATA) + geom_point(aes(x = PC1, y = PC2, color = Var),size = 2) + 
  labs(x = paste0("PC1 (",round(ll[1],2),"%)"),
       y = paste0("PC2 (",round(ll[2],2),"%)"),
       color = "Variety") +
  theme(axis.title = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(vjust = 0.5),
        axis.text = element_text(size=12),
        legend.text = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.line = element_line(colour = "black")
  )  + stat_ellipse(data=subset(DATA, DATA$Var %in% c("Broccoli","Cabbage","Kohlrabi","Cauliflower")),
                                aes(x = PC1, y = PC2, color = Var))

png("pca0.5_0.1_scale.png",res=400,height = 1900, width=3000)
p
dev.off()

