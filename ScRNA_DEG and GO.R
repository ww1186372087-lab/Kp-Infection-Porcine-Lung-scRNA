##differentially expressed genes (DEG)

scRNA$Cellstyle.Group <- paste(scRNA$Cellstyle, scRNA$Group, sep = "_")
unique(scRNA$Cellstyle.Group)
scRNA$Celltype <- Idents(scRNA)
Idents(scRNA) <- "Cellstyle.Group"

##For example neutrophil cell
mydeg <- FindMarkers(scRNA_NEW,ident.1 = 'Neu_KP',ident.2 = 'Neu_NC', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)


if(!require(rvcheck))devtools::install_version("rvcheck",version="0.1.8",repos="http://cran.us.r-pro")
if(!require(clusterProfiler))BiocManager::install("clusterProfiler")
if(!require(org.Ss.eg.db))BiocManager::install("org.Ss.eg.db")
library(org.Ss.eg.db)
library(dplyr)
library(clusterProfiler)
library(ggplot2)

mygene <- Neu%>% top_n(100, wt = avg_log2FC) %>% rownames()

gene.df <- bitr(mygene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Ss.eg.db")


goBP <- enrichGO(gene.df$ENTREZID, org.Ss.eg.db, keyType = "ENTREZID", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)

data2plot<-goBP@result
data2plot <- data2plot[order(data2plot$qvalue,decreasing = F)[1:10],]
#算提供的基因占当前通路的比例
data2plot$BgRatio<-
  apply(data2plot,1,function(x){
    as.numeric(strsplit(x[3],'/')[[1]][1])
  })/apply(data2plot,1,function(x){
    as.numeric(strsplit(x[4],'/')[[1]][1])
  })

write.csv(goBP@result,'C:/Users/2023/Desktop/seven/Neu_BP.csv')




