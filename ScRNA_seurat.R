########## Seurat analysis  
library(dplyr)
library(tibble)
library(cowplot)
library(Seurat)
library(ggplot2)
library(tidyverse)

dir.create('cluster1')
dir.create('cluster2')
dir.create('cluster3')
set.seed(123)  

dir = c('D:/06_scRNA/control/60filtered_feature_bc_matrix', 
        'D:/06_scRNA/control/61filtered_feature_bc_matrix',
        'D:/06_scRNA/control/62filtered_feature_bc_matrix',
        'D:/06_scRNA/kb/02KB',
        'D:/06_scRNA/kb/09KB',
        'D:/06_scRNA/kb/15KB')

names(dir) = c('NC1', 'NC2', 'NC3', 'KP1', 'KP2', 'KP3')
counts <- Read10X(data.dir = dir)

scRNA = CreateSeuratObject(counts)


###  Load datasets and create Seurat objects 
gf=gzfile('D:/06_scRNA/kb/02KB/features.tsv.gz','rt')
data<-read.table(gf)
ensemblQ=data$V2[!grepl('ENSSSCG', data$V2)]
counts_matrix = GetAssayData(scRNA, slot="counts")
scRNA_sub <- subset(scRNA, features = ensemblQ)
counts_matrixQ = GetAssayData(scRNA_sub, slot="counts")

table(grepl("^ENSSS",rownames(scRNA_sub)))
scRNA <- scRNA_sub

MT.genes <- c("ND1","ND2","ND3","ND4" ,"ND5","ND6","COX1","COX2","COX3","ATP6","ATP8","ND4L","CYTB")
MT_m <- match(MT.genes, rownames(scRNA@assays$RNA)) 
MT.genes <- rownames(scRNA@assays$RNA)[MT_m] 
MT.genes <- MT.genes[!is.na(MT.genes)] 
scRNA[["percent.mt"]]<-PercentageFeatureSet(scRNA, features=MT.genes) 
head(scRNA@meta.data)

HB.genes <- c("HBB","HGB")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
head(scRNA@meta.data)
col.num <- length(levels(as.factor(scRNA@meta.data$orig.ident)))
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)

### filter low quality cells
scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 30 )
scRNA <- NormalizeData(scRNA,normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst",nfeatures = 2000)
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
scRNA <- JackStraw(scRNA, num.replicate = 100)
scRNA <- ScoreJackStraw(scRNA, dims = 1:20)
JackStrawPlot(scRNA, dims = 1:20)
ElbowPlot(scRNA)
scRNA<- FindNeighbors(scRNA, dims = 1:20)
scRNA<- FindClusters(scRNA, resolution = 0.2)
scRNA <- RunUMAP(scRNA, dims = 1:50)

####Remove double cells
library(DoubletFinder)
library(tidyverse)
library(Seurat)
library(patchwork)

sweep.res.list <- paramSweep_v3(scRNA, PCs = c(1:20), sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
print(bcmvn)
pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
DoubleRate=ncol(scRNA)*8*1e-6
homotypic.prop <- modelHomotypic(scRNA$seurat_clusters)  
nExp_poi <- round(DoubleRate*ncol(scRNA))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
scRNA<- doubletFinder_v3(scRNA, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                         nExp = nExp_poi.adj, reuse.pANN = F, sct = T)


cols <-c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF" ,"#F7B6D2FF", "#E377C2FF", "#9EDAE5FF" ,"#BCBD22FF" ,"#17BECFFF",
         "#AEC7E8FF" ,"#FFBB78FF", "#98DF8AFF", "#FF9896FF" ,"#C5B0D5FF" ,"#C49C94FF")
      
DimPlot(scRNA, reduction = "umap", cols = cols, pt.size = 0.1,group.by='Cellstyle',label = T,label.box  = T,repel = T) 
















