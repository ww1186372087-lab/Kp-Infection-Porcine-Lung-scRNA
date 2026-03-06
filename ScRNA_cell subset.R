
###T cell subsets analyzed
T = scRNA[,scRNA@meta.data$seurat_clusters  %in% c(0,9,11)]
T <- NormalizeData(T,normalization.method = "LogNormalize", scale.factor = 10000)
T <- FindVariableFeatures(T, selection.method = "vst",nfeatures = 2000)
T <- ScaleData(T, features = rownames(T))
T <-RunPCA(T ,npcs = 20,verbose = FALSE)
T <- RunPCA(T, features = VariableFeatures(object = T))
T <- FindNeighbors(T , dims = 1:20)
T <- FindClusters(T , resolution = 0.1)
T <- RunUMAP(T , dims = 1:20)
T <- RunTSNE(T , dims = 1:20)
DimPlot(T, reduction = "umap",group.by = "Group")+
  DimPlot(T,reduction = "umap",label = TRUE)


s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
T<- CellCycleScoring(T, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
T@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()
DimPlot(T,reduction = "umap")

T<-ScaleData(T,vars.to.regress = c("S.Score","G2M.Score"))

T<-RunPCA(T,npcs = 20,verbose = FALSE)
T<- FindNeighbors(T, dims = 1:20)
T<- FindClusters(T, resolution = 0.2)
T<- RunUMAP(T, dims = 1:20)
DimPlot(T, reduction = "umap",group.by = "Phase")
DimPlot(T, reduction = "umap", label=TRUE)


mycolors <- c( '#E59CC4',  '#57C3F3','#D6E7A3','#FFBB78FF', '#E63863','#91D0BE','#F7F398')
my2colors <- c( "#2372A9","#CA2A28")
DimPlot(T, reduction = "umap", label = TRUE,cols = mycolors, pt.size = 0.1)+  
  DimPlot(T, reduction = "umap", label = TRUE,cols = my2colors, pt.size = 0.1,group.by='Group') 



###AM cell subsets analyzed

AM <- NormalizeData(AM,normalization.method = "LogNormalize", scale.factor = 10000)
AM <- FindVariableFeatures(AM, selection.method = "vst",nfeatures = 2000)
AM <- ScaleData(AM, features = rownames(AM))
AM<- RunPCA(AM, features = VariableFeatures(object = AM))
AM <- FindNeighbors(AM, dims = 1:20)
AM <- FindClusters(AM, resolution = 0.1)
AM <- RunUMAP(AM, dims = 1:20)
DimPlot(AM, reduction = "umap",label = T)+
  DimPlot(AM, reduction = "umap",group.by='Group')

mycolors <- c( '#4DBBD5FF', '#6778AE', '#E4C755','#F39B7FFF','#00A087FF','#91D1C2FF','#E59CC4','#9467BDFF')
my2colors <- c( "#2372A9","#CA2A28")

DimPlot(AM, reduction = "umap", label = TRUE,cols = mycolors) +
  DimPlot(AM, reduction = "umap", label = TRUE,cols = my2colors, group.by='Group') 
