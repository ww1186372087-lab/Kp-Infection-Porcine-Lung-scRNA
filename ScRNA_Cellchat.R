library(CellChat)
library(patchwork)
library(cowplot)
library(Seurat)
library(ggplot2)
library(DT)

data.input <-scRNA[["RNA"]]@data
meta <- scRNA@meta.data
unique(meta$Cellstyle)

cell.use <- rownames(meta)[meta$Cellstyle == "T"]
data.input <- data.input[, cell.use]
meta = meta[cell.use, ]
unique(meta$labels)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Cellstyle")
groupSize <- as.numeric(table(cellchat@idents))
groupSize

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat,features = NULL)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = T)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
class(df.net)

DT::datatable(df.net)
write.csv(df.net,"01.df.net.csv")

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize

par(mfrow = c(1,3), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 targets.use = 'T')

mat <- cellchat@net$weight
par(mfrow = c(3,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- df.net$pathway_name

par(mfrow = c(1,2), xpd=TRUE)
vertex.receiver = seq(1,4) 
netVisual_aggregate(cellchat, signaling =pathways.show[1],  
                    vertex.receiver = vertex.receiver,layout = 'hierarchy')

netVisual_aggregate(cellchat, signaling = pathways.show[1], layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show[1], layout = "chord")

length(pathways.show)

group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
names(group.cellType) <- levels(cellchat@idents)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_cell(cellchat, signaling = pathways.show[1],
                     group = group.cellType,
                     title.name = paste0(pathways.show[1], " signaling network"))

netVisual_chord_cell(cellchat, signaling = pathways.show[1],
                     title.name = paste0(pathways.show[1], " signaling network"))

p1 <- netAnalysis_contribution(cellchat, signaling = pathways.show[1],
                               title =  pathways.show[1])
p2 <- netAnalysis_contribution(cellchat, signaling = pathways.show)
cowplot::plot_grid(p1, p2, align = "h",ncol=2)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show[1],
                                 geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,]
vertex.receiver = seq(1,4)
p1<-netVisual_individual(cellchat, signaling = pathways.show,  
                         pairLR.use = LR.show, vertex.receiver = vertex.receiver,
                         layout = 'hierarchy')

vertex.receiver = seq(1,7)
p2<-netVisual_individual(cellchat, signaling = pathways.show,  
                         pairLR.use = LR.show, vertex.receiver = vertex.receiver,
                         layout = 'hierarchy')

netVisual_individual(cellchat, signaling = pathways.show[1], 
                     pairLR.use = LR.show, layout = "circle")

netVisual_individual(cellchat, signaling = pathways.show,
                     pairLR.use = LR.show, layout = "chord")

pathways.show.all <- cellchat@netP$pathways

levels(cellchat@idents)

vertex.receiver = seq(1,4)
dir.create('04.pathwat.show')
for (i in 1:length(pathways.show.all)) {
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0('04.pathwat.show/',pathways.show.all[i], "_L-R_contribution.pdf"),
         plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
gg

par(mfrow = c(1,3), xpd=TRUE)

netVisual_bubble(cellchat, sources.use = 4, 
                 targets.use = c(1:16), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = 4, 
                 targets.use = c("T","DC"), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = 4, 
                 targets.use = c(5:11), signaling = c("CCL","CXCL"), 
                 remove.isolate = FALSE)

cowplot::plot_grid(p1, p2,p3, align = "h",ncol=4)

par(mfrow = c(1,3), xpd=TRUE)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
pairLR.use

netVisual_bubble(cellchat, sources.use = c(4), targets.use = c(5:8), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE)

netVisual_bubble(cellchat, sources.use = 1, targets.use = c(5:11), 
                 signaling = c("CCL","CXCL"), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), 
                 signaling = c("CCL","CXCL"), remove.isolate = T)