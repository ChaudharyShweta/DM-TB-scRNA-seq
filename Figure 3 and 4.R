library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(harmony)
library(rliger)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(scCustomize)
library(viridis)
options(Seurat.object.assay.version = "v5")

d2 <- Read10X("F:/scRNA seq files/S2/sample_filtered_feature_bc_matrix/")
d3 <- Read10X("F:/scRNA seq files/S3/sample_filtered_feature_bc_matrix/")

srat_d2 <- CreateSeuratObject(d2$`Gene Expression`, project = "Healthy")
srat_d3 <- CreateSeuratObject(d3$`Gene Expression`, project = "DM")

srat_d2[["percent.mt"]] <- PercentageFeatureSet(srat_d2, pattern = "^mt-")
srat_d3[["percent.mt"]] <- PercentageFeatureSet(srat_d3, pattern = "^mt-")

summary(srat_d2@meta.data)
summary(srat_d3@meta.data)

srat_d2 <- subset(x=srat_d2, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & nCount_RNA > 1500 & nCount_RNA <20000 & percent.mt<5)
srat_d3 <- subset(x=srat_d3, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & nCount_RNA > 1500 & nCount_RNA <20000 & percent.mt<5)

data_list <- list()
data_list[["TB"]] <- srat_d2
data_list[["DM-TB"]] <- srat_d3


for (i in 1:length(data_list)) {
  data_list[[i]] <- NormalizeData(data_list[[i]], verbose = F)
  data_list[[i]] <- FindVariableFeatures(data_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = F)
}

data_anchors <- FindIntegrationAnchors(object.list = data_list, dims = 1:30)

# create list of common genes to keep
to_integrate <- lapply(data_list, row.names) %>% Reduce(intersect, .)

data_seurat <- IntegrateData(anchorset = data_anchors, dims = 1:30, features.to.integrate = to_integrate)

all.genes <- rownames(data_seurat)


#After Integration
DefaultAssay(data_seurat) <- "integrated"
data_seurat <- ScaleData(data_seurat, verbose = FALSE,features = all.genes)
data_seurat <- RunPCA(data_seurat)
DimPlot(data_seurat, reduction = "pca")
ElbowPlot(data_seurat)
data_seurat <- RunUMAP(data_seurat, reduction = "pca", dims = 1:30, verbose = F)
DimPlot(data_seurat, reduction = "umap", label = T)
DimPlot(data_seurat, reduction = "umap", split.by = "orig.ident") + NoLegend()

data_seurat <- FindNeighbors(data_seurat, dims = 1:30, verbose = T)
data_seurat <- FindClusters(data_seurat, verbose = T, resolution = 0.75)

# png(file="after integration.png", res=300, width=2000, height=1500)
DimPlot(data_seurat, label = T, label.size = 4)
# dev.off()
# 
# png(file="after integration split.png", res=300, width=3000, height=1500)
# DimPlot(data_seurat, label = F,split.by = "orig.ident") + NoLegend()
# dev.off()

clustercounts <- table(data_seurat@meta.data$seurat_clusters, data_seurat@meta.data$orig.ident)
clustercounts

# saveRDS(data_seurat, "integated_d2d3.rds") 

#### Lymphoid cell subclustering ####

lymphoid_cluster <- subset(data_seurat, idents =
                             c(0,1,2,3,4,5,6,7,9,13,15))
all.genes2 <- rownames(lymphoid_cluster)
lymphoid_cluster <- ScaleData(lymphoid_cluster, features = all.genes2)
lymphoid_cluster <- RunPCA(lymphoid_cluster, npcs = 30,features = VariableFeatures(object=lymphoid_cluster))
lymphoid_cluster <- RunUMAP(lymphoid_cluster, reduction = "pca", dims = 1:30)
lymphoid_cluster <- FindNeighbors(lymphoid_cluster, dims = 1:30, k.param = 50)
lymphoid_cluster <- FindClusters(lymphoid_cluster, resolution = c(0.72))


########################################################################################

lapply(c("dplyr","Seurat","HGNChelper"), library, character.only=TRUE)
library(ggforce)
#source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
#source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
source("F:/scRNA seq files/gene_sets_prepare.R")
source("F://scRNA seq files/sctype_score_.R")
# get cell-type-specific gene sets from our in-built database (DB)
#gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
gs_list = gene_sets_prepare("F:/scRNA seq files/all groups/13022024/Major_markers.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain


#Cell Type by Matrix
scRNAseqData = data_seurat
es.max = sctype_score(scRNAseqData = data_seurat[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
#Merge By Cluster
cL_resutls = do.call("rbind", lapply(unique(data_seurat@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(data_seurat@meta.data[data_seurat@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data_seurat@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

print(sctype_scores[,1:3], n = 25)

#Plotting in the graph
data_seurat@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  data_seurat@meta.data$customclassif[data_seurat@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

####### Figure 3B #######

png(file="major clusters.png", res=300, width=3000, height=2000)
DimPlot(data_seurat, reduction = "umap", label = F, repel = T,  label.size = 18, group.by = 'customclassif')+ theme(legend.text=element_text(size=30),text = element_text(size = 20),axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))
dev.off()

clustercounts <- table(data_seurat@meta.data$customclassif, data_seurat@meta.data$orig.ident)
clustercounts


#All cell types
gs_list = gene_sets_prepare("F:/scRNA seq files/all groups/13022024/all cell types.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain


scRNAseqData = data_seurat
es.max = sctype_score(scRNAseqData = data_seurat[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
#Merge By Cluster
cL_resutls = do.call("rbind", lapply(unique(data_seurat@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(data_seurat@meta.data[data_seurat@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data_seurat@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

#Low conf to unknown
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3], n = 25)

#Plotting in the graph
data_seurat@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  data_seurat@meta.data$customclassif[data_seurat@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

#saveRDS(data_seurat,  file = "integrated_d2d3")

####### Figure 3C merged #######

png(file="all cells clusters.png", res=300, width=2500, height=2000)
DimPlot(data_seurat, reduction = "umap", label = T, repel = T, label.size = 6, group.by = 'customclassif') + NoLegend() + theme(legend.text=element_text(size=30),text = element_text(size = 20),axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))
dev.off()

####### Figure 3C Split #######

png(file="all cells annotated split.png", res=300, width=3000, height=1500)
DimPlot(data_seurat, reduction = "umap", label = F, repel = T, group.by = 'customclassif',split.by = "orig.ident")+ theme(legend.text=element_text(size=12),text = element_text(size = 20),axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))
dev.off()

clustercounts <- table(data_seurat@meta.data$customclassif, data_seurat@meta.data$orig.ident)
clustercounts
write.csv(clustercounts, "Figure 3E.csv")

####### Figure 3E #######

source("F:/scRNA seq files/custom_seurat_functions_customclassif_H.R")
png(file="integrated cluster metrics.png", res=300, width=4000, height=2000)
plot_integrated_clusters(data_seurat) +  theme(
  panel.grid.major = element_blank(),  # Remove major grid lines
  panel.grid.minor = element_blank()   # Remove minor grid lines
)
dev.off()

####### Figure 3D #######
Idents(data_seurat) <- "customclassif"
Idents(data_seurat) <- factor(data_seurat@active.ident, sort(levels(data_seurat@active.ident)))
png(file="Dotplot all cell types.png", res=300, width=5500, height=1750)
DotPlot(data_seurat, features = c("Cd14","Itgax","Marco","Siglecf","Cd3d","Cd4","Mki67","Gzmk","Gzma","Cd8a","Nkg7","Klrd1","Ciita","Xcr1", "Cd24a", "Clec9a","Fcgr3","Treml4", "Cd36","Cd68", "Tlr2","Nos2","Ifi30","Flot1","Psmb9","Sell","Ccr7","Ncr1","Cd160","Klrk1" ,"Itgam","Siglech","Bst2","Maf","Bcl6","Icos","Tbx21","Ifng","Ccr5","Rorc","Il17a","Il23r","Foxp3","Ctla4","Stat5b"),cols = c("RdYlGn"),group.by = 'customclassif') + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + RotatedAxis() +theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA),
  text = element_text(size = 12),
  axis.text.x = element_text(size = 16), 
  axis.text.y = element_text(size = 18),
  panel.grid.major.x = element_line(color = "lightgrey"),
  panel.grid.major.y = element_line(color = "lightgrey") 
  
)
dev.off()

#### Cell Type annotation - lymphoid cluster ####

gs_list = gene_sets_prepare("H:/scRNA seq files/all groups/13022024/Finer_cell_types_n.xlsx", "Immune system")

es.max = sctype_score(scRNAseqData = lymphoid_cluster[["integrated"]]@scale.data, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
#Merge By Cluster
cL_resutls = do.call("rbind", lapply(unique(lymphoid_cluster@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(lymphoid_cluster@meta.data[lymphoid_cluster@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(lymphoid_cluster@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

#Low conf to unknown
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


##### Data pseudobulking and DESeq2 analysis#####

library(ExperimentHub)
library(DESeq2)
library(tidyverse)

mtdata <- lymphoid_cluster@meta.data

abdata_final = data.frame()
d2ab <- CreateSeuratObject(counts = d2$`Antibody Capture`)
abdata <- t(data.frame(d2ab@assays$RNA$counts))
fetchrep <- colnames(abdata)[max.col(abdata)]

row.names(abdata) <- gsub(x=row.names(abdata), pattern = "\\.", replacement = "-")
changereq <- grep("^[A-Z]", row.names(abdata))
row.names(abdata)[changereq] <- paste0(row.names(abdata)[changereq],"_1")
abdata <- data.frame(abdata)
abdata['Rep'] <- fetchrep
#Execute this command only once, Command it after
abdata_final = rbind(abdata_final,abdata)

#abdata_final2 = data.frame()
d3ab <- CreateSeuratObject(counts = d3$`Antibody Capture`)
abdata <- t(data.frame(d3ab@assays$RNA$counts))
fetchrep <- colnames(abdata)[max.col(abdata)]

row.names(abdata) <- gsub(x=row.names(abdata), pattern = "\\.", replacement = "-")
changereq <- grep("^[A-Z]", row.names(abdata))
row.names(abdata)[changereq] <- paste0(row.names(abdata)[changereq],"_2")
abdata <- data.frame(abdata)
abdata['Rep'] <- fetchrep
#Execute this command only once, Command it after
abdata_final = rbind(abdata_final,abdata)

abdata_for_merge <- select(abdata_final,Rep)


#mtdata['Rep'] <- ""
mtdatamerge = merge(mtdata, abdata_for_merge, by="row.names")
mtdatamerge1 = mtdatamerge[,-1]
rownames(mtdatamerge1) <- mtdatamerge[,1]
lymphoid_cluster@meta.data <- mtdatamerge1

#Pseudobulk analysis
lymphoid_cluster$samples <- paste0(lymphoid_cluster$orig.ident, lymphoid_cluster$Rep)

cts <- AggregateExpression(lymphoid_cluster,
                           group.by = c("customclassif","samples"),
                           assays = "RNA",
                           slot = "counts",
                           return.seurat = F)


#cts <- AggregateExpression(lymphoid_cluster, group.by = c("customclassif","samples"))$RNA
cts <- cts$RNA
cts <- data.frame(cts)
cts.t <- t(cts)
cts.t <- as.data.frame(cts.t)

SplitRows <- gsub('_.*', '', rownames(cts.t))

cts.split <- split.data.frame(cts.t,
                              f= factor(SplitRows))

cts.split.modi <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', "\\1", rownames(x))
  t(x)
})
#cts.split.modi$`Alveolar macrophages`
#gsub('.*_(.*)', "\\1", )

#Perform DESeq2
countdata <- cts.split.modi$Naïve.CD4..T.cells
colData <- data.frame(samples = colnames(countdata))
condi <- c("DM","DM","DM","DM","Healthy","Healthy","Healthy","Healthy")
colData['condition'] <- condi
colData <- column_to_rownames(colData, "samples")

dds <- DESeqDataSetFromMatrix(countdata, colData, design = ~condition)
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

dds$condition <- relevel(dds$condition, ref = "Healthy")

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name = "condition_DM_vs_Healthy")
res
write.csv(res, file="DE_results_DM_vs_H_Naive_CD4.csv")

######## Figure 4A #######

library(EnhancedVolcano)
png(file="Naive CD4+ T cells.png", res=300, width=3000, height=4000)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue', 
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                title = 'Naive CD4+ T-cells',
                drawConnectors = TRUE)
dev.off()



###### CellChat analysis #######

library(data.table)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

data.input = data_seurat[["RNA"]]$'data.DM' # normalized data matrix

#labels <- Idents(lymphoid_cluster)
labels <- data_seurat$customclassif
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
meta$condition <- ifelse(endsWith(rownames(meta), "_1"), "Healthy", "DM")
cell.use = rownames(meta)[meta$condition == "DM"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]

cellChat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.mouse

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB to use all CellChatDB for cell-cell communication analysis
# set the used database in the object
cellChat@DB <- CellChatDB.use

cellchat <- subsetData(cellChat) # This step is necessary even if using the whole database

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "cellchat_DM_all.rds")

#################################################################################################
data.input = data_seurat[["RNA"]]$'data.Healthy' # normalized data matrix

#labels <- Idents(lymphoid_cluster)
labels <- data_seurat$customclassif
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
meta$condition <- ifelse(endsWith(rownames(meta), "_1"), "Healthy", "DM")
cell.use = rownames(meta)[meta$condition == "Healthy"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]

cellChat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.mouse

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB to use all CellChatDB for cell-cell communication analysis
# set the used database in the object
cellChat@DB <- CellChatDB.use

cellchat <- subsetData(cellChat) # This step is necessary even if using the whole database

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "cellchat_Healthy_all.rds")

#CellChat comparative analysis
library(patchwork)
ptm = Sys.time()
updateCellChat

#Comparative cell chat analysis for all clusters
cellchat.Healthy <- readRDS("F:/scRNA seq files/all groups/29022023_h_dm/final figures/cellchat_Healthy_all.rds")
cellchat.DM <- readRDS("F:/scRNA seq files/all groups/29022023_h_dm/final figures/cellchat_DM_all.rds")

cellchat.Healthy <- netAnalysis_computeCentrality(cellchat.Healthy, slot.name = "netP") ## dataset 1
cellchat.DM <- netAnalysis_computeCentrality(cellchat.DM, slot.name = "netP") ## dataset2

object.list <- list(Healthy = cellchat.Healthy, DM= cellchat.DM)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use = c("#0072B2","#009E73"))+theme(text = element_text(size = 20))
gg1

png("Interactions Healthy vs DM all.png")
print(gg1)
dev.off()


######## Figure 4C #########

png("diffInteraction weight individual Healthy vs DM all.png",width=1200, height=600)
par(mfrow = c(1,2), xpd=TRUE)
png("diffInteraction weight Healthy vs DM all.png",width=650, height=650)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", vertex.label.cex = 1.5,vertex.size.max = 10)
dev.off()

######## Figure 4D #########

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, font.size = 11,color.use = c("#0072B2","#009E73")) 
png("information flow Healthy vs DM all stacked.png")
print(gg1)
dev.off()
write.csv(gg1$data, "Figure 4D.csv")

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Healthy"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, do.fast = F) 

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in Healthy
net.up <- subsetCommunication(cellchat, net = net, datasets = "Healthy",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "DM",ligand.logFC = -0.05, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

write.csv(net.down, "Downregulated interactions DM.csv")
write.csv(net.up, "Upregulated interactions DM.csv")

library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 7, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 7, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

png("outgoing signalling healthy.png",res=300,width=2300, height=2300)
print(ht1)
dev.off()

png("outgoing signalling DM.png",res=300,width=2300, height=2300)
print(ht2)
dev.off()

####### Figure 4E ########

pathways.show <- c("IFN-II")
# pathways.show <- c("FLT3")
# pathways.show <- c("GALECTIN")
# pathways.show <- c("IL16")
# pathways.show <- c("VISFATIN")
# pathways.show <- c("MIF")
# pathways.show <- c("LIGHT")
#pathways.show <- c("TGFb")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
png("IFN-II signaling pathway differential.png",res=300,width=3000, height=2300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

