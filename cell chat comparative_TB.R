library(CellChat)
library(patchwork)
ptm = Sys.time()
updateCellChat

#Comparative cell chat analysis for all clusters
cellchat.TB <- readRDS("F:/scRNA seq files/all groups/14022024_tb_tbdm/02052024/final figures/cellchat_TB_all.rds")
cellchat.DMTB <- readRDS("F:/scRNA seq files/all groups/14022024_tb_tbdm/02052024/final figures/cellchat_DMTB_all.rds")
cellchat.TB <- netAnalysis_computeCentrality(cellchat.TB, slot.name = "netP") ## dataset 1
cellchat.DMTB <- netAnalysis_computeCentrality(cellchat.DMTB, slot.name = "netP") ## dataset2
object.list <- list(TB = cellchat.TB, DMTB = cellchat.DMTB)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use = c("#CC79A7","#E69F00"))+theme(text = element_text(size = 20))
gg1
png("Interactions TB vs DMTB all.png")
print(gg1)
dev.off()

gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",size.text = 20)
gg2
png("Interactions weight TB vs DMTB all.png")
print(gg1)
dev.off()


png("diffInteraction weight TB vs DMTB all.png",width=650, height=650)
gg1<-netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", vertex.label.cex = 1.5,vertex.size.max = 10)
dev.off()

png("diffInteraction weight individual TB vs DMTB all.png")
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, font.size = 11,color.use = c("#CC79A7","#E69F00")) 
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE, font.size = 11,color.use = c("#CC79A7","#E69F00"))

gg1 + gg2

#Figure 5D
png("information flow TB vs DMTB all stacked.png")
print(gg1)
dev.off()
write.csv(gg1$data, "Figure 5D.csv")

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "TB"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, do.fast = F) 

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in TB
net.up <- subsetCommunication(cellchat, net = net, datasets = "TB",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "DMTB",ligand.logFC = -0.05, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

write.csv(net.down, "Downregulated interactions DM-TB.csv")
write.csv(net.up, "Upregulated interactions DM-TB.csv")


library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 7, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 7, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

png("outgoing signalling TB.png",res=300,width=2300, height=2300)
print(ht1)
dev.off()

png("outgoing signalling DMTB.png",res=300,width=2300, height=2300)
print(ht2)
dev.off()

#Figure 5F
pathways.show <- c("IFN-II")
# pathways.show <- c("CXCL")
# pathways.show <- c("GALECTIN")
# pathways.show <- c("IL16")
# pathways.show <- c("VISFATIN")
# pathways.show <- c("MIF")
# pathways.show <- c("CXCL")
# pathways.show <- c("TGFb")
# pathways.show <- c("PARs")
# pathways.show <- c("TNF")
# pathways.show <- c("COMPLEMENT")
# pathways.show <- c("CD40")
# pathways.show <- c("FASLG")
# pathways.show <- c("XCR")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
png("VISFATIN signaling pathway differential.png",res=300,width=3000, height=2300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

