rm(list = ls())
# Some functions to be used
auc_thresh_kmeans <- function(regulonAUC){
  
  ## Iterate over each regulon in the AUC matrix
  regulons <- rownames(regulonAUC)
  
  kmeans_thresholds <- list()
  
  print("Processing regulon distributions...")
  
  for(regulon_no in 1:length(regulons)){
    
    # print(regulon_no)
    
    regulon <- regulons[regulon_no]
    
    df <- data.frame("auc" = as.numeric(regulonAUC[regulon,]),
                     "cells" = as.character(colnames(regulonAUC)),
                     "regulon" = regulon)
    
    ## Remove cells with 0 AUC as they interfere with kmeans clustering
    df <- df %>%
      subset(auc > 0)
    
    if(dim(df)[1]<=2){
      kmeans_thresholds[[regulon]] <- median(df$auc)
      next
    }
    
    km <- kmeans(df$auc,centers=2)
    df$cluster <- as.factor(km$cluster)
    # print("kmean_run")
    
    cluster1_max <- max(subset(df,cluster == 1)$auc)
    cluster2_max <- max(subset(df,cluster == 2)$auc)
    
    if(cluster1_max > cluster2_max){
      
      df <- df %>%
        mutate("cluster" = gsub(2,3,cluster)) %>%
        mutate("cluster" = gsub(1,2,cluster)) %>%
        mutate("cluster" = gsub(3,1,cluster))
    }
    
    df <- df %>%
      arrange(desc(auc))
    
    df_sub <- df %>%
      subset(cluster == 1)
    
    auc_thresholds <- df_sub[1,]$auc
    
    kmeans_thresholds[[regulon]] <- auc_thresholds
    
  }
  
  
  print("Done evaluating thresholds...")
  return(kmeans_thresholds)
  
}
binarize_regulons <- function(regulonAUC,
                              thresholds){
  
  binary_regulon_list <- list()
  binary_regulon_df <- data.frame("cells" = c())
  
  for(regulon_no in 1:length(names(thresholds))){
    
    # print(regulon_no)
    
    regulon <- names(thresholds)[regulon_no]
    
    auc_df <-  data.frame("auc" = as.numeric(regulonAUC[regulon,]),
                          "cells"= as.character(names(regulonAUC[regulon,])))
    # colnames(auc_df) = c("cells","auc")
    auc_df
    auc_df <- auc_df %>%
      mutate("regulon"= if_else(auc >= thresholds[regulon],1,0)) %>%
      select(-auc)
    
    colnames(auc_df) <- c("cells",regulon)
    
    binary_regulon_list[[regulon]] <- auc_df
  }
  joined_bin_reg <- binary_regulon_list %>% 
    purrr::reduce(left_join,by="cells") %>% 
    tibble::column_to_rownames("cells")
  return(as.matrix(t(joined_bin_reg)))
  
  # return(binary_regulon_list)
}

# read_data and set flags (manually)
seurat_data <- readRDS("integrated_d1d4.rds")
auc_1 = read.csv("auc_mtx.srat_d1.csv", header = FALSE)
auc_2 = read.csv("auc_mtx.srat_d4.csv", header = FALSE)
# set to false if looms were created using the non-integrated seurat object
IS_COUNT_DATA_EXTRACTED_FROM_INTEGRATED_RDS = TRUE
NUMBER_OF_GENES_TO_PLOT = 25


##remaining is automated
# save counts to loom

# SCopeLoomR::build_loom(file.name = "srat_d1.loom", dgem = seurat_data@assays$RNA$counts.TB)
# SCopeLoomR::build_loom(file.name = "srat_d4.loom", dgem = seurat_data@assays$RNA$`counts.DM-TB`)


#set data
seurat_meta <- seurat_data@meta.data
table(seurat_meta$orig.ident)
labels.table = table(sapply(rownames(seurat_meta), function(i) unlist(strsplit(i, split = "_")[[1]][2])))

rownames_1 = as.character(auc_1[, 1])
rownames_2 = as.character(auc_2[, 1])
rownames_1[2]
rownames_2[2]

# run only if rownames do not contain _1 and _2 at the end of cell name 
if(!IS_COUNT_DATA_EXTRACTED_FROM_INTEGRATED_RDS){
  
if(length(rownames_1)-1==labels.table[1]){
  rownames_1 = sapply(rownames_1, function(i){paste(i,names(labels.table)[1], sep = "_")})
  rownames_2 = sapply(rownames_2, function(i){paste(i,names(labels.table)[2], sep = "_")})
} else if(length(rownames_2)-1==labels.table[1]){
  rownames_2 = sapply(rownames_2, function(i){paste(i,names(labels.table)[1], sep = "_")})
  rownames_1 = sapply(rownames_1, function(i){paste(i,names(labels.table)[2], sep = "_")})
}
}

## rownames here should have _1 and _2 at their ends
rownames_1[2]
rownames_2[2]


#Step 1: Load and manage data; read(d) => t(d)

colnames(auc_1) = auc_1[1, ]
colnames(auc_2) = auc_2[1, ]
auc_1 = auc_1[-1, ]
auc_2 = auc_2[-1, ]
auc_1 = auc_1[, -1]
auc_2 = auc_2[, -1]
auc_1 = as.data.frame(apply(auc_1, 2, as.numeric))
auc_2 = as.data.frame(apply(auc_2, 2, as.numeric))
row.names(auc_1) = rownames_1[-1]
rownames(auc_2) = rownames_2[-1]

auc_1_t = t(auc_1)
auc_2_t = t(auc_2)

library(tidyverse)
thresholds_1 = auc_thresh_kmeans(regulonAUC = auc_1_t)
auc_1_t.binary = binarize_regulons(regulonAUC = auc_1_t, thresholds = thresholds_1)

thresholds_2 = auc_thresh_kmeans(regulonAUC = auc_2_t)
auc_2_t.binary = binarize_regulons(regulonAUC = auc_2_t, thresholds = thresholds_2)



# Step 2: Perform mann-whitney test on each common regulon taking data from both groups
# rownames(auc_1_t) %in% rownames(auc_2_t)
common_regulons = intersect(rownames(auc_1_t), rownames(auc_2_t))

regulon_pvalue = c()
for (regulon in common_regulons) {
  wilcox.result = wilcox.test(as.numeric(auc_1_t[regulon, ]), as.numeric(auc_2_t[regulon, ]))
  regulon_pvalue = rbind(regulon_pvalue, c(regulon, wilcox.result$p.value))
}
colnames(regulon_pvalue) = c("Regulon", "p-value")
regulon_pvalue = as.data.frame(regulon_pvalue)
regulon_pvalue$`p-value` = as.numeric(regulon_pvalue$`p-value`)

regulon_binary_pvalue = c()
for (regulon in common_regulons) {
  # print(regulon)
  # wilcox.result = wilcox.test(auc_1_t.binary[regulon, ], auc_2_t.binary[regulon, ])
  input_matrix = as.table(rbind(table(auc_1_t.binary[regulon, ]), table(auc_2_t.binary[regulon, ])))
  dimnames(input_matrix) = list(sample=c("D2","D3"), binary_auc=c("0","1"))
  chisq.result = chisq.test(input_matrix)
  regulon_binary_pvalue = rbind(regulon_binary_pvalue, c(regulon, chisq.result$p.value))
}

colnames(regulon_binary_pvalue) = c("Regulon", "p-value")
regulon_binary_pvalue = as.data.frame(regulon_binary_pvalue)
regulon_binary_pvalue$`p-value` = as.numeric(regulon_binary_pvalue$`p-value`)


# Step 3: filter to all combinations less than 0.05 p-value => Significant differential regulons

signif.different.regulons = regulon_pvalue[which(regulon_pvalue$`p-value` <= 0.001 &
                                                           regulon_pvalue$`p-value` > 0),]
signif.different.regulons = signif.different.regulons[order(signif.different.regulons$`p-value`, decreasing = FALSE),][1:NUMBER_OF_GENES_TO_PLOT,1]

auc_1_t.signif = auc_1_t[signif.different.regulons, ]
auc_2_t.signif = auc_2_t[signif.different.regulons, ]


combined_auc.diff = cbind(auc_1_t.signif, auc_2_t.signif)

signif.different.binary_regulons = regulon_binary_pvalue[which(regulon_binary_pvalue$`p-value` <= 0.001 &
                                                           regulon_binary_pvalue$`p-value` > 0),]
signif.different.binary_regulons = signif.different.binary_regulons[order(signif.different.binary_regulons$`p-value`, decreasing = FALSE),][1:NUMBER_OF_GENES_TO_PLOT,1]
auc_1_t.binary_signif = auc_1_t.binary[signif.different.binary_regulons, ]
auc_2_t.binary_signif = auc_2_t.binary[signif.different.binary_regulons, ]


combined_auc.binary_diff = cbind(auc_1_t.binary_signif, auc_2_t.binary_signif)
# combined_auc.diff.scaled = t(apply(combined_auc.diff, 1, scale))
# combined_auc.diff.binary = combined_auc.diff > 0.1
# Step 4: plot all data as a clustermap after combining and grouping

seurat_data_final = data.frame(cbind(celltype=seurat_meta$customclassif, 
                                     seurat_cluster=seurat_meta$seurat_clusters, 
                                     identifier=seurat_meta$orig.ident))


combined_auc.diff.t = t(combined_auc.diff)
combined_auc.binary_diff.t = t(combined_auc.binary_diff)
#associate celltypes to the cells
plot_data_auc = merge(x=seurat_meta, y=combined_auc.diff.t, all=TRUE, by="row.names")
plot_data_auc.binary = merge(x=seurat_meta, y=combined_auc.binary_diff.t, all=TRUE, by="row.names")

# install.packages("pheatmap")
# install.packages("magick")
library(RColorBrewer)
library(ComplexHeatmap)

ident_1 = which(plot_data_auc$orig.ident==unique(plot_data_auc$orig.ident)[1])
ident_2 = which(plot_data_auc$orig.ident==unique(plot_data_auc$orig.ident)[2])
pdf(paste0(plot_data_auc$orig.ident[ident_1[1]],"_",plot_data_auc$orig.ident[ident_2[1]],".pdf"))

dend1 = cluster_within_group(t(plot_data_auc[ident_1,9:dim(plot_data_auc)[2]]), plot_data_auc$customclassif[ident_1])
dend2 = cluster_within_group(t(plot_data_auc[ident_2,9:dim(plot_data_auc)[2]]), plot_data_auc$customclassif[ident_2])
dend = merge(dend1, dend2)

# Create a palette with 16 colors
custom_colors_cellname <- colorRampPalette(brewer.pal(8, "Accent"))(16)

# Ensure the custom palette is named appropriately
custom_colors_cellname <- setNames(custom_colors_cellname[seq_along(names(table(plot_data_auc.binary$customclassif)))], 
                                   nm = names(table(plot_data_auc.binary$customclassif)))

# Modify the HeatmapAnnotation to use the custom colors for cellname
ha_1 <- HeatmapAnnotation(
  identifier = plot_data_auc$orig.ident[ident_1], 
  cellname = plot_data_auc$customclassif[ident_1], show_annotation_name = FALSE,
  col = list(
    identifier = setNames(c("#E69F00", "#CC79A7"), nm = unique(plot_data_auc.binary$orig.ident)),
    cellname = custom_colors_cellname
  )
)

hm_1 = Heatmap(
  t(plot_data_auc[ident_1,9:dim(plot_data_auc)[2]]),
  name = "AUC",
  # column_split = paste(plot_data_auc$orig.ident, plot_data_auc$customclassif),
  cluster_columns = dend1,
  # cluster_column_slices = FALSE,
  use_raster = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE, row_title = "Regulons", column_title = plot_data_auc$orig.ident[ident_1[1]],
  top_annotation = ha_1
)


ha_2 = HeatmapAnnotation(
  identifier = plot_data_auc$orig.ident[ident_2], 
  cellname = plot_data_auc$customclassif[ident_2], show_annotation_name = TRUE,
  col = list(identifier = setNames(c("#E69F00","#CC79A7"), nm = unique(plot_data_auc.binary$orig.ident)),
             cellname = custom_colors_cellname
  )
)

hm_2 = Heatmap(
  t(plot_data_auc[ident_2,9:dim(plot_data_auc)[2]]),
  name = "AUC",
  # column_split = paste(plot_data_auc$orig.ident, plot_data_auc$customclassif),
  cluster_columns = dend2,
  # cluster_column_slices = FALSE,
  use_raster = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE, row_title = "Regulons", column_title = plot_data_auc$orig.ident[ident_2[1]],
  top_annotation = ha_2
)


# top_hm = Heatmap(rbind(celltype = plot_data_auc$customclassif), 
#                  column_split = plot_data_auc$customclassif, 
#                  name = "Cell Types", 
#                  use_raster = FALSE,
#                  show_row_dend = TRUE)

draw(hm_1 + hm_2)
# main_hm
dev.off()

hm_1+hm_2

