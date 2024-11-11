plot_integrated_clusters = function (srat) { 
  ## take an integrated Seurat object, plot distributions over orig.ident
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  
  
  count_table <- table(srat@meta.data$customclassif, srat@meta.data$orig.ident)
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)

  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
  
  #sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  sorted_labels <- paste(sort(as.character(levels(cluster_size$cluster)),decreasing = T))
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
  
  colnames(melt_mtx)[2] <- "dataset"
  
  
  p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
    theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("")
  p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
    geom_bar(position="dodge", stat="identity") + theme_bw() + coord_flip() + 
    scale_fill_manual(values=c("#009E73", 
                               "#0072B2")) +
    ylab("Cells in each dataset") + xlab("Cluster") + theme(legend.position="top")+ theme(legend.text=element_text(size=12),text = element_text(size = 14),axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 12))
  
  p2 + theme(plot.margin = unit(c(1, 3, 1, 1), "cm")) 
  
  }