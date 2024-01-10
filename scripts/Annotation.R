Annotation <- function(Data, ref_singler_dir = NULL, ref_cell_dex = NULL, ref_markers = NULL, 
                       re = FALSE, figure_format = "png"){
  results_dir <- paste0(Data@tools$results_dir,"/Annotation/")
  if(!dir.exists(results_dir)) dir.create(results_dir)
  results_dir <- paste0(results_dir,Data@tools$data_name)
  
  metadata_name <- NULL
  if (!is.null(ref_singler_dir)) {
    cat("\n %%%%% run SingleR %%%%% \n")
    ref <- ref_cell_dex
    
    annotation_dir <- paste0(results_dir,"_annotation_",ref,".rds")
    metadata_name <- append(
      metadata_name,
      values = paste0(c("cell_","cluster_"),ref)
    )
    if (file.exists(annotation_dir) & re == FALSE) {
      Data <- AddMetaData(Data, metadata = readRDS(file = annotation_dir))
    }else{
      # read ref
      ref_celldex <- readRDS(paste0(ref_singler_dir, ref, ".rds"))
      # cell annotation------
      cat("\n %%%%% run annotation by cell %%%%% \n")
      cell_annotation <- SingleR(test = Data@assays$RNA@data, ref = ref_celldex,
                                 assay.type.test = 1, labels = ref_celldex$label.main)
      Data <- AddMetaData(Data,metadata = cell_annotation$labels,col.name = metadata_name[1])
      # run cluster annotation------
      cat("\n %%%%% run annotation by cluster %%%%% \n")
      cluster_annotation <- SingleR(test = Data@assays$RNA@data, ref = ref_celldex, 
                                    clusters = Data$seurat_clusters, # 按类群注释
                                    assay.type.test = 1, labels = ref_celldex$label.fine)
      clusters <- Data$seurat_clusters
      levels(clusters) <- cluster_annotation$labels
      clusters <- factor(clusters)
      Data <- AddMetaData(Data,metadata = clusters,col.name = metadata_name[2])
      
      # save
      saveRDS(subset(Data@meta.data, select = metadata_name),
              file = annotation_dir)
    }
  }
  
  if (!is.null(ref_markers)) {
    print("根据已知markers注释------------------------------------------------")
    metadata_name <- append(
      metadata_name,
      values = "ref_markers_clusters"
    )
    
    CellTypeScoreMatrix <- matrix(
      nrow = length(levels(Data$seurat_clusters)),
      ncol = length(ref_markers),
      dimnames = list(levels(Data$seurat_clusters), names(ref_markers))
    )
    
    Data.all.markers <- readRDS(
      paste0(
        Data@tools$results_dir,
        "/DifferentialExpressing/",
        Data@tools$data_name,
        "_all_markers.rds"
      )
    )
    
    for (i in levels(Data$seurat_clusters)) {
      cluster_markers <- Data.all.markers[Data.all.markers$cluster == i,]
      cluster_markers <- arrange(cluster_markers, desc(abs(avg_log2FC)), p_val_adj)
      for (j in 1:length(ref_markers)) {
        scores <- sapply(ref_markers[[j]], function(x) {
          if (x %in% cluster_markers$gene) {
            which(cluster_markers$gene == x)/nrow(cluster_markers)
          } else {
            0
          }
        })
        CellTypeScoreMatrix[i,j] <- sum(scores)
      }
    }
    ClusterCellType <- apply(CellTypeScoreMatrix, 1, function(x) {names(which.max(x))})
    clusters <- Data$seurat_clusters
    levels(clusters) <- ClusterCellType
    clusters <- factor(clusters)
    Data <- AddMetaData(
      Data,
      metadata = clusters,
      col.name = "ref_markers_clusters"
    )
  }
  
  # 合并一级细胞类型相同的细胞
  celltypes <- levels(Data@meta.data[,metadata_name[2]])
  levels(Data@meta.data[,metadata_name[2]]) <- sapply(
    strsplit(celltypes, ":"), function(x) x[1]
  )
  n_celltypes <- length(levels(Data@meta.data[,metadata_name[2]]))
  
  
  # plot annotation umap--------------------------------------------------------
  cat("\n %%%%% run annotation by cells %%%%% \n")
  plot_umap <- DimPlot(Data,group.by = metadata_name[1], reduction = "umap",
                       split.by = "orig.ident")+
    theme(legend.position = "bottom")
  # plot_tsne <- DimPlot(Data,group.by = metadata_name[i], reduction = "tsne",
  #                      label = T , repel = T, label.size = 5)
  # plot_cluster <- plot_umap+plot_tsne+plot_layout(guides = "collect")
  ggsave(paste0(results_dir,"_",metadata_name[1],"_samples.",figure_format),
         plot_umap,width = 9*length(unique(Data$orig.ident)),height = 9)
  plot_umap_all <- DimPlot(Data,group.by = metadata_name[1], reduction = "umap")+
    theme(legend.position = "bottom")
  ggsave(paste0(results_dir,"_",metadata_name[1],".",figure_format),
         plot_umap_all, width = 9, height = 9)
  
  cat("\n %%%%% run annotation by clusters %%%%% \n")
  plot_umap <- DimPlot(Data,group.by = metadata_name[2], reduction = "umap",
                       split.by = "orig.ident",
                       label = T , repel = T, label.size = 5)+
    theme(legend.position = "bottom")
  # plot_tsne <- DimPlot(Data,group.by = metadata_name[i], reduction = "tsne",
  #                      label = T , repel = T, label.size = 5)
  # plot_cluster <- plot_umap+plot_tsne+plot_layout(guides = "collect")
  ggsave(paste0(results_dir,"_",metadata_name[2],"_samples.",figure_format),
         plot_umap,width = 9*length(unique(Data$orig.ident)),height = 9)
  plot_umap_all <- DimPlot(Data,group.by = metadata_name[2], reduction = "umap",
                       label = T , repel = T, label.size = 5)+
    theme(legend.position = "bottom")
  ggsave(paste0(results_dir,"_",metadata_name[2],".",figure_format),
         plot_umap_all, width = 9, height = 9)
  
  # 堆积图----------------------------------------------------------------------
  # 样本和细胞类型交叉表
  mydata <- table(Data$orig.ident,Data@meta.data[,metadata_name[2]])
  mydata <- mydata/rowSums(mydata)
  data <- as.data.frame(mydata)
  
  # 创建堆积图
  plot_bar <- ggplot(data=data, aes(Var1, Freq, fill=Var2)) +
    geom_bar(stat="identity", position="stack", color="black", width=0.7, linewidth=0.25) +
    labs(x = "", y = "Proportion %") +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(
      panel.background=element_rect(fill="white", colour="black", linewidth = 0.25),
      axis.line=element_line(colour="black", linewidth = 0.25),
      axis.title=element_text(size=13, color="black"),
      axis.text = element_text(size=12, color="black"),
      axis.text.x = element_text(angle = 30,vjust = 0.5),
      # legend.position= "bottom",
      legend.title = element_blank()
    )
  if (n_celltypes < 13) {
    plot_bar <- plot_bar+scale_fill_manual(
      values = RColorBrewer::brewer.pal(n_celltypes, "Set3")
    )
  }else{
    plot_bar <- plot_bar+scale_fill_manual(
      values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_celltypes)
    )
  }
  ggsave(paste0(results_dir,"_barplot.",figure_format),
         plot_bar,width = 3+1*length(unique(Data$orig.ident)),height = 5)
  
  
  # 条形图----------------------------------------------------------------------
  # metadata_name <- "cell_MonacoImmuneData"
  # ncluster <- as.data.frame(table(Data[[metadata_name]]))
  # names(ncluster)[1] <- "cell_type"
  # plot_bar <- ggplot(ncluster,aes(x = cell_type,y = Freq,fill = cell_type)) +
  #   geom_bar(stat = 'identity') +
  #   geom_text(aes(label = Freq),hjust = 1) +
  #   coord_flip()+
  #   labs(x=NULL,y=NULL,title = metadata_name)+
  #   NoLegend()
  # ggsave(paste0(results_dir,"_bar",".",figure_format),
  #        plot_bar,width = 10,height = 5)
  
    # barplot for comparison------
  # if (length(unique(Data$orig.ident)) > 1) {
  #   table_ba <- table(Data$orig.ident,Data$cell_annotation)
  #   percentage_cell <- table_ba/as.vector(table(Data$orig.ident))
  #   ab <- c(rep(rownames(table_ba)[1],ncol(table_ba)),rep(rownames(table_ba)[2],ncol(table_ba)))
  #   composition <- data.frame(cell_type=colnames(table_ba),p=as.vector(t(percentage_cell)),
  #                             ab=ab)
  #   composition$cell_type <- factor(composition$cell_type,levels = colnames(table_ba),ordered = T)
  #   composition$ab <- factor(composition$ab,levels = rownames(table_ba),ordered = T)
  #   p <- ggplot(composition)+
  #     geom_bar(aes(x = cell_type,y = p,fill = ab),
  #              stat = 'identity',position = 'dodge')+
  #     scale_y_continuous(limits = c(0,1))+
  #     labs(x = '细胞类型',y = '比例',fill='前后')+
  #     theme_set(theme_gray(base_size = 18))
  #   ggsave(filename = paste0(results_dir,"_composition_MonacoImmune.",figure_format)
  #          ,plot = p,width = 15,height = 8)
  # }
  
  # saveRDS(Data,file = paste0(results_dir, "_annotation.rds"))
  cat("\n --------------------Annotation finished-------------------- \n")
  return(Data)
}
