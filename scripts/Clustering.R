Clustering <- function(Data, resolution = 0.5, figure_format = "png"){
  if(!dir.exists(paste0(Data@tools$results_dir,"/Clustering/")))
    dir.create(paste0(Data@tools$results_dir,"/Clustering/"))
  results_dir <- paste0(Data@tools$results_dir,"/Clustering/",Data@tools$data_name)

  nsample <- length(unique(Data$orig.ident))
  
  # 判断是否进行批次效应
  if (nsample == 1) {
    cat("\n 1个样本 \n")
    # run cluster-----
    cat("\n %%%%% run FindNeighbors %%%%% \n")
    Data <- FindNeighbors(Data, dims = 1:50)
    cat("\n %%%%% run FindClusters %%%%% \n")
    Data <- FindClusters(Data,resolution = resolution)
    
    # plot umap-----
    plot_umap <- DimPlot(Data, reduction = "umap")
    # plot_tsne <- DimPlot(Data, reduction = "tsne")
    # plot_cluster <- plot_umap+plot_layout(guides = "collect")
    ggsave(paste0(results_dir,"_cluster.",figure_format),plot_umap,width = 3,height = 3,scale = 3)
  } else {
    # run umap before correct batch
    cat("\n %%%%% run UMAP %%%%% \n")
    Data <- RunUMAP(Data, dims = 1:50)
    plot_umap <- DimPlot(Data, reduction = "umap", group.by = "orig.ident")
    ggsave(
      filename = paste0(results_dir,"_before_correct_batch_effect.",figure_format),
      plot = plot_umap,
      width = 6, height = 6, scale = 3
    )
    
    cat("\n ", nsample ,"个样本 \n 矫正批次效应 \n")
    # run correct batch-----
    cat("\n %%%%% run Harmony %%%%% \n")
    Data <- suppressWarnings(harmony::RunHarmony(Data, "orig.ident"))
    
    # run UMAP based on harmony-----
    cat("\n %%%%% run UMAP %%%%% \n")
    Data <- RunUMAP(Data, dims = 1:50, reduction = "harmony")
    # TSNE
    # Data <- RunTSNE(Data, dims = 1:10, reduction = "harmony")
    cat("\n %%%%% run Neighbors %%%%% \n")
    Data <- FindNeighbors(Data, dims = 1:50, reduction = "harmony")
    cat("\n %%%%% run Clusters %%%%% \n")
    Data <- FindClusters(Data,resolution = resolution)
    
    # plot umap and tsne-----
    plot_samples <- DimPlot(Data, reduction = "umap", group.by = "orig.ident")
    ggsave(
      paste0(results_dir,"_after_correct_batch_effect.",figure_format),
      plot = plot_samples, 
      width = 6, height = 6, scale = 3
    )
  }
  
  # Data@assays$RNA@scale.data <- matrix()
  saveRDS(Data,file = paste0(results_dir, "_cluster.rds"))
  cat("\n----------Clustering finished----------\n")
  return(Data)
}