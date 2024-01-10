DimensionalityReduction <- function(Data, pc_num = 50, figure_format = "png"){
  results_dir <- paste0(Data@tools$results_dir,"/DimensionalityReduction/",Data@tools$data_name)
  
  # run Normalize-----
  cat("\n %%%%% run Normalize %%%%% \n")
  Data <- NormalizeData(Data)
  
  # run find VG-----
  cat("\n %%%%% run FindVariableFeatures %%%%% \n")
  Data <- FindVariableFeatures(Data)
  
  top10 <- head(VariableFeatures(Data), 10)
  # Demonstration of highly variable genes
  plot_hvg <- VariableFeaturePlot(Data) + theme(legend.position = "bottom")
  plot_hvg <- LabelPoints(
    plot = plot_hvg, points = top10, repel = TRUE,xnudge = 0,ynudge = 0
  )
  suppressWarnings(ggsave(
    paste0(results_dir,"_hvg.",figure_format),
    plot_hvg, width = 5, height = 5, bg = "white"
  ))
  
  # run Scale-----
  cat("\n %%%%% run Scale %%%%% \n")
  Data <- ScaleData(Data, features = rownames(Data))
  
  
  # run PCA----
  cat("\n %%%%% run PCA %%%%% \n")
  Data <- suppressMessages(
    RunPCA(Data, features = VariableFeatures(object = Data), npcs = pc_num)
  )
  
  # plot elbow----
  plot_elbow <- ElbowPlot(Data)
  ggsave(paste0(results_dir,"_elbow.",figure_format),plot_elbow,width = 6,height = 3,bg = "white")
  
  # plot dimloading----
  plot_dimloading <- VizDimLoadings(Data, dims = 1:16, reduction = "pca")
  ggsave(paste0(results_dir,"_dimloading.",figure_format),plot_dimloading,width = 20,height = 20)
  
  # plot pc1 - pc2----
  plot_pca <- DimPlot(Data, reduction = "pca",group.by = "orig.ident")+NoLegend()
  ggsave(paste0(results_dir,"_pca.",figure_format),plot_pca,width = 6,height = 5)
  
  # plot 热图----
  png( 
    filename = paste0(results_dir,"_dimheatmap.",figure_format),
    width = 1920,
    height = 1080,
    units = "px",
    bg = "white",
    res = 100)
  DimHeatmap(Data, dims = 1:9, cells = 500)
  dev.off()
  
  # Data@assays$RNA@scale.data <- matrix()
  # saveRDS(Data,file = paste0(results_dir, "_dr.rds"))
  cat("\n ----------Dimensionality reduction finished---------- \n")
  return(Data)
}