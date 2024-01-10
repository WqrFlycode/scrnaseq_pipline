QualityControl <- function(Data, figure_format = "png"){
  results_dir <- paste0(
    Data@tools$results_dir,
    "/QualityControl/",
    Data@tools$data_name
  )
  
  
  number_before <- dim(Data) # 细胞和基因的数量
  qc_index <- c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")
  # compute the proportion of mito-genes expression
  Data[["percent.mito"]] <- PercentageFeatureSet(Data, pattern = "^MT-") 
  # compute the proportion of ribo-genes expression
  Data[["percent.ribo"]] <- PercentageFeatureSet(Data, pattern = "^RP[SL]")
  # results is percentage times 100
  # plot violion before QC-----
  plot_qc <- suppressWarnings(
    VlnPlot(
      Data, 
      features = qc_index, 
      ncol = length(qc_index),
      group.by = "orig.ident"
    )
  )
  ggsave(
    paste0(results_dir,"_before_qc.",figure_format),
    plot_qc, width = 12, height = 6
  )
  
  # QC-----
  cat("\n %%%%% run quality control %%%%% \n")
  Data <- subset(
    Data,
    subset = nFeature_RNA > 100 & 
      nCount_RNA > 100 & # nCount_RNA < 1e5 &
      percent.mito < 30 &
      percent.ribo < 50,
    features = row.names(Data)[rowSums(GetAssayData(Data,slot = "counts")>0)>0]
  )
  meta.data <- Data@meta.data
  write.csv(meta.data, paste0(results_dir,"_qc_metadata.csv"))
  
  # output-----
  number_after <- dim(Data)
  difference <- number_before - number_after
  cat("\n 剔除细胞",difference[1],"个，基因",difference[2],"个 \n")
  number <- matrix(c(number_before,number_after,difference),ncol = 2,byrow = T)
  rownames(number) <- c("before","after","discard")
  colnames(number) <- c("features","cells")
  write.csv(number,file = paste0(results_dir,"_difference.csv"))
  # plot violion after QC-----
  plot_qc <- suppressWarnings(
    VlnPlot(
      Data,
      features = c(qc_index),
      ncol = length(qc_index),
      group.by = "orig.ident"
    )
  )
  ggsave(
    paste0(results_dir,"_before_qc.",figure_format),
    plot_qc,width = 12,height = 6
  )
  # QC-index scatter-----
  suppressWarnings(
    plot_FeatureScatter <- 
      FeatureScatter(
        Data,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",
        group.by = "orig.ident"
      )+NoLegend()+
      FeatureScatter(
        Data,feature1 = "nCount_RNA",feature2 = "percent.mito",
        group.by = "orig.ident"
      )+NoLegend()+
      FeatureScatter(
        Data,feature1 = "nCount_RNA",feature2 = "percent.ribo",
        group.by = "orig.ident"
      )+NoLegend()
  )
  
  ggsave(paste0(results_dir,"_qc_index_scatter.",figure_format),
         plot_FeatureScatter,width = 12,height = 6)
  
  
  # saveRDS(Data,file = paste0(results_dir, "_qc.rds"))
  cat("\n----------Quality control finished----------\n")
  return(Data)
}
