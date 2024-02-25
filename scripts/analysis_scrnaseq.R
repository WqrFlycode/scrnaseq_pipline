analysis_scrnaseq <- function(Data,
                              min_nFeature = 100, min_nCount = 100, 
                              max_percent_mito = 30, max_percent_ribo = 50,max_percent_hb = 30,
                              pc_num = 50, resolution = 0.8,
                              ref_singler_dir = NULL, ref_cell_dex = NULL, ref_markers = NULL){
  info <- Data@tools$info
  results_path <- paste0(info$results_dir, info$data_name)
  sink(file = paste0(results_path, "_analysis_log.txt"),split = TRUE)
  
  info$raw_dim <- dim(Data)
  names(info$raw_dim) <- c("gene", "cell")
  info$qc_index <- c(
    min_nFeature, min_nCount, max_percent_mito, max_percent_ribo
  )
  names(info$qc_index) <- c(
    "min_nFeature", "min_nCount", "max_percent_mito", "max_percent_ribo"
  )
  info$pc_num <- pc_num
  info$cluster_resolution <- resolution
  info$ref_singler_dir <- ref_singler_dir
  info$ref_cell_dex <- ref_cell_dex
  info$ref_markers <- ref_markers
  
  # QC
  # qc_index <- c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")
  # compute the proportion of mito-genes expression
  if(! "percent.mito" %in% names(Data@meta.data)){
    Data[["percent.mito"]] <- PercentageFeatureSet(Data, pattern = "^MT-") 
  }
  # compute the proportion of ribo-genes expression
  if(! "percent.ribo" %in% names(Data@meta.data)){
    Data[["percent.ribo"]] <- PercentageFeatureSet(Data, pattern = "^RP[SL]")
  }
  # 计算红血细胞基因比例
  if(! "percent.hb" %in% names(Data@meta.data)){
    Data[["percent.hb"]] <- PercentageFeatureSet(Data, pattern = "^Hb[^(p)]")
  }
  # percentage times 100
  
  cat("\n %%%%% run quality control %%%%% \n")
  Data <- subset(
    Data,
    subset = nFeature_RNA > min_nFeature & 
      nCount_RNA > min_nCount & # nCount_RNA < 1e5 &
      percent.mito < max_percent_mito &
      percent.ribo < max_percent_ribo &
      percent.hb < max_percent_hb,
    features = row.names(Data)[rowSums(GetAssayData(Data,slot = "counts")>0)>0]
  )
  info$filter_dim <- dim(Data)
  names(info$filter_dim) <- c("gene", "cell")
  
  # run Normalize-----
  cat("\n %%%%% run Normalize %%%%% \n")
  Data <- NormalizeData(Data)
  
  # run find VG-----
  cat("\n %%%%% run FindVariableFeatures %%%%% \n")
  Data <- FindVariableFeatures(Data)
  
  # run Scale-----
  cat("\n %%%%% run Scale %%%%% \n")
  Data <- ScaleData(Data)
  
  # run PCA----
  cat("\n %%%%% run PCA %%%%% \n")
  Data <- suppressMessages(
    RunPCA(Data, features = VariableFeatures(object = Data), npcs = pc_num)
  )
  
  nsample <- length(unique(Data$orig.ident))
  # 判断是否进行批次效应
  if (nsample == 1) {
    cat("\n 1个样本 \n")
    # run cluster-----
    cat("\n %%%%% run FindNeighbors %%%%% \n")
    Data <- FindNeighbors(Data, dims = 1:pc_num)
    cat("\n %%%%% run FindClusters %%%%% \n")
    Data <- FindClusters(Data, resolution = resolution)
    cat("\n %%%%% run UMAP %%%%% \n")
    Data <- RunUMAP(Data, dims = 1:pc_num)
    cat("\n %%%%% run TSNE %%%%% \n")
    Data <- RunTSNE(Data)
  } else {
    # run correct batch-----
    cat("\n ", nsample ,"个样本 \n 矫正批次效应 \n")
    cat("\n %%%%% run Harmony %%%%% \n")
    Data <- suppressWarnings(harmony::RunHarmony(Data, "orig.ident"))
    
    # run UMAP based on harmony-----
    cat("\n %%%%% run UMAP %%%%% \n")
    Data <- RunUMAP(Data, dims = 1:pc_num, reduction = "harmony")
    cat("\n %%%%% run TSNE %%%%% \n")
    Data <- RunTSNE(Data, dims = 1:pc_num, reduction = "harmony")
    cat("\n %%%%% run Neighbors %%%%% \n")
    Data <- FindNeighbors(Data, dims = 1:pc_num, reduction = "harmony")
    cat("\n %%%%% run Clusters %%%%% \n")
    Data <- FindClusters(Data, resolution = resolution)
  }
  
  # run annotation-----
  cat("\n %%%%% run annotation %%%%% \n")
  Data <- Annotation(
    Data = Data,
    ref_singler_dir = ref_singler_dir, 
    ref_cell_dex = ref_cell_dex, 
    ref_markers = ref_markers, 
    figure_format = "png"
  )
  # levels(Data$singler_by_cluster) <- sapply(
  #   strsplit(levels(Data$singler_by_cluster), ":"), 
  #   function(x) x[1]
  # )
  info$celltypes <- list(
    by_seurat = levels(Data$seurat_clusters),
    by_cell = unique(Data$singler_by_cell),
    by_cluster = levels(Data$singler_by_cluster)
  )
  
  # run find marker-----
  all.markers <- list()
  cat("\n %%%%% run FindAllMarkers by seurat cluster %%%%% \n")
  Data.all.markers.by.cluster <- FindAllMarkers(
    Data, 
    features = Data@assays$RNA@var.features
  )
  all.markers$by.cluster <- Data.all.markers.by.cluster
  
  Idents(Data) <- Data@meta.data[,"singler_by_cluster"]
  cat("\n %%%%% run FindAllMarkers by cell types %%%%% \n")
  Data.all.markers.by.celltype <- FindAllMarkers(
    Data, 
    features = Data@assays$RNA@var.features
  )
  all.markers$by.celltype <- Data.all.markers.by.celltype
  
  all.markers_dir <- paste0(results_path, "_all.markers.rds")
  info$all.markers_dir <- all.markers_dir
  saveRDS(all.markers, all.markers_dir)
  Data@tools$all.markers <- all.markers
  
  Idents(Data) <- Data@meta.data$seurat_clusters
  
  info$by_cluster <- "seurat_clusters"
  # run TrajectoryAnalysis-----
  tryCatch({
    cds <- TrajectoryAnalysis(
      Data = Data, by_cluster = info$by_cluster
    )
    cds_dir <- paste0(results_path, "_cds.rds")
    info$cds_dir <- cds_dir
    saveRDS(cds, cds_dir)
    Data@tools$cds <- cds
  }, error = function(e){
    message("%%% TrajectoryAnalysis error %%%")
  })
  
  
  # CellCommunication-----
  tryCatch({
    cellchat <- CellCommunication(
      Data = Data, by_cluster = "singler_by_cluster"
    )
    cellchat_dir <- paste0(results_path, "_cellchat.rds")
    info$cellchat_dir <- cellchat_dir
    saveRDS(cellchat, cellchat_dir)
    Data@tools$cellchat <- cellchat
  }, error = function(e){
    message("%%% CellChat error %%%")
  })
  
  # save info
  saveRDS(info, paste0(results_path, "_info.rds"))
  Data@tools$info <- info
  
  # save Data
  saveRDS(Data, paste0(results_path, "_results.rds"))
  cat("\n %%%%% analysis finished %%%%% \n")
  sink()
  return(Data)
}
