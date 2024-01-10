analysis_scrnaseq <- function(seurat_object_dir,
                              pc_num = 50, resolution = 0.5,
                              ref_singler_dir = NULL, ref_cell_dex = NULL, ref_markers = NULL){
  cat("\n %%%%% read data %%%%% \n")
  Data <- readRDS(seurat_object_dir)
  
  qc_index <- c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")
  # compute the proportion of mito-genes expression
  Data[["percent.mito"]] <- PercentageFeatureSet(Data, pattern = "^MT-") 
  # compute the proportion of ribo-genes expression
  Data[["percent.ribo"]] <- PercentageFeatureSet(Data, pattern = "^RP[SL]")
  # results is percentage times 100
  
  cat("\n %%%%% run quality control %%%%% \n")
  Data <- subset(
    Data,
    subset = nFeature_RNA > 100 & 
      nCount_RNA > 100 & # nCount_RNA < 1e5 &
      percent.mito < 30 &
      percent.ribo < 50,
    features = row.names(Data)[rowSums(GetAssayData(Data,slot = "counts")>0)>0]
  )
  
  # run Normalize-----
  cat("\n %%%%% run Normalize %%%%% \n")
  Data <- NormalizeData(Data)
  
  # run find VG-----
  cat("\n %%%%% run FindVariableFeatures %%%%% \n")
  Data <- FindVariableFeatures(Data)
  
  # run Scale-----
  cat("\n %%%%% run Scale %%%%% \n")
  Data <- ScaleData(Data, features = rownames(Data))
  
  
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
    Data <- FindNeighbors(Data, dims = 1:50)
    cat("\n %%%%% run FindClusters %%%%% \n")
    Data <- FindClusters(Data,resolution = resolution)
  } else {
    # run correct batch-----
    cat("\n ", nsample ,"个样本 \n 矫正批次效应 \n")
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
  }
  
  metadata_name <- NULL
  if (!is.null(ref_singler_dir)) {
    cat("\n %%%%% run SingleR %%%%% \n")
    ref <- ref_cell_dex
    
    annotation_dir <- paste0(results_dir,"_annotation_",ref,".rds")
    metadata_name <- append(
      metadata_name,
      values = paste0(c("cell_","cluster_"),ref)
    )
    
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
  saveRDS(Data, results_dir)
}