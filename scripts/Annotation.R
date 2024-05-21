Annotation <- function(Data, metaname, ref_singler_dir = NULL, ref_cell_dex = NULL, ref_markers = NULL){
  metadata_name <- NULL
  if (!is.null(ref_singler_dir)) {
    ref <- ref_cell_dex
    ref_celldex <- readRDS(paste0(ref_singler_dir, ref, ".rds"))
    cat("celldex:",ref)
    metadata_name <- append(
      metadata_name,
      values = paste0("singler_by_", c("cell","cluster"))
    )
    
    # cell annotation-----------------------------------------------------------
    cat("\n %%%%% run SingleR by cell %%%%% \n")
    cell_annotation <- SingleR(
      test = Data@assays$RNA@data,
      ref = ref_celldex,
      assay.type.test = 1, 
      labels = ref_celldex$label.fine
    )
    Data <- AddMetaData(
      Data,metadata = cell_annotation$labels,col.name = metadata_name[1]
    )
    
    # run cluster annotation----------------------------------------------------
    cat("\n %%%%% run SingleR by cluster: main %%%%% \n")
    cat("cluster:",metaname)
    cluster_annotation <- SingleR(
      test = Data@assays$RNA@data,
      ref = ref_celldex,
      clusters = Data@meta.data[,metaname],
      assay.type.test = 1,
      labels = ref_celldex$label.main
    )
    clusters <- Data@meta.data[,metaname]
    levels(clusters) <- cluster_annotation$labels
    Data <- AddMetaData(
      Data,
      metadata = clusters,
      col.name = paste0(metadata_name[2],"_main")
    )
    
    cat("\n %%%%% run SingleR by cluster: fine %%%%% \n")
    cluster_annotation <- SingleR(
      test = Data@assays$RNA@data,
      ref = ref_celldex,
      clusters = Data@meta.data[,metaname], # 按类群注释
      assay.type.test = 1,
      labels = ref_celldex$label.fine
    )
    clusters <- Data@meta.data[,metaname]
    levels(clusters) <- cluster_annotation$labels
    Data <- AddMetaData(
      Data,
      metadata = clusters,
      col.name = paste0(metadata_name[2],"_fine")
    )
    # save
    # saveRDS(subset(Data@meta.data, select = metadata_name),
    #         file = annotation_dir)
    
  }
  
  if (!is.null(ref_markers)) {
    print("annotation by markers----------------------------------------------")
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
        Data@tools$info$results_dir,
        "/DifferentialExpressing/",
        Data@tools$info$data_name,
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
  
  # saveRDS(Data,file = paste0(results_dir, "_annotation.rds"))
  cat("\n --------------------Annotation finished-------------------- \n")
  return(Data)
}

Annotation_markers <- function(seurat_object, markerList) {
  if(sum(duplicated(unlist(markerList))) > 0)
    stop("gene duplicate in different cell type")
  
  genes <- rownames(seurat_object)
  markerList <- lapply(markerList,function(x) {
    x[x %in% genes]
  })
  markerData <- seurat_object[genes %in% unlist(markerList),]
  seurat_clusters <- markerData$seurat_clusters
  Clusters <- levels(seurat_clusters)
  markers <- rownames(markerData)
  
  R <- do.call(
    rbind,
    lapply(Clusters, function(x) seurat_clusters == x)
  )
  rownames(R) <- Clusters
  
  # delete markers without differential expression in clusters
  meanmatrix <- GetAssayData(markerData,slot = "data") %*% t(1/rowSums(R) * R)
  mean_range <- as.numeric(matrix(c(-1,1),1,2) %*% apply(meanmatrix, 1, range))
  names(mean_range) <- markers
  markers <- markers[-which(mean_range < 0.5)]
  markerData <- markerData[markers,]
  markerList <- lapply(markerList,function(x) {
    x[x %in% markers]
  })
  markerList <- markerList[which(sapply(markerList, length) > 0)]
  
  L <- do.call(rbind,lapply(markerList,function(x) markers %in% x))
  
  meanmatrix <- 1/rowSums(L) * L %*% GetAssayData(markerData,slot = "data") %*% t(1/rowSums(R) * R)
  
  cellpercentmatrix <- 1/rowSums(L) * L %*% (GetAssayData(markerData,slot = "data")>0) %*% t(1/rowSums(R) * R)
  
  a <- apply(meanmatrix, 1, function(x) {
    meansort <- sort(x,decreasing = TRUE)
    (meansort[1]-meansort[2:3])/(meansort[1]-tail(meansort,1))
  })
  a <- matrix(
    t(matrix(c(0.5,0.5),1,2) %*% a),
    length(markerList),
    length(Clusters)
  )
  dimnames(a) <- list(rownames(meanmatrix),Clusters)
  
  CellTypeScoreMatrix <- 0.2*a+0.3*cellpercentmatrix+0.5*meanmatrix/10
  ClusterCellType <- apply(CellTypeScoreMatrix, 2, function(x) {names(which.max(x))})
  gc()
  return(ClusterCellType)
}
