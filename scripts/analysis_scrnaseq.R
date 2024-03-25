analysis_scrnaseq <- function(Data,
                              min_nFeature = 100, min_nCount = 100, 
                              max_percent_mito = 30, max_percent_ribo = 50, max_percent_hb = 30,
                              pc_num = 50, resolution = 0.8,
                              ref_singler_dir = NULL, ref_cell_dex = NULL, ref_markers = NULL,
                              ncl = 1){
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
  
  # QC-----
  # qc_index <- c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")
  # compute the proportion of mito-genes expression
  if(! "percent.mito" %in% names(Data@meta.data)){
    if (info$Species == "mouse") {
      pattern_mt <- "^mt-"
    }else{
      pattern_mt <- "^MT-"
    }
    Data[["percent.mito"]] <- PercentageFeatureSet(Data, pattern = pattern_mt)
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
  
  cat("\n %%%%% run quality control %%%%% \n") # 3sigma, mad
  # box
  outlier_range <- function(x, n = 1.5){
    l <- sum(fivenum(x)*c(0,1+n,0,-n,0))
    u <- sum(fivenum(x)*c(0,-n,0,1+n,0))
    return(c(l,u))
  }
  nFeature_lu <- outlier_range(Data$nFeature_RNA,n = 3) # lower & upper
  nCount_lu <- outlier_range(Data$nCount_RNA,n = 3)
  Data <- subset(
    Data,
    subset =
      nFeature_RNA > max(1e2, nFeature_lu[1]) & nFeature_RNA < nFeature_lu[2] &
      nCount_RNA > max(1e2, nCount_lu[1]) & nCount_RNA < nCount_lu[2] &
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
    RunPCA(Data, features = VariableFeatures(object = Data))
  )
  pc_num <- as.numeric(findPC::findPC(Data@reductions$pca@stdev))
  info$pc_num <- pc_num

  cat("\n %%%%% run UMAP %%%%% \n")
  Data <- RunUMAP(Data, dims = 1:pc_num)
  cat("\n %%%%% run TSNE %%%%% \n")
  Data <- RunTSNE(Data, dims = 1:pc_num)

  # run Cluster-----------------------------------------------------------------
  cat("\n %%%%% run FindNeighbors %%%%% \n")
  Data <- FindNeighbors(Data, dims = 1:pc_num)
  cat("\n %%%%% run FindClusters %%%%% \n")
  # select the number of clusters
  seq_res <- seq(0.3,1.5,0.1)
  Data <- FindClusters(Data, resolution = seq_res, verbose = FALSE)
  clustree_plt <- clustree::clustree(
    Data, prefix = paste0(DefaultAssay(Data),"_snn_res.")
  )
  Data@tools$clustree <- clustree_plt
  cell_dists <- dist(Embeddings(Data)[,1:pc_num])
  cluster_info <- Data@meta.data[,grep("snn", names(Data@meta.data))] %>% mutate_all(as.numeric)
  silhouette_res <- apply(cluster_info, 2, function(x){
    si <- cluster::silhouette(x, cell_dists) # 轮廓系数越接近1，表示分类越正确
    mean(si[,"sil_width"])
  })
  info$silhouette_res <- silhouette_res
  Data <- AddMetaData(
    Data,
    metadata = Data@meta.data[,names(which.max(silhouette_res))],
    col.name = "seurat_clusters"
  )
  Data@meta.data <- dplyr::select(
    Data@meta.data, 
    -grep("RNA_snn_res", names(Data@meta.data))
  )
  Idents(Data) <- Data$seurat_clusters
  
  # nsample <- length(unique(Data$orig.ident))
  # 判断是否进行批次效应
  # if (nsample == 1) {
  #   cat("\n 1个样本 \n")
  # } else {
  #   # run correct batch-----
  #   cat("\n ", nsample ,"个样本 \n 矫正批次效应 \n")
  #   cat("\n %%%%% run Harmony %%%%% \n")
  #   Data <- suppressWarnings(harmony::RunHarmony(Data, "orig.ident"))
  #
  #   # run UMAP based on harmony-----
  #   cat("\n %%%%% run UMAP %%%%% \n")
  #   Data <- RunUMAP(Data, dims = 1:pc_num, reduction = "harmony")
  #   cat("\n %%%%% run TSNE %%%%% \n")
  #   Data <- RunTSNE(Data, dims = 1:pc_num, reduction = "harmony")
  #   cat("\n %%%%% run Neighbors %%%%% \n")
  #   Data <- FindNeighbors(Data, dims = 1:pc_num, reduction = "harmony")
  #   cat("\n %%%%% run Clusters %%%%% \n")
  #   Data <- FindClusters(Data, resolution = resolution)
  # }
  
  # run annotation-----
  cat("\n %%%%% run annotation %%%%% \n")
  Data <- Annotation(
    Data = Data,
    ref_singler_dir = ref_singler_dir,
    ref_cell_dex = ref_cell_dex,
    ref_markers = ref_markers
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
  
  # save celltype count to xlsx
  celltype_count <- table(Data$singler_by_cluster, Data$orig.ident)
  celltype_percent <- apply(celltype_count, 2, prop.table)
  celltype_percent <- round(celltype_percent, digits = 4)*1e2
  celltype_count <- tidyr::spread(
    as.data.frame(celltype_count), key = Var2, value = Freq
  )
  names(celltype_count)[1] <- "celltypes"
  celltype_percent <- data.frame(
    celltypes = rownames(celltype_percent),
    as.data.frame(celltype_percent)
  )
  rownames(celltype_percent) <- NULL
  xlsx_name <- paste0(results_path,"_celltype_count.xlsx")
  writexl::write_xlsx(
    list(count = celltype_count, percent = celltype_percent),
    path = xlsx_name
  )
  rm(celltype_percent,celltype_count,xlsx_name)

  # run find marker-------------------------------------------------------------
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
  
  clusters <- levels(Data$seurat_clusters)
  names(clusters) <- paste0("cluster_",clusters)
  cluster_markers_list <- lapply(clusters, function(x){
    subset(
      Data.all.markers.by.cluster,
      cluster == x,
      select = c(avg_log2FC,p_val_adj,gene)
    )
  })
  genes <- lapply(cluster_markers_list, function(x){
    dplyr::filter(
      x,
      p_val_adj < 0.05,
      avg_log2FC > 1 | avg_log2FC < -1
    )$gene
  })
  saveRDS(genes, paste0(results_path, "_deg.rds"))

  # run enrich------------------------------------------------------------------
  enrichment <- Enrich(
    marker_list = genes,
    Species = info$Species,
    ncl = ncl # 线程数
  )

  # save enrichment
  saveRDS(enrichment, paste0(results_path, "_enrichment.rds"))
  enrichment_dataframe <- lapply(enrichment, function(x){
    lapply(x, function(y){
      dplyr::filter(y@result, p.adjust < 0.05)
    })
  })
  enrich_dir <- paste0(results_path, "_enrich_table/")
  if(!dir.exists(enrich_dir)) dir.create(enrich_dir)
  for(go in names(enrichment_dataframe)){
    writexl::write_xlsx(
      enrichment_dataframe[[go]],
      path = paste0(enrich_dir,info$data_name,"_",go,".xlsx")
    )
  }
  rm(enrichment,enrichment_dataframe,enrich_dir)

  # run TrajectoryAnalysis------------------------------------------------------
  info$by_cluster <- "seurat_clusters"
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
  cat("\n %%%%% all analysis finished %%%%% \n")
  sink()
  return(Data)
}
