analysis_scrnaseq <- function(
    Data,
    pc_num = 50, resolution = 0.8,
    ref_singler_dir = NULL, ref_cell_dex = NULL, ref_markers = NULL,
    ncl = 1,
    run_enrich =TRUE, run_TA = TRUE, run_CC = TRUE
  ){
  info <- Data@tools$info
  data_name <- info$data_name
  
  # create rds dir
  info$dir$rds <- paste0(data_name,"_rds_files/")
  rds_dir <- paste0(info$dir$dir, info$dir$rds)
  if (!dir.exists(rds_dir)) dir.create(rds_dir)
  
  # create result direction
  info$dir$results <- paste0(data_name, "_results/")
  results_dir <- paste0(info$dir$dir, info$dir$results)
  if (!dir.exists(results_dir)) dir.create(results_dir)
  
  results_path <- paste0(results_dir,data_name,"_")
  sink(file = paste0(results_path, "analysis_log.txt"),split = TRUE)
  
  info$pc_num <- pc_num
  info$ref_singler_dir <- ref_singler_dir
  info$ref_cell_dex <- ref_cell_dex
  info$ref_markers <- ref_markers
  
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

  nsample <- length(unique(Data$orig.ident))
  # 判断是否进行批次效应
  if (nsample == 1) {
    cat("\n 1个样本 \n")
    reduction_used <- "pca"
  } else {
    # run correct batch
    cat("\n ", nsample ,"个样本 \n 矫正批次效应 \n")
    cat("\n %%%%% run Harmony %%%%% \n")
    Data <- suppressWarnings(harmony::RunHarmony(Data, "orig.ident"))
    reduction_used <- "harmony"
  }
  
  cat("\n %%%%% run UMAP %%%%% \n")
  Data <- RunUMAP(Data, dims = 1:pc_num, reduction = reduction_used)
  cat("\n %%%%% run TSNE %%%%% \n")
  Data <- RunTSNE(Data, dims = 1:pc_num, reduction = reduction_used)
  
  # run Cluster-----------------------------------------------------------------
  cat("\n %%%%% run FindNeighbors %%%%% \n")
  Data <- FindNeighbors(Data, dims = 1:pc_num, reduction = reduction_used)
  cat("\n %%%%% run FindClusters %%%%% \n")
  # select the number of clusters
  seq_res <- seq(0.3,1.5,0.1)
  Data <- FindClusters(
    Data, resolution = seq_res, verbose = FALSE
  )
  clustree_plt <- clustree::clustree(
    Data, prefix = paste0(DefaultAssay(Data),"_snn_res.")
  )
  stability <- sapply(
    X = 1:13,
    FUN = function(i) {
      clustree_plt$data %>% 
        filter(RNA_snn_res. == seq_res[i]) %>% 
        select(sc3_stability) %>% 
        sum
    }
  )
  ggsave(
    paste0(results_dir,info$data_name,"_clustree.svg"),
    clustree_plt,
    width = 10,height = 10
  )
  rm(clustree_plt)
  resolution <- seq_res[which.max(stability)]
  info$cluster_resolution <- resolution
  res <- paste0("RNA_snn_res.", resolution)
  Data <- AddMetaData(
    Data,
    metadata = Data@meta.data[,res],
    col.name = "seurat_clusters"
  )
  Idents(Data) <- Data$seurat_clusters
  
  # run annotation--------------------------------------------------------------
  tryCatch({
    cat("\n %%%%% run annotation %%%%% \n")
    Data <- Annotation(
      Data = Data,
      ref_singler_dir = ref_singler_dir,
      ref_cell_dex = ref_cell_dex,
      ref_markers = ref_markers
    )
    
    info$ref_cell_dex <- ref_cell_dex
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
    xlsx_name <- paste0(results_path,"celltype_count.xlsx")
    writexl::write_xlsx(
      list(count = celltype_count, percent = celltype_percent),
      path = xlsx_name
    )
    rm(celltype_percent,celltype_count,xlsx_name)
  }, error = function(e){
    message("%%% annotation error %%%")
  })
  
  
  # run find marker-------------------------------------------------------------
  tryCatch({
    cat("\n %%%%% run FindAllMarkers by seurat cluster %%%%% \n")
    all.markers <- list()
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
    
    info$filename$all_markers <- paste0(data_name,"_all.markers.rds")
    all_markers_path <- paste0(rds_dir,info$filename$all_markers)
    saveRDS(all.markers, all_markers_path)
    cat("\nsave all.markers to: \n", all_markers_path)
    rm(all_markers_path)
    
    writexl::write_xlsx(
      all.markers,
      path = paste0(results_path,"_all.markers.xlsx")
    )
    
    Idents(Data) <- Data$seurat_clusters
    
    # save DEG list
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
    info$filename$deg <- paste0(data_name,"_deg.rds")
    deg_path <- paste0(rds_dir,info$filename$deg)
    saveRDS(genes, deg_path)
    cat("\nsave deg to: \n", deg_path)
    rm(deg_path)
  }, error = function(e){
    run_enrich <- FALSE
    run_TA <- FALSE
    run_CC <- FALSE
    message("%%% find DEG error %%%")
  })
  
  # run enrich------------------------------------------------------------------
  if(run_enrich == TRUE) {
    tryCatch({
      enrichment <- Enrich(
        marker_list = genes,
        Species = info$Species,
        ncl = ncl # 线程数
      )
      
      # save enrichment
      info$filename$enrich <- paste0(data_name,"_enrichment.rds")
      enrich_path <- paste0(rds_dir,info$filename$enrich)
      saveRDS(enrichment, enrich_path)
      cat("\nsave enrichment to: \n", enrich_path)
      rm(enrich_path)
      enrichment_dataframe <- lapply(enrichment, function(x){
        lapply(x, function(y){
          dplyr::filter(y@result, p.adjust < 0.05)
        })
      })
      enrich_dir <- paste0(results_path, "enrich_table/")
      if(!dir.exists(enrich_dir)) dir.create(enrich_dir)
      for(go in names(enrichment_dataframe)){
        writexl::write_xlsx(
          enrichment_dataframe[[go]],
          path = paste0(enrich_dir,data_name,"_",go,".xlsx")
        )
      }
      rm(enrichment,enrichment_dataframe,enrich_dir)
    }, error = function(e){
      message("%%% enrich error %%%")
    })
  }
  
  # run TrajectoryAnalysis------------------------------------------------------
  if(run_TA == TRUE) {
    info$by_cluster <- "seurat_clusters"
    tryCatch({
      cds <- TrajectoryAnalysis(
        Data = Data, by_cluster = info$by_cluster
      )
      info$filename$cds <- paste0(data_name,"_cds.rds")
      cds_path <- paste0(rds_dir,info$filename$cds)
      saveRDS(cds, cds_path)
      cat("\nsave cds to: \n", cds_path)
      rm(cds, cds_path);gc()
    }, error = function(e){
      message("%%% TrajectoryAnalysis error %%%")
    })
  }
  
  # CellCommunication-----------------------------------------------------------
  if(run_CC == TRUE) {
    tryCatch({
      cellchat <- CellCommunication(
        Data = Data, by_cluster = "singler_by_cluster"
      )
      info$filename$cellchat <- paste0(data_name,"_cellchat.rds")
      cellchat_path <- paste0(rds_dir,info$filename$cellchat)
      saveRDS(cellchat, cellchat_path)
      cat("\nsave cellchat to: \n", cellchat_path)
      rm(cellchat, cellchat_path);gc()
    }, error = function(e){
      message("%%% CellChat error %%%")
    })
  }
  
  # save info
  Data@tools$info <- info
  info_path <- paste0(info$dir$dir,info$dir$seurat,info$filename$info)
  saveRDS(info, info_path)
  cat("\nsave info to: \n", info_path)
  # save result Data
  info$filename$results <- paste0(data_name, "_results_seurat.rds")
  results_data_path <- paste0(rds_dir, info$filename$results)
  saveRDS(Data, results_data_path)
  cat("\nsave results data to: \n", results_data_path)
  
  cat("\n %%%%% all analysis finished %%%%% \n")
  sink()
  return(Data)
}
