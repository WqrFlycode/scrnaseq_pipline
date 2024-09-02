ElementaryPipline <- function(Data, ngenes = 2000, pc_num = NULL, ncl = 1) {
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
  
  # run Normalize-----
  cat("\n %%%%% run Normalize %%%%% \n")
  Data <- NormalizeData(Data)
  
  # run find VG-----
  cat("\n %%%%% run FindVariableFeatures %%%%% \n")
  Data <- FindVariableFeatures(Data,nfeatures = ngenes)
  
  # run Scale-----
  cat("\n %%%%% run Scale %%%%% \n")
  Data <- ScaleData(Data)
  
  # run PCA----
  cat("\n %%%%% run PCA %%%%% \n")
  Data <- suppressMessages(
    RunPCA(Data, features = VariableFeatures(object = Data))
  )
  
  
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
  
  cat("\n %%%%% plot UMAP for n pc %%%%% \n")
  if(is.null(pc_num)) {
    pc_seq <- seq(10,50,10)
  } else {
    pc_seq <- pc_num
  }
  names(pc_seq) <- paste(pc_seq, "pc")
  
  cl <- makeCluster(ncl)
  clusterExport(
    cl,
    c("pc_seq", "Data","reduction_used"),
    envir = environment()
  )
  umaplist <- pbapply::pblapply(
    X = pc_seq,
    FUN = function(n) {
      Data <- Seurat::RunUMAP(Data, dims = 1:n, reduction = reduction_used)
      umapplot <- Seurat::DimPlot(
        Data, reduction = "umap",label = TRUE,raster = FALSE
      ) + ggplot2::labs(title = paste(n, "pc"))
      return(umapplot)
    },
    cl = cl
  )
  stopCluster(cl)
  pc_dir <- paste0(info$dir$dir,info$dir$results,"umap_npc/")
  mkdirs(pc_dir)
  ggsave(
    paste0(info$dir$dir,info$dir$results,"umap_npc/ElbowPlot.png"),
    ElbowPlot(Data,ndims = 50,reduction = reduction_used),
    width = 6,height = 6,bg = "white"
  )
  for(npc in pc_seq) {
    umap_path <- paste0(pc_dir,info$data_name,"_",npc,"_pc_","_umap.png")
    ggsave(
      umap_path,
      umaplist[[paste(npc,"pc")]],
      width = 6,height = 6
    )
  }
  cat("\nsave umap to\n",pc_dir)
  
  Data@tools$info <- info
  return(Data)
}

DetermineCluster <- function(Data, pc_num, cl_res = NULL, ncl = 1){
  info <- Data@tools$info
  info$pc_num <- pc_num
  if("harmony" %in% Reductions(Data)){
    reduction_used <- "harmony"
  }else{
    reduction_used <- "pca"
  }
  
  cat("\n %%%%% run UMAP %%%%% \n")
  Data <- RunUMAP(Data, dims = 1:pc_num, reduction = reduction_used)
  cat("\n %%%%% run TSNE %%%%% \n")
  Data <- RunTSNE(Data, dims = 1:pc_num, reduction = reduction_used)
  
  # run Cluster-----------------------------------------------------------------
  cat("\n %%%%% run FindNeighbors %%%%% \n")
  Data <- FindNeighbors(Data, dims = 1:pc_num, reduction = reduction_used)
  cat("\n %%%%% run FindClusters %%%%% \n")
  # select resolution
  if(is.null(cl_res)) {
    seq_res <- seq(0.1,1.5,0.1)
  } else {
    seq_res <- cl_res
  }
  names(seq_res) <-paste0("RNA_snn_res.",seq_res)
  
  cl <- makeCluster(ncl)
  clusterExport(
    cl,
    c("seq_res", "Data"),
    envir = environment()
  )
  RNA_snn_res <- pbapply::pblapply(
    X = seq_res,
    FUN = function(res) {
      Data <- Seurat::FindClusters(
        Data, resolution = res, verbose = FALSE
      )
      RNA_snn_res <- Data@meta.data[[paste0("RNA_snn_res.",res)]]
      return(RNA_snn_res)
    },
    cl = cl
  )
  stopCluster(cl)
  Data@meta.data <- data.frame(Data@meta.data,RNA_snn_res)
  rm(RNA_snn_res)
  
  res_dir <- paste0(info$dir$dir,info$dir$results,"all_res/")
  mkdirs(paste0(res_dir,"umap_res_all/"))
  clustree_plt <- clustree::clustree(
    Data, prefix = paste0(DefaultAssay(Data),"_snn_res.")
  )
  ggsave(
    paste0(res_dir, info$data_name, "_clustree.svg"),
    clustree_plt,
    width = 10,height = 10
  )
  res <- names(Data@meta.data)[grep("RNA_snn_res.",names(Data@meta.data))]
  for(i in 1:15) {
    umap_path <- paste0(res_dir,"umap_res_all/",info$data_name,"_res_",i,"_umap.png")
    Idents(Data) <- Data@meta.data[[res[i]]]
    ggsave(
      umap_path,
      DimPlot(Data, reduction = "umap",label = TRUE,raster = FALSE),
      width = 6,height = 6
    )
  }
  cat("\nsave umap to\n",res_dir)
  
  stability <- sapply(
    X = 1:13, function(i) {
      clustree_plt$data %>% 
        filter(RNA_snn_res. == seq_res[i]) %>% 
        select(sc3_stability) %>% 
        sum
    }
  )
  rm(clustree_plt)
  resolution <- seq_res[which.max(stability)]
  cat("\nAutomatic selection of resolution:", resolution)
  info$cluster_resolution <- resolution
  
  info$filename$results <- paste0(info$data_name, "_results_seurat.rds")
  # Data <- SaveData(Data, info)
  Data@tools$info <- info
  return(Data)
}

AdvancedAnalysis <- function(
    Data, resolution,
    metanames = c("seurat_clusters","singler_by_cluster_fine","singler_by_cluster_main"),
    ref_singler_dir = NULL, ref_cell_dex = NULL, ref_markers = NULL,
    ncl = 1,
    run_enrich = FALSE, run_TA = FALSE, run_CC = FALSE
  ){
  info <- Data@tools$info
  data_name <- info$data_name
  rds_dir <- paste0(info$dir$dir, info$dir$rds)
  results_dir <- paste0(info$dir$dir, info$dir$results)
  results_path <- paste0(results_dir,data_name,"_")
  
  info$ref_singler_dir <- ref_singler_dir
  info$ref_cell_dex <- ref_cell_dex
  info$ref_markers <- ref_markers
  
  # set resolution--------------------------------------------------------------
  cat("\n %%%%% set resolution %%%%% \n")
  Data <- AddMetaData(
    Data,
    metadata = Data[[paste0("RNA_snn_res.",resolution)]],
    col.name = "seurat_clusters"
  )
  Idents(Data) <- Data$seurat_clusters
  info$cluster_resolution <- resolution
  
  # run annotation--------------------------------------------------------------
  tryCatch({
    cat("\n %%%%% run annotation %%%%% \n")
    Data <- Annotation(
      Data = Data,
      metaname = "seurat_clusters",
      ref_singler_dir = ref_singler_dir,
      ref_cell_dex = ref_cell_dex,
      ref_markers = ref_markers
    )
    
    info$ref_cell_dex <- ref_cell_dex
    info$celltypes <- list(
      by_seurat = levels(Data$seurat_clusters),
      by_cell = unique(Data$singler_by_cell),
      by_cluster_main = levels(Data$singler_by_cluster_main),
      by_cluster_fine = levels(Data$singler_by_cluster_fine)
    )
    
    # save celltype count to xlsx
    CountList <- list()
    for(metaname in metanames) {
      celltype_count <- table(Data@meta.data[,metaname], Data$orig.ident)
      celltype_percent <- apply(celltype_count, 2, prop.table)
      celltype_percent <- round(celltype_percent, digits = 4)*1e2
      celltype_count <- tidyr::spread(
        as.data.frame(celltype_count), key = Var2, value = Freq
      )
      names(celltype_count)[1] <- "celltypes"
      CountList[[paste0(metaname,"_count")]] <- celltype_count
      celltype_percent <- data.frame(
        celltypes = rownames(celltype_percent),
        as.data.frame(celltype_percent)
      )
      rownames(celltype_percent) <- NULL
      CountList[[paste0(metaname,"_percent")]] <- celltype_percent
    }
    
    xlsx_name <- paste0(results_path,"celltype_count.xlsx")
    writexl::write_xlsx(
      CountList,
      path = xlsx_name
    )
    rm(celltype_percent,celltype_count,xlsx_name, CountList)
  }, error = function(e){
    run_enrich <- FALSE
    run_TA <- FALSE
    run_CC <- FALSE
    message("%%% annotation error %%%")
  })
  
  # run find marker-------------------------------------------------------------
  tryCatch({
    metanameUse <- sapply(
      metanames, 
      FUN = function(x) length(levels(Data@meta.data[[x]])) > 1
    )
    metanames <- metanames[metanameUse]
    all.markers <- FindDEG(seurat_object = Data,metanames = metanames,ncl = ncl)
    Idents(Data) <- Data$seurat_clusters
    
    info$filename$all_markers <- paste0(data_name,"_all.markers.rds")
    all_markers_path <- paste0(rds_dir,info$filename$all_markers)
    saveRDS(all.markers, all_markers_path)
    cat("\nsave all.markers to: \n", all_markers_path)
    rm(all_markers_path)
    
    writexl::write_xlsx(
      all.markers,
      path = paste0(results_path,"all.markers.xlsx")
    )
    
    # save DEG list
    clusters <- levels(Data$seurat_clusters)
    names(clusters) <- paste0("cluster_",clusters)
    cluster_markers_list <- lapply(clusters, function(x){
      subset(
        all.markers$seurat_clusters,
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
    rm(deg_path, all.markers, clusters)
  }, error = function(e){
    run_enrich <- FALSE
    message("%%% find DEG error %%%")
  })
  
  # run enrich------------------------------------------------------------------
  info$filename$enrich <- paste0(data_name,"_enrichment.rds")
  if(run_enrich == TRUE) {
    tryCatch({
      enrichment <- Enrich(
        marker_list = genes,
        Species = info$Species,
        ncl = ncl # 线程数
      )
      
      # save enrichment
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
  info$filename$cds <- paste0(data_name,"_cds.rds")
  if(run_TA == TRUE) {
    info$by_cluster <- "seurat_clusters"
    tryCatch({
      cds <- TrajectoryAnalysis(
        Data = Data, by_cluster = info$by_cluster
      )
      cds_path <- paste0(rds_dir,info$filename$cds)
      saveRDS(cds, cds_path)
      cat("\nsave cds to: \n", cds_path)
      rm(cds, cds_path);gc()
    }, error = function(e){
      message("%%% TrajectoryAnalysis error %%%")
    })
  }
  
  # CellCommunication-----------------------------------------------------------
  info$filename$cellchat <- paste0(data_name,"_cellchat.rds")
  if(run_CC == TRUE) {
    tryCatch({
      cellchat <- CellCommunication(
        Data = Data, by_cluster = "singler_by_cluster_fine"
      )
      cellchat_path <- paste0(rds_dir,info$filename$cellchat)
      saveRDS(cellchat, cellchat_path)
      cat("\nsave cellchat to: \n", cellchat_path)
      rm(cellchat, cellchat_path);gc()
    }, error = function(e){
      message("%%% CellChat error %%%")
    })
  }
  
  Data <- SaveData(Data, info)
  cat("\n %%%%% all analysis finished %%%%% \n")
  return(Data)
}

SaveData <- function(Data, info) {
  # save info
  info_path <- paste0(info$dir$dir,info$dir$rds,info$filename$info)
  saveRDS(info, info_path)
  cat("\nsave info to: \n", info_path)
  # save result Data
  Data@tools$info <- info
  results_data_path <- paste0(info$dir$dir, info$dir$rds, info$filename$results)
  saveRDS(Data, results_data_path)
  cat("\nsave results data to: \n", results_data_path)
  return(Data)
}

# AutomaticAnalysis <- funcion()

