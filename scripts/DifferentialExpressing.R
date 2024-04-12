DifferentialExpressing <- function(Data, re = FALSE, figure_format = "png"){
  results_dir <- paste0(Data@tools$results_dir,"/DifferentialExpressing/")
  if(!dir.exists(results_dir)) dir.create(results_dir)
  results_dir <- paste0(results_dir,Data@tools$data_name)
  
  # 按聚类结果
  if (!file.exists(paste0(results_dir,"_all_markers_by_cluster.csv"))) {
    cat("\n %%%%% run FindAllMarkers %%%%% \n")
    Data.all.markers.by.cluster <- FindAllMarkers(Data)
    write.csv(
      Data.all.markers.by.cluster,
      file = paste0(results_dir,"_all_markers_by_cluster.csv")
    )
  }
  
  # 按注释结果
  if (file.exists(paste0(results_dir,"_all_markers_by_annotation.csv")) & re == FALSE) {
    Data.all.markers <- read.csv(paste0(results_dir,"_all_markers_by_annotation.csv"))
  }else{
    ident_index <- grep(pattern = "cluster_", x = names(Data@meta.data))
    cluster_name <- names(Data@meta.data)[ident_index]
    Idents(Data) <- Data@meta.data[,cluster_name]
    cat("\n %%%%% run FindAllMarkers %%%%% \n")
    Data.all.markers <- FindAllMarkers(Data)
    write.csv(
      Data.all.markers,
      file = paste0(results_dir,"_all_markers_by_annotation.csv")
    )
  }
  
  cat("----------DifferentialExpressing finished----------")
  return(Data)
}

# clusters <- info$celltypes$by_seurat
# lapply(1:length(clusters), function(i) {
#   c1 <- clusters[i]
#   deglist <- lapply(clusters[-i], function(c2){
#     FindMarkers(
#       Data,
#       ident.1 = c1,
#       ident.2 = c2,
#       features = Data@assays$RNA@var.features
#     )
#   })
#   names(deglist) <- paste0("c",clusters[-i])
#   deg <- lapply(deglist, function(x) {
#     x %>% filter(p_val < 0.05, p_val_adj < 0.05, avg_log2FC > 0) %>% rownames
#   })
#   common <- Reduce(intersect, deg)
# })
# 
# deg <- FindMarkers(
#   Data,
#   ident.1 = clusters[1],
#   ident.2 = clusters[2],
#   features = Data@assays$RNA@var.features
# )
