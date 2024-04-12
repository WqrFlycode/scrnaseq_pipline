# seuratWrappers
TrajectoryAnalysis <- function(Data, by_cluster = "singler_by_cluster"){
  cat("\n %%%%% transfer seurat to cds %%%%% \n")
  cds <- SeuratWrappers::as.cell_data_set(Data)
  cat("\n %%%%% cluster_cells %%%%% \n")
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cat("\n %%%%% learn_graph %%%%% \n")
  cds <- learn_graph(cds,use_partition = TRUE)
  cat("\n %%%%% Order the cells in pseudotime %%%%% \n")
  cds <- order_cells(cds, root_pr_nodes = "Y_1")
  
  cat("\n ----------TrajectoryAnalysis finished---------- \n")
  return(cds)
}

# 重建cds
# TrajectoryAnalysis <- function(Data, by_cluster = "singler_by_cluster"){
#   cat("\n %%%%% create cell data set %%%%% \n")
#   data <- GetAssayData(Data, assay = "RNA", slot = "counts")
#   cell_metadata <- Data@meta.data
#   gene_annotation <- data.frame(gene_short_name = rownames(data))
#   rownames(gene_annotation) <- rownames(data)
#   cds <- new_cell_data_set(
#     data,
#     cell_metadata = cell_metadata,
#     gene_metadata = gene_annotation
#   )
#   rm(data, cell_metadata, gene_annotation)
# 
#   # cds <- as.cell_data_set(Data)
#   cat("\n %%%%% preprocess %%%%% \n")
#   # cds <- preprocess_cds(cds, num_dim = 50)
#   cat("\n %%%%% reduce_dimension %%%%% \n")
#   cds@int_colData$reducedDims$PCA <- Data@reductions$pca@cell.embeddings
#   cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
#   cat("\n %%%%% cluster_cells %%%%% \n")
#   cds <- cluster_cells(cds, reduction_method = "UMAP")
#   cat("\n %%%%% learn_graph %%%%% \n")
#   cds <- learn_graph(cds,use_partition = TRUE)
# 
#   cat("\n %%%%% Order the cells in pseudotime %%%%% \n")
#   root_cell <- levels(colData(cds)[,by_cluster])[1]
#   root_pr_nodes <- get_earliest_principal_node(
#     cds, my_select = root_cell, by_cluster = by_cluster
#   )
#   cds <- order_cells(cds, root_pr_nodes=root_pr_nodes)
# 
#   cds@assays <- NULL
# 
#   rm(root_cell, root_pr_nodes)
#   cat("\n ----------TrajectoryAnalysis finished---------- \n")
#   return(cds)
# }


get_earliest_principal_node <- function(cds, my_select, by_cluster){
  cell_ids <- which(colData(cds)[,by_cluster] == my_select)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))
  ]
  root_pr_nodes
}
