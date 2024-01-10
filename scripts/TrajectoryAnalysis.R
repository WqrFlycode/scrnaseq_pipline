TrajectoryAnalysis <- function(Data, by_cluster = "seurat_clusters", figure_format = "png"){
  results_dir <- paste0(Data@tools$results_dir,"/TrajectoryAnalysis/",Data@tools$data_name)
  
  cat("----------create cell data set--------------------\n")
  data <- GetAssayData(Data, assay = "RNA", slot = "counts")
  cell_metadata <- Data@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  cds <- new_cell_data_set(data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  
  # cds <- as.cell_data_set(Data)
  
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- reduce_dimension(cds, preprocess_method = "PCA")
  
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds,use_partition = TRUE)
  
  plot_trojectory_cluster <- plot_cells(
    cds,
    color_cells_by = by_cluster, 
    label_groups_by_cluster = FALSE, 
    label_leaves = FALSE, 
    label_branch_points = FALSE,
    group_label_size = 4
  )
  
  cat("Order the cells in pseudotime----------\n")
  get_earliest_principal_node <- function(cds, my_select){
    cell_ids <- which(colData(cds)[,by_cluster] == my_select)
    closest_vertex <-
      cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
      igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                (which.max(table(closest_vertex[cell_ids,]))))]
    root_pr_nodes
  }
  
  root_cell <- levels(colData(cds)[,by_cluster])[1]
  root_pr_nodes <- get_earliest_principal_node(cds,my_select = root_cell)
  cds <- order_cells(cds, root_pr_nodes=root_pr_nodes)
  plot_trojectory_pseudotime <- plot_cells(cds,
                                           color_cells_by = "pseudotime", 
                                           label_groups_by_cluster = FALSE, 
                                           label_leaves = TRUE, label_branch_points = TRUE,
                                           group_label_size = 4)
  plot_trojectory <- plot_trojectory_cluster + plot_trojectory_pseudotime
  ggsave(paste0(results_dir,"_trojectory.",figure_format),
         plot_trojectory,width = 10,height = 5)
  
  # cat("pseudotime of each cell type----------\n")
  # celltypes <- levels(colData(cds)[,by_cluster])
  # n <- length(celltypes)
  # for (i in 1:n) {
  #   cell_subset <- which(colData(cds)[,by_cluster] == celltypes[i])
  #   cds_subset <- cds[,cell_subset]
  #   cds_subset <- cluster_cells(cds_subset)
  #   cds_subset <- learn_graph(cds_subset)
  #   root_pr_nodes <- get_earliest_principal_node(cds_subset,my_select = celltypes[i])
  #   cds_subset <- order_cells(cds_subset, root_pr_nodes=root_pr_nodes)
  #   assign(paste0("plot_pseudotime_",i),
  #          plot_cells(cds_subset,
  #                     color_cells_by = "pseudotime",
  #                     label_groups_by_cluster = FALSE, 
  #                     label_leaves = FALSE, label_branch_points = FALSE,
  #                     group_label_size = 4)+
  #            ggtitle(celltypes[i]))
  # }
  # plot_text <- paste0("plot_pseudotime <- ggarrange(",
  #                     paste("plot_pseudotime_",1:n,sep = "",collapse = ", "),
  #                     ",ncol = 3,nrow = ",ceiling(n/3),
  #                     ",common.legend = TRUE, legend = \"right\")"
  # )
  # eval(parse(text = plot_text))
  # ggsave(paste0(results_dir,"_pseudotime.",figure_format),
  #        plot_pseudotime,width = 3*3,height = 3*ceiling(n/3))
  
  # cat("拟时轨迹差异基因----------\n")
  # if (file.exists(paste0(results_dir,"_pseudotime.csv"))) {
  #   Track_genes <- read.csv(paste0(results_dir,"_pseudotime.csv"))
  # }else{
  #   Track_genes <- graph_test(cds, neighbor_graph =  "principal_graph", cores = 6)
  #   Track_genes <- Track_genes[,c(2,3,4,1,5)] %>% filter(q_value < 1e-3)
  #   Track_genes <- head(data.frame(gene_short_name = rownames(Track_genes), Track_genes))
  #   write.csv(Track_genes, paste0(results_dir,"_pseudotime.csv"), 
  #             row.names = F)
  # }
  # Track_genes_sig <- Track_genes %>% 
  #   top_n(n=10, morans_I) %>% 
  #   pull(gene_short_name) %>% 
  #   as.character()
  # plot_jitter <- plot_genes_in_pseudotime(
  #   cds[Track_genes_sig,], label_by_short_name = FALSE,
  #   color_cells_by = by_cluster,
  #   min_expr = 0.5,
  #   ncol = 2
  # )+guides(color=guide_legend(title = "cell_type"))
  # ggsave(
  #   paste0(results_dir,"_jitter.",figure_format),
  #   plot_jitter, width = 8, height = 6
  # )
  
  # cat("寻找共表达基因模块----------\n")
  # genelist <- pull(Track_genes, gene_short_name) %>% as.character()
  # gene_module <- find_gene_modules(cds[genelist,],resolution = 1e-2)
  
  # saveRDS(Data,file = paste0(results_dir, "_trajectory_cds.rds"))
  print("----------TrajectoryAnalysis finished----------")
  return(Data)
}