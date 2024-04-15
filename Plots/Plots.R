Plot_seurat <- function(Data){
  info <- Data@tools$info
  Data_name <- info$data_name
  plots_dir <- paste0(paste0(info$dir$dir, info$dir$results), "/plots/")
  rds_dir <- paste0(info$dir$dir, info$dir$rds)
  if(!dir.exists(plots_dir)) dir.create(plots_dir)
  
  # DimReduction----------------------------------------------------------------
  cat("\n %%%%% plot DR %%%%% \n")
  results_dir <- paste0(plots_dir, "2_DimensionalityReduction/", Data_name)
  # Demonstration of highly variable genes
  cat("\n %%%%% plot DR HVG %%%%% \n")
  top10 <- head(VariableFeatures(Data), 10)
  plot_hvg <- VariableFeaturePlot(Data) + theme(legend.position = "bottom")
  plot_hvg <- LabelPoints(
    plot = plot_hvg, points = top10, repel = TRUE,xnudge = 0,ynudge = 0
  )
  suppressWarnings(ggsave(
    paste0(results_dir,"_hvg.png"),
    plot_hvg, width = 5, height = 5, bg = "white"
  ))
  cat("\n %%%%% plot DR elbow %%%%% \n")
  plot_elbow <- ElbowPlot(Data)
  ggsave(
    paste0(results_dir,"_elbow.svg"),
    plot_elbow,
    width = 6,height = 3,bg = "white"
  )
  cat("\n %%%%% plot DR dimloading %%%%% \n")
  plot_dimloading <- VizDimLoadings(Data, dims = 1:16, reduction = "pca")
  ggsave(
    paste0(results_dir,"_dimloading.png"),
    plot_dimloading,
    width = 20,height = 20
  )
  cat("\n %%%%% plot DR pc1-pc2 %%%%% \n")
  plot_pca <- DimPlot(
    Data, reduction = "pca",group.by = "orig.ident",raster = FALSE
  )+NoLegend()
  ggsave(
    paste0(results_dir,"_pca.png"),
    plot_pca,
    width = 6,height = 5
  )
  rm(plot_dimloading, plot_elbow, plot_hvg, plot_pca, top10)
  
  # Clustering------------------------------------------------------------------
  cat("\n %%%%% plot Clustering %%%%% \n")
  results_dir <- paste0(plots_dir, "3_Clustering/")
  if(!dir.exists(results_dir)) dir.create(results_dir)
  results_dir <- paste0(results_dir, Data_name)
  
  # 判断是否进行批次效应
  if (length(unique(Data$orig.ident)) == 1) {
    # plot umap
    plot_umap <- DimPlot(Data, reduction = "umap")
    # plot_cluster <- plot_umap+plot_layout(guides = "collect")
    ggsave(
      paste0(results_dir,"_cluster_umap.png"),
      plot_umap,
      width = 3,height = 3,scale = 3
    )
    plot_tsne <- DimPlot(Data, reduction = "tsne")
    ggsave(
      paste0(results_dir,"_cluster_tsne.png"),
      plot_tsne,
      width = 3,height = 3,scale = 3
    )
  } else {
    # plot umap and tsne
    plot_samples <- DimPlot(
      Data, 
      reduction = "umap", 
      group.by = "orig.ident",
      raster = FALSE
    )
    ggsave(
      paste0(results_dir,"_after_correct_batch_effect.png"),
      plot = plot_samples, 
      width = 3, height = 3, scale = 3
    )
    plot_umap <- DimPlot(Data, reduction = "umap", raster = FALSE)
    # plot_cluster <- plot_umap+plot_layout(guides = "collect")
    ggsave(
      paste0(results_dir,"_cluster_umap.png"),
      plot_umap,
      width = 3,height = 3,scale = 3
    )
    plot_tsne <- DimPlot(Data, reduction = "tsne",raster = FALSE)
    ggsave(
      paste0(results_dir,"_cluster_tsne.png"),
      plot_tsne,
      width = 3,height = 3,scale = 3
    )
    rm(plot_samples)
  }
  rm(plot_umap,plot_tsne)
  
  # Annotation------------------------------------------------------------------
  results_dir <- paste0(plots_dir,"4_Annotation/")
  if(!dir.exists(results_dir)) dir.create(results_dir)
  results_dir <- paste0(results_dir,Data_name)
  # 合并一级细胞类型相同的细胞
  # celltypes <- levels(Data@meta.data[,metadata_name[2]])
  # levels(Data@meta.data[,metadata_name[2]]) <- sapply(
  #   strsplit(celltypes, ":"), function(x) x[1]
  # )
  # n_celltypes <- length(levels(Data@meta.data[,metadata_name[2]]))
  # 
  # 
  # plot annotation umap
  cat("\n %%%%% plot annotation by cells %%%%% \n")
  if ("singler_by_cell" %in% names(Data@meta.data)) {
    plot_umap_all <- DimPlot(
      Data,group.by = "singler_by_cell", reduction = "umap", raster = FALSE
    )+
      theme(legend.position = "bottom",legend.text = element_text(size = 8))+ 
      guides(color = guide_legend(ncol = 3, override.aes = list(size = 4)))
    ggsave(
      paste0(results_dir,"_", "singler_by_cell.png"),
      plot_umap_all, 
      width = 10, height = 10+length(unique(Data$singler_by_cell))/(3*5)
    )
    if(length(unique(Data$orig.ident)) > 1){
      plot_umap <- DimPlot(
        Data,group.by = "singler_by_cell", reduction = "umap",
        split.by = "orig.ident", ncol = 3, 
        raster = FALSE
      )+
        theme(legend.position = "bottom") + 
        guides(color = guide_legend(ncol = 3, override.aes = list(size = 5)))
      ggsave(
        paste0(results_dir,"_", "singler_by_cell", "_samples.png"),
        plot_umap,
        width = 9*3, height = 9*ceiling(length(unique(Data$orig.ident))/3)+9,
        limitsize = FALSE
      )
    }
    # plot_tsne <- DimPlot(Data,group.by = metadata_name[i], reduction = "tsne",
    #                      label = T , repel = T, label.size = 5)
    # plot_cluster <- plot_umap+plot_tsne+plot_layout(guides = "collect")
  }
  
  cat("\n %%%%% plot annotation by clusters %%%%% \n")
  if ("singler_by_cluster" %in% names(Data@meta.data)){
    plot_umap_all <- DimPlot(
      Data,group.by = "singler_by_cluster", reduction = "umap",
      label = T , repel = T, label.size = 3, 
      raster = FALSE
    )+
      theme(legend.position = "bottom",legend.text = element_text(size = 10)) + 
      guides(color = guide_legend(ncol = 3, override.aes = list(size = 4)))
    ggsave(
      paste0(results_dir,"_", "singler_by_cluster.png"),
      plot_umap_all, 
      width = 10, height = 10+length(unique(Data$singler_by_cluster))/(3*5),
      limitsize = FALSE
    )
    if(length(unique(Data$orig.ident)) > 1){
      plot_umap <- DimPlot(
        Data,group.by = "singler_by_cluster", reduction = "umap",
        split.by = "orig.ident", ncol = 2,
        label = T , repel = T, label.size = 3, 
        raster = FALSE
      )+
        theme(legend.position = "bottom")
      # plot_tsne <- DimPlot(Data,group.by = metadata_name[i], reduction = "tsne",
      #                      label = T , repel = T, label.size = 5)
      # plot_cluster <- plot_umap+plot_tsne+plot_layout(guides = "collect")
      ggsave(
        paste0(results_dir,"_", "singler_by_cluster","_samples.png"),
        plot_umap,
        width = 10*2, height = 9*ceiling(length(unique(Data$orig.ident))/2) + 2, 
        limitsize = FALSE
      )
    }
  }
  rm(plot_umap_all)
  
  ## bar-----
  # 样本和细胞类型交叉表
  mydata <- table(Data$orig.ident, Data$singler_by_cluster)
  mydata <- mydata/rowSums(mydata)
  data <- as.data.frame(mydata)
  rm(mydata)
  
  # 创建堆积图
  plot_bar <- ggplot(data=data, aes(Var1, Freq, fill=Var2)) +
    geom_bar(
      stat="identity", position="stack", color="black", 
      width=0.7, linewidth=0.25
    ) +
    labs(x = "", y = "Proportion %") +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(
      panel.background = element_rect(
        fill="white", colour="black", linewidth = 0.25
      ),
      axis.line=element_line(colour="black", linewidth = 0.25),
      axis.title=element_text(size=13, color="black"),
      axis.text = element_text(size=12, color="black"),
      axis.text.x = element_text(angle = 30,vjust = 0.5),
      # legend.position= "bottom",
      legend.title = element_blank()
    )
  n_celltypes <- length(unique(data$Var2))
  if (n_celltypes < 13) {
    plot_bar <- plot_bar+scale_fill_manual(
      values = RColorBrewer::brewer.pal(n_celltypes, "Set3")
    )
  }else{
    plot_bar <- plot_bar+scale_fill_manual(
      values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_celltypes)
    )
  }
  ggsave(
    paste0(results_dir,"_barplot.svg"),
    plot_bar,width = 3+1*length(unique(Data$orig.ident)),height = 5
  )
  rm(data, n_celltypes, plot_bar)
  
  # Markers---------------------------------------------------------------------
  allmarkers_path <- paste0(rds_dir,info$filename$all_markers)
  if(file.exists(allmarkers_path)) {
    cat("\n %%%%% plot differential expressing genes %%%%% \n")
    results_dir <- paste0(plots_dir, "5_DifferentialExpressing/")
    if(!dir.exists(results_dir)) dir.create(results_dir)
    results_dir <- paste0(results_dir,Data_name)
    
    all.markers <- readRDS(allmarkers_path)
    
    ## bubble-----
    cluster_markers <- all.markers$by.cluster
    cluster_max <- cluster_markers %>%
      filter(p_val < 0.05, p_val_adj < 0.05, avg_log2FC > 0) %>%
      group_by(cluster) %>%
      slice_max(n = 3, order_by = avg_log2FC)
    plot_de_dot <- DotPlot(
      object = Data,
      features = unique(cluster_max$gene)
    ) + coord_flip() +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        plot.background = element_rect(fill = "white")
      )
    ggsave(
      paste0(results_dir,"_DEG_dots_cluster.svg"),
      plot_de_dot,width = 9,height = 9
    )
    
    celltype_markers <- all.markers$by.celltype
    celltype_markers %>%
      group_by(cluster) %>%
      slice_max(n = 3, order_by = avg_log2FC) -> celltype_max
    plot_de_dot <- DotPlot(
      object = Data,
      features = unique(celltype_max$gene),
      group.by = "singler_by_cluster"
    )+coord_flip()+
      theme(
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        plot.background = element_rect(fill = "white")
      )
    ggsave(
      paste0(results_dir,"_DEG_dots_celltype.svg"),
      plot_de_dot,width = 9,height = 9
    )
    rm(plot_de_dot)
    
    ## violin-----
    clusters <- unique(cluster_max$cluster)
    ncluster <- length(clusters)
    for (i in 1:ncluster) {
      genes <- cluster_max$gene[cluster_max$cluster == clusters[i]]
      # assign(paste0("plot_de_violin_",i),
      #        VlnPlot(Data, features = genes, slot = "counts", log = TRUE)+
      #          patchwork::plot_annotation(title = paste0("cluster_",i))
      # )
      ggsave(
        paste0(results_dir,"_deg_violin_cluster_", i,".png"),
        VlnPlot(
          Data, features = genes, slot = "counts", log = TRUE, raster = FALSE
        ),
        width = ncluster*2,height = 5
      )
    }
    
    ## umap-----
    for (i in 1:ncluster) {
      genes <- cluster_max$gene[cluster_max$cluster == clusters[i]]
      assign(
        paste0("plot_de_umap_",i),
        FeaturePlot(
          Data,
          features = genes,
          reduction = "umap",
          ncol = 3, 
          raster = FALSE
        )
        
      )
    }
    plot_de_umap <- eval(parse(text = paste(
      "plot_de_umap_",
      1:ncluster,sep = "",collapse = "/"
    )))
    ggsave(
      paste0(results_dir,"_de_genes_umap.png"),
      plot_de_umap,width = 12,height = 3*ncluster,limitsize = FALSE
    )
    rm(
      i, plot_de_umap, genes,
      list = paste("plot_de_umap_",1:ncluster,sep = "")
    )
    
    ## heatmap-----
    plot_heatmap <- DoHeatmap(
      Data, features = cluster_max$gene, raster = FALSE
    ) + NoLegend()
    ggsave(
      paste0(results_dir,"_heatmap.png"),
      plot_heatmap,
      width = ncluster, height = 2*ncluster,
      scale = 1,
      limitsize = FALSE
    )
    rm(plot_heatmap)
    
    ## volcano-----
    # volcano figure
    cluster_markers_list <- list()
    for (i in 1:ncluster) {
      cluster_markers_list[[i]] <- subset(
        cluster_markers,
        cluster == clusters[i],
        select = c(avg_log2FC,p_val_adj,gene)
      )
      names(cluster_markers_list)[i] <- paste0("cluster_",i)
    }
    
    volcano_list <- list()
    for (i in 1:length(cluster_markers_list)) {
      volcano_list[[i]] <- plot_volcano(
        data = cluster_markers_list[[i]], 
        FC = 1, 
        PValue = 0.05,
        volcano_title = names(cluster_markers_list)[i]
      )
    }
    for (i in 1:length(volcano_list)) {
      ggsave(
        paste0(results_dir,"_volcano_cluster",i,".png"),
        volcano_list[[i]], 
        width = 6, height = 6, scale = 1, 
        bg = "white"
      )
    }
    rm(volcano_list, clusters)
  }
  rm(allmarkers_path)
  
  # Enrich----------------------------------------------------------------------
  enrichment_path <- paste0(rds_dir,info$filename$enrich)
  if(file.exists(enrichment_path)){
    cat("\n %%%%% plot Enrichment %%%%% \n")
    results_dir <- paste0(plots_dir, "6_Enrich/")
    if(!dir.exists(results_dir)) dir.create(results_dir)
    results_dir <- paste0(results_dir, Data_name)
    
    enrich_plot <- plot_enrich(enrichment_path, results_dir)
  }
  rm(enrichment_path)
  
  # Trajectory------------------------------------------------------------------
  cds_path <- paste0(rds_dir,info$filename$cds)
  if(file.exists(cds_path)){
    cat("\n %%%%% Trajectory %%%%% \n")
    results_dir <- paste0(plots_dir, "7_Trajectory/")
    if(!dir.exists(results_dir)) dir.create(results_dir)
    results_dir <- paste0(results_dir,Data_name)
    
    cds <- readRDS(cds_path)
    cat("\n %%%%% plot trojectory %%%%% \n")
    plot_trojectory_cluster <- plot_cells(
      cds,
      color_cells_by = info$by_cluster,
      # label_groups_by_cluster = FALSE,
      label_cell_groups = FALSE,
      label_leaves = FALSE,
      label_branch_points = FALSE,
      group_label_size = 4
    )
    ggsave(
      paste0(results_dir,"_trojectory_cluster.png"),
      plot_trojectory_cluster,
      width = 8, height = 5
    )
    cat("\n %%%%% plot pseudotime %%%%% \n")
    plot_trojectory_pseudotime <- plot_cells(
      cds,
      color_cells_by = "pseudotime",
      label_cell_groups = FALSE,
      # label_groups_by_cluster = FALSE,
      label_leaves = TRUE, 
      label_branch_points = TRUE,
      group_label_size = 4
    )
    ggsave(
      paste0(results_dir,"_trojectory_pseudotime.png"),
      plot_trojectory_pseudotime,
      width = 8, height = 5
    )
    rm(plot_trojectory_cluster, plot_trojectory_pseudotime)
  }
  rm(cds_path)
  
  # CellCommunication-----------------------------------------------------------
  cellchat_path <- paste0(rds_dir,info$filename$cellchat)
  if(file.exists(cellchat_path)){
    cat("\n %%%%% plot CellCommunication %%%%% \n")
    results_dir <- paste0(plots_dir, "8_CellCommunication/")
    if(!dir.exists(results_dir)) dir.create(results_dir)
    results_dir <- paste0(results_dir,Data_name)
    cellchat <- readRDS(cellchat_path)
    group_size <- as.numeric(table(cellchat@idents))
    # 细胞通讯网络圈图
    filename_interactions_n <- paste0(results_dir,"_interactions_number.png")
    filename_interactions_w <- paste0(results_dir,"_interactions_weights.png")
    tryCatch({
      png(
        filename = filename_interactions_n,
        width = 1080,
        height = 1080,
        units = "px",
        bg = "white",
        res = 100
      )
      CellChat::netVisual_circle(
        cellchat@net$count,
        vertex.weight = group_size,
        weight.scale = T,
        label.edge = F,
        title.name = "number of interactions"
      )
      dev.off() # close
      
      png(
        filename = filename_interactions_w,
        width = 1080,
        height = 1080,
        units = "px",
        bg = "white",
        res = 100)
      CellChat::netVisual_circle(
        cellchat@net$weight,
        vertex.weight = group_size,
        weight.scale = T,
        label.edge = F,
        title.name = "Interaction weights/strength"
      )
      dev.off() # close
    }, error = function(e){
      dev.off()
      if(file.exists(filename_interactions_n)) file.remove(filename_interactions_n)
      if(file.exists(filename_interactions_w)) file.remove(filename_interactions_w)
      rm(filename_interactions_n, filename_interactions_w)
      message("%%% netVisual_circle error %%%")
    })
    
    # hierarchy figure
    filename_hierarchy <- paste0(results_dir,"_",1:5,"_hierarchy.svg")
    tryCatch({
      for (p in 1:5) {
        svg(
          filename = filename_hierarchy[p],
          bg = "white" ,height = 7,width = 15
          )
        CellChat::netVisual_aggregate(
          cellchat,
          signaling = cellchat@netP$pathways[p],
          vertex.receiver = 1:floor(length(info$celltypes$by_cluster)/2),
          layout =  "hierarchy",
          title.space = 3
        )
        dev.off()
      }
    }, error = function(e){
      replicate(length(dev.list()),dev.off())
      if(any(file.exists(filename_hierarchy))) {
        file.remove(filename_hierarchy[file.exists(filename_hierarchy)])
      }
      message("%%% hierarchy error %%%")
    })
    rm(filename_hierarchy)
  }
  rm(cellchat_path)
  
  cat("\n %%%%% all plots finished %%%%% \n")
}


# data: data$p_val_adj  data$avg_log2FC data$genes
plot_volcano <- function(data, FC = 1, PValue = 0.05, volcano_title) {
  if (!all(c("p_val_adj", "avg_log2FC", "gene") %in% names(data)))
    stop("colnames must contain p_val_adj, avg_log2FC, gene")
  
  # 判断每个基因的上下调
  data$sig[ (-1*log10(data$p_val_adj) < -1*log10(PValue)|data$p_val_adj=="NA")|(data$avg_log2FC < FC)& data$avg_log2FC > -FC] <- "NotSig"
  data$sig[-1*log10(data$p_val_adj) >= -1*log10(PValue) & data$avg_log2FC >= FC] <- "Up"
  data$sig[-1*log10(data$p_val_adj) >= -1*log10(PValue) & data$avg_log2FC <= -FC] <- "Down"
  
  # 绘制火山图
  volcano_plot <- ggplot(data, aes(avg_log2FC, -1*log10(p_val_adj))) +
    geom_point(aes(color = sig)) +
    labs(title="Volcano plot", x="log2 (FC)", y="-log10 (PValue)") +
    geom_hline(yintercept=-log10(PValue), linetype=2) +
    geom_vline(xintercept=c(-FC, FC), linetype=2) +
    scale_color_manual(values=c("Up" = "#ff4757", "Down" = "#546de5", "NotSig" = "#d2dae2")) +
    # 基因标签
    ggrepel::geom_text_repel(
      data = data[order(data$p_val_adj, decreasing = FALSE)[1:10], ],
      aes(label = gene),
      size = 3
    ) +
    ggtitle(volcano_title)
  return(volcano_plot)
}

plot_enrichment <- function(genes, Species){
  if(Species == "human") {
    orgdb <- "org.Hs.eg.db"
  } else if(Species == "mouse") {
    orgdb <- "org.Mm.eg.db"
  } else if(Species == "pig") {
    orgdb <- "org.Ss.eg.db"
  }
  #将SYMBOL格式转为ENSEMBL和ENTERZID格式 
  genes_entrezid <- clusterProfiler::bitr(
    genes, fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = orgdb
  )
  
  enrichment_plot <- list(
    GO_MF = NULL, GO_BP = NULL, GO_CC = NULL, kegg = NULL, reactome = NULL
  )
  cat("\n %%%%% enrich GO %%%%% \n")
  GO <- c("MF","BP","CC")
  for (i in 1:3) {
    go <- GO[i]
    enrich_go <- enrichGO(
      gene = genes_entrezid$ENTREZID,
      OrgDb = orgdb,
      keyType = "ENTREZID",
      ont = go,
      readable = TRUE
    )
    ## 条形图，按p从小到大排，绘制前20个Term
    bar_plot <- barplot(
      enrich_go, showCategory = 20, title = paste0("EnrichmentGO_",go,"_bar")
    )
    ## 点图，按富集的数从大到小的
    dot_plot <- dotplot(
      enrich_go, showCategory = 20, title = paste0("EnrichmentGO_",go,"_dot")
    )
    if (nrow(enrich_go) > 1) {
      ## acyclic_graph
      acyclic_plot <- goplot(enrich_go)
      enrichment_plot[[paste0("GO_", go)]] <- (bar_plot + dot_plot)/acyclic_plot
    }else{
      enrichment_plot[[paste0("GO_", go)]] <- (bar_plot + dot_plot)
    }
  }
  
  cat("\n %%%%% enrich KEGG %%%%% \n")
  if(Species == "human") {
    Organisms <- "hsa" # Homo sapiens (human)
  } else if(Species == "mouse") {
    Organisms <- "mmu" # Mus musculus (house mouse)
  } else if(Species == "pig") {
    Organisms <- "org.Ss.eg.db"
  }
  tryCatch({
    enrich_kegg <- enrichKEGG(
      gene = genes_entrezid$ENTREZID,
      organism = Organisms,
      pvalueCutoff = 0.05
    )
    KEGG_bar_plot <- barplot(enrich_kegg, title="Enrichment KEGG bar")
    KEGG_dot_plot <- dotplot(enrich_kegg, title="Enrichment KEGG dot")
    enrichment_plot[["kegg"]] <- KEGG_bar_plot + KEGG_dot_plot
  }, error = function(e){
    message("%%% Failed to download KEGG data %%%")
  })
  
  cat("\n %%%%% enrich Reactome %%%%% \n")
  tryCatch({
    enrich_reactome <- enrichPathway(
      gene = genes_entrezid$ENTREZID,
      pvalueCutoff = 0.05,
      organism = Species,
      readable = TRUE
    )
    reactome_bar_plot <- barplot(enrich_reactome,title="Enrichment Reactome bar")
    reactome_dot_plot <- dotplot(enrich_reactome,title="Enrichment Reactome dot")
    enrichment_plot[["reactome"]] <- reactome_bar_plot + reactome_dot_plot
  }, error = function(e){
    message("%%% Failed to enrich reactome %%%")
  })
  
  return(enrichment_plot)
}

# genes_in_pseudotime
plot_gene_pseudotime <- function(cds, gene, metadata_name, celltype, color_cells_by){
  cds_genes <- cds[gene, cds@colData[,metadata_name] %in% celltype]
  root_pr_nodes <- get_earliest_principal_node(
    cds_genes,
    my_select = celltype[1],
    by_cluster = metadata_name
  )
  cds_genes <- order_cells(cds_genes, root_pr_nodes=root_pr_nodes)
  plot_genes_in_pseudotime(
    cds_genes,
    color_cells_by = color_cells_by,
    min_expr=0.5
  )
}
# 
# plot_gene_pseudotime(
#   cds = Data@tools$cds, 
#   gene = Data@tools$all.markers.by.cluster$gene[1], 
#   metadata_name = "seurat_clusters", 
#   celltype = levels(Data$seurat_clusters)[c(2,5:7)],
#   color_cells_by = names(Data@meta.data)[9]
# )
