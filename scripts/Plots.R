Plot_seurat <- function(Data, figure_format = "png"){
  info <- Data@tools$info
  Results_dir <- info$results_dir
  Data_name <- info$data_name
  if(!dir.exists(paste0(Results_dir, "/plots/"))) {
    dir.create(paste0(Results_dir, "/plots/"))
  }
  sink(
    file = paste0(Results_dir, "/plots/", Data_name, "_plots_log.txt"),
    append = FALSE,
    split = TRUE
  )
  
  # QC-----
  cat("\n %%%%% plot QC %%%%% \n")
  results_dir <- paste0(Results_dir, "/plots/1_QualityControl/", Data_name)
  qc_index <- c(
    "nFeature_RNA", "nCount_RNA", "percent.mito", 
    "percent.ribo", "percent.redcell"
  )
  qc_index <- qc_index[qc_index %in% names(Data@meta.data)]
  
  # plot violion after QC
  plot_qc <- suppressWarnings(
    VlnPlot(
      Data,
      features = qc_index,
      group.by = "orig.ident",
      ncol = 4,
      raster = FALSE
    )
  )
  ggsave(
    paste0(results_dir,"_after_qc.",figure_format),
    plot_qc,
    width = 4*5*length(unique(Data$orig.ident)),
    height = 6 # *ceiling(length(qc_index)/3)
    # QC-index scatter
  )
  suppressWarnings(
    plot_FeatureScatter <- 
      FeatureScatter(
        Data,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",
        group.by = "orig.ident"
      )+NoLegend()+
      FeatureScatter(
        Data,feature1 = "nCount_RNA",feature2 = "percent.mito",
        group.by = "orig.ident"
      )+NoLegend()+
      FeatureScatter(
        Data,feature1 = "nCount_RNA",feature2 = "percent.ribo",
        group.by = "orig.ident"
      )+NoLegend()
  )
  ggsave(paste0(results_dir,"_qc_index_scatter.",figure_format),
         plot_FeatureScatter,width = 12,height = 6)
  rm(plot_FeatureScatter, plot_qc, qc_index)
  
  # DimReduction-----
  cat("\n %%%%% plot DR %%%%% \n")
  results_dir <- paste0(
    Results_dir,
    "/plots/2_DimensionalityReduction/",
    Data_name
  )
  # Demonstration of highly variable genes
  cat("\n %%%%% plot DR HVG %%%%% \n")
  top10 <- head(VariableFeatures(Data), 10)
  plot_hvg <- VariableFeaturePlot(Data) + theme(legend.position = "bottom")
  plot_hvg <- LabelPoints(
    plot = plot_hvg, points = top10, repel = TRUE,xnudge = 0,ynudge = 0
  )
  suppressWarnings(ggsave(
    paste0(results_dir,"_hvg.",figure_format),
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
  ggsave(paste0(results_dir,"_dimloading.",figure_format),plot_dimloading,width = 20,height = 20)
  cat("\n %%%%% plot DR pc1-pc2 %%%%% \n")
  plot_pca <- DimPlot(
    Data, reduction = "pca",group.by = "orig.ident",raster = FALSE
  )+NoLegend()
  ggsave(paste0(results_dir,"_pca.",figure_format),plot_pca,width = 6,height = 5)
  ## 热图-----
  # cat("\n %%%%% plot DR heatmap %%%%% \n")
  # png( 
  #   filename = paste0(results_dir,"_dimheatmap.",figure_format),
  #   width = 1920,
  #   height = 1080,
  #   units = "px",
  #   bg = "white",
  #   res = 100)
  # DimHeatmap(Data, dims = 1:9, cells = 500)
  # dev.off()
  rm(plot_dimloading, plot_elbow, plot_hvg, plot_pca, top10)
  
  # Clustering-----
  cat("\n %%%%% plot Clustering %%%%% \n")
  results_dir <- paste0(Results_dir,"/plots/3_Clustering/")
  if(!dir.exists(results_dir)) dir.create(results_dir)
  results_dir <- paste0(results_dir, Data_name)
  
  # 判断是否进行批次效应
  if (length(unique(Data$orig.ident)) == 1) {
    # plot umap
    plot_umap <- DimPlot(Data, reduction = "umap")
    # plot_cluster <- plot_umap+plot_layout(guides = "collect")
    ggsave(
      paste0(results_dir,"_cluster_umap.",figure_format),
      plot_umap,
      width = 3,height = 3,scale = 3
    )
    plot_tsne <- DimPlot(Data, reduction = "tsne")
    ggsave(
      paste0(results_dir,"_cluster_tsne.",figure_format),
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
      paste0(results_dir,"_after_correct_batch_effect.",figure_format),
      plot = plot_samples, 
      width = 6, height = 6, scale = 3
    )
    plot_umap <- DimPlot(Data, reduction = "umap", raster = FALSE)
    # plot_cluster <- plot_umap+plot_layout(guides = "collect")
    ggsave(
      paste0(results_dir,"_cluster_umap.",figure_format),
      plot_umap,
      width = 3,height = 3,scale = 3
    )
    plot_tsne <- DimPlot(Data, reduction = "tsne",raster = FALSE)
    ggsave(
      paste0(results_dir,"_cluster_tsne.",figure_format),
      plot_tsne,
      width = 3,height = 3,scale = 3
    )
    rm(plot_samples)
  }
  rm(plot_umap,plot_tsne)
  
  # Annotation-----
  results_dir <- paste0(Results_dir,"/plots/4_Annotation/")
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
      paste0(results_dir,"_", "singler_by_cell", ".png"),
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
        paste0(results_dir,"_", "singler_by_cell", "_samples.",figure_format),
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
      paste0(results_dir,"_", "singler_by_cluster", ".png"),
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
        paste0(results_dir,"_", "singler_by_cluster","_samples.",figure_format),
        plot_umap,
        width = 10*2, height = 9*ceiling(length(unique(Data$orig.ident))/2) + 2, 
        limitsize = FALSE
      )
    }
  }
  rm(plot_umap_all)
  
  ## 堆积图---------------------------------------------------------------------
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
  cat("\n %%%%% plot differential expressing genes %%%%% \n")
  results_dir <- paste0(Results_dir,"/plots/5_DifferentialExpressing/")
  if(!dir.exists(results_dir)) dir.create(results_dir)
  results_dir <- paste0(results_dir,Data_name)
  Data.all.markers <- Data@tools$all.markers$by.cluster
  ## 气泡图-----
  Data.all.markers %>%
    group_by(cluster) %>%
    slice_max(n = 3, order_by = avg_log2FC) -> cluster_max
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
  plot_de_dot <- DotPlot(
    object = Data,
    features = unique(cluster_max$gene),
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
  
  ## 小提琴图-----
  # Data.all.markers %>%
  #   group_by(cluster) %>%
  #   slice_max(n = 3, order_by = avg_log2FC) -> cluster_max
  clusters <- unique(cluster_max$cluster)
  ncluster <- length(clusters)
  for (i in 1:ncluster) {
    genes <- cluster_max$gene[cluster_max$cluster == clusters[i]]
    # assign(paste0("plot_de_violin_",i),
    #        VlnPlot(Data, features = genes, slot = "counts", log = TRUE)+
    #          patchwork::plot_annotation(title = paste0("cluster_",i))
    # )
    ggsave(
      paste0(results_dir,"_deg_violin_cluster_", i,".",figure_format),
      VlnPlot(
        Data, features = genes, slot = "counts", log = TRUE, raster = FALSE
      ),
      width = ncluster,height = 3
    )
  }
  # plot_de_violin <- eval(parse(
  #   text = paste("plot_de_violin_",1:ncluster,sep = "",collapse = "/")
  # ))
  # ggsave(paste0(results_dir,"_differential_expression_genes.",figure_format),
  #        plot_de_violin,width = 9,height = 3*ncluster)
  # rm(plot_de_violin, genes, list = paste("plot_de_violin_",1:ncluster,sep = ""))
  
  ## umap图-----
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
    paste0(results_dir,"_de_genes_umap.",figure_format),
    plot_de_umap,width = 12,height = 3*ncluster,limitsize = FALSE
  )
  rm(
    i, plot_de_umap, genes, cluster_max,
    list = paste("plot_de_umap_",1:ncluster,sep = "")
  )
  
  ## 热图-----
  Data.all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  plot_heatmap <- DoHeatmap(Data, features = top10$gene, raster = FALSE) + 
    NoLegend()
  ggsave(
    paste0(results_dir,"_heatmap.",figure_format),
    plot_heatmap,
    width = ncluster, height = 2*ncluster,
    scale = 1,
    limitsize = FALSE
  )
  rm(plot_heatmap, top10)
  
  ## 火山图-----
  # volcano figure
  names(Data.all.markers)[7] <- "genes"
  clusters <- unique(Data.all.markers$cluster)
  ncluster <- length(unique(Data.all.markers$cluster))
  cluster_markers_list <- list()
  for (i in 1:ncluster) {
    cluster_markers_list[[i]] <- subset(
      Data.all.markers,
      cluster == clusters[i],
      select = c(avg_log2FC,p_val_adj,genes)
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
      paste0(results_dir,"_volcano_cluster",i,".",figure_format),
      volcano_list[[i]], 
      width = 6, height = 6, scale = 1, 
      bg = "white"
    )
  }
  rm(volcano_list, clusters)
  
  # Enrich----------------------------------------------------------------------
  cat("\n %%%%% plot Enrichment %%%%% \n")
  results_dir <- paste0(Results_dir, "/plots/6_Enrich/")
  if(!dir.exists(results_dir)) dir.create(results_dir)
  results_dir <- paste0(results_dir, Data_name)
  for (i in 1:ncluster) {
    enrichment_plot <- plot_enrichment(
      genes = filter(
        cluster_markers_list[[i]],
        p_val_adj < 0.05,
        avg_log2FC > 1 | avg_log2FC < -1
      )$genes,
      Species = info$Species
    )
    for (j in 1:5) {
      db <- names(enrichment_plot)
      if(!is.null(enrichment_plot[[j]])){
        enrichment_plot[[j]] <-  enrichment_plot[[j]] + 
          patchwork::plot_annotation(paste0("cluster_", i)) & 
          theme(plot.title = element_text(hjust = 0.5))
        tryCatch({
          ggsave(
            paste0(
              results_dir, "_", names(cluster_markers_list)[i], "_", db[j], 
              ".svg"
            ),
            enrichment_plot[[j]],
            width = 12,height = 6,scale = 3,
            bg = "white"
          )
        }, error = function(e){
          message("%%% failed to plot ", db[j], " %%%")
        })
      }
    }
  }
  rm(Orgdb, db, enrichment_plot, Data.all.markers, cluster_markers_list)
  
  # Trajectory------------------------------------------------------------------
  if("cds" %in% names(Data@tools)){
    cat("\n %%%%% Trajectory %%%%% \n")
    results_dir <- paste0(Results_dir,"/plots/7_Trajectory/")
    if(!dir.exists(results_dir)) dir.create(results_dir)
    results_dir <- paste0(results_dir,Data_name)
    cds <- Data@tools$cds
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
      paste0(results_dir,"_trojectory_cluster.",figure_format),
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
      paste0(results_dir,"_trojectory_pseudotime.",figure_format),
      plot_trojectory_pseudotime,
      width = 8, height = 5
    )
    rm(plot_trojectory_cluster, plot_trojectory_pseudotime)
  }
  
  # CellCommunication-----------------------------------------------------------
  if("cellchat" %in% names(Data@tools)){
    cat("\n %%%%% plot CellCommunication %%%%% \n")
    results_dir <- paste0(Results_dir,"/plots/8_CellCommunication/")
    if(!dir.exists(results_dir)) dir.create(results_dir)
    results_dir <- paste0(results_dir,Data_name)
    cellchat <- Data@tools$cellchat
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
    
    # signal from a type of cell
    tryCatch({
      mat <- cellchat@net$count
      # plotdim <- c(ceiling(nrow(mat)/3),3)
      # par(mfrow = plotdim, xpd=TRUE, mar = c(1,1,1,1))
      filename_each_cell <- paste0(results_dir,"_each_cell_interaction_",i,".png")
      for (i in 1:nrow(mat)) {
        png(
          filename = filename_each_cell,
          width = 540,
          height = 540,
          units = "px",
          bg = "white"
          # res = 450
        )
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i,] <- mat[i,]
        netVisual_circle(
          mat2,
          vertex.weight = group_size,
          weight.scale = TRUE,
          label.edge = F,
          # arrow.width = 0.2,
          # arrow.size = 0.1,
          edge.weight.max = max(mat),
          title.name = rownames(mat)[i]
        )
        dev.off()
      }
    }, error = function(e){
      dev.off()
      if(file.exists(filename_each_cell)) file.remove(filename_each_cell)
      rm(filename_each_cell)
      message("%%% cell type circle error %%%")
    })
    
    # hierarchy figure
    filename_hierarchy <- paste0(results_dir,"_hierarchy.png")
    tryCatch({
      png(
        filename = filename_hierarchy,
        width = 1920*2,
        height = 1920*(3/4)*2,
        units = "px",
        bg = "white",
        res = 450)
      CellChat::netVisual_aggregate(
        cellchat,
        signaling = cellchat@netP$pathways,
        vertex.receiver = info$celltypes$by_seurat,
        layout = "hierarchy"
      )
      dev.off()
    }, error = function(e){
      dev.off()
      if(file.exists(filename_hierarchy)) file.remove(filename_hierarchy)
      rm(filename_hierarchy)
      message("%%% hierarchy error %%%")
    })
    
    # rm(mat, mat2, clusters, ncluster, plotdim, i, group_size)
    sink()
  }
  
  cat("\n %%%%% all plots finished %%%%% \n")
  Data@tools$all.markers <- NULL
  Data@tools$cds <- NULL
  Data@tools$cellchat <- NULL
  saveRDS(Data, paste0(info$results_dir, info$data_name, "_results.rds"))
}


# data: data$p_val_adj  data$avg_log2FC data$genes
plot_volcano <- function(data, FC = 1, PValue = 0.05, volcano_title) {
  if (!all(c("p_val_adj", "avg_log2FC", "genes") %in% names(data)))
    stop("colnames must contain p_val_adj, avg_log2FC, genes ")
  
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
      aes(label = genes),
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
    ## 条形图
    bar_plot <- barplot(
      enrich_go, showCategory=20, title = paste0("EnrichmentGO_",go,"_bar")
    ) #条状图，按p从小到大排，绘制前20个Term
    ## 点图
    dot_plot <- dotplot(
      enrich_go,title=paste0("EnrichmentGO_",go,"_dot")
    )#点图，按富集的数从大到小的
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
