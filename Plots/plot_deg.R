DegTwo <- function(Data, id1,id2 = NULL, metaname = "orig.ident") {
  deg_results <- FindMarkers(
    Data,
    ident.1 = id1, ident.2 = id2,
    group.by = metaname
  )
  deg_results$abs_avg_log2FC <- abs(deg_results$avg_log2FC)
  deg_results$gene <- rownames(deg_results)
  info <- Data@tools$info
  deg_results_dir <- paste0(info$dir$dir,info$dir$results,"deg_results/")
  if(!dir.exists(deg_results_dir)) dir.create(deg_results_dir)
  if(is.null(id2)) id2 <- "other"
  DegTitle <- paste0(id1,"_and_",id2)
  writexl::write_xlsx(
    deg_results,
    path = paste0(deg_results_dir,DegTitle,".xlsx")
  )
  deg_results$gene <- rownames(deg_results) 
  volcano_plot <- plot_volcano(
    data = deg_results, 
    FC = 1, 
    PValue = 0.05,
    volcano_title = DegTitle
  )
  volcano_path <- paste0(deg_results_dir,id1,"_and_",id2,".png")
  ggsave(
    volcano_path,
    volcano_plot,
    width = 6, height = 6, scale = 1, 
    bg = "white"
  )
  cat("\nsave volcano to\n", volcano_path)
}

plot_heatmap <- function(Data, genes, metaname, plot_dir = NULL, slot_use = "scale.data", suffix = NULL) {
  if(!metaname %in% names(Data@meta.data)) stop(metaname, " is not in meta.data")
  heatmap_plot <- DoHeatmap(
    Data,
    features = genes,
    slot = slot_use,
    raster = FALSE,
    group.by = metaname
  ) + 
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5,size = 15)
    ) + 
    ggtitle(suffix)
  ncluster <- length(unique(Data@meta.data[,metaname]))
  
  info <- Data@tools$info
  if(is.null(plot_dir)) {
    plot_dir <- paste0(
      info$dir$dir,
      info$dir$results,
      "plots/",
      "5_DifferentialExpressing/"
    )
  }
  plot_path <- paste0(
    plot_dir,
    info$data_name,"_",metaname,"_",slot_use,"_",suffix,"_heatmap.png"
  )
  ggsave(
    plot_path,
    heatmap_plot,
    width = 3+ncluster, height = length(genes)+3,
    scale = 1,
    limitsize = FALSE
  )
  cat("\nsave heatmap to\n",plot_path)
}

# data: data$p_val_adj  data$avg_log2FC data$gene
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
      data = data[which(abs(data$avg_log2FC) > FC), ],
      aes(label = gene),
      size = 3
    ) +
    ggtitle(volcano_title)
  return(volcano_plot)
}

plot_markers <- function(Data, markerlist, metaname) {
  info <- Data@tools$info
  gene <- rownames(Data)
  markerlist <- lapply(
    markerlist,function(x) {
      x[x %in% gene]
    }
  )
  markerlist <- markerlist[sapply(markerlist,length) != 0]
  plotdir <- paste0(info$dir$dir,info$dir$results)
  # umap
  if(!dir.exists(plotdir))
    stop("dir does not exist")
  plotdir <- paste0(plotdir,info$data_name,"_markers_plots/")
  umap_path <- paste0(
    plotdir,
    info$data_name,"_",metaname,"_umap.png"
  )
  # if(!identical(Idents(Data),Data$seurat_clusters))
  #   stop("idents != seurat_clusters")
  # if(!identical(Data@meta.data[["seurat_clusters"]],Data@meta.data[[metaname]]))
  #   stop("seurat_clusters != metaname")
  Idents(Data) <- Data@meta.data[[metaname]]
  ggsave(
    umap_path,
    DimPlot(Data, reduction = "umap",label = TRUE,raster = FALSE),
    width = 3,height = 3,scale = 3
  )
  cat("\nsave umap to\n",umap_path)
  # heatmap
  for(celltype in names(markerlist)) {
    plot_heatmap(
      Data = Data,
      genes = markerlist[[celltype]],
      metaname = metaname,
      plot_dir = paste0(plotdir,"heatmap/"),
      slot_use = "data",
      suffix = celltype
    )
  }
  # gene umap
  markers <- unlist(markerlist,use.names = FALSE)
  markers <- markers[markers %in% rownames(Data)]
  for(gene in markers) {
    ggsave(
      paste0(plotdir,"markerumap/",gene,"_umap.png"),
      FeaturePlot(
        Data,
        features = gene,
        reduction = "umap",order = TRUE,pt.size = 0.01,
        raster = FALSE
      ),
      width = 5,height = 5
    )
    cat("\nsave marker umap to\n",paste0(plotdir,"markerumap/",gene,"_umap.png"))
    ggsave(
      paste0(plotdir,"markerviolin/",gene,"_violin.png"),
      VlnPlot(
        Data,
        features = gene,
        # log = TRUE,
        raster = FALSE
      ),
      width = 5,height = 5
    )
    cat("\nsave violin to\n",paste0(plotdir,"markerviolin/",gene,"_violin.png"))
  }
}
