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
  
  # 气泡图-----
  Data.all.markers %>%
    group_by(cluster) %>%
    slice_max(n = 3, order_by = avg_log2FC) -> cluster_max
  plot_de_dot <- DotPlot(
    object = Data, 
    features = unique(cluster_max$gene), 
    group.by = unique(cluster_max$cluster)
  )+coord_flip()+
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
      plot.background = element_rect(fill = "white")
    )
  ggsave(paste0(results_dir,"_differential_expression_genes_dots.",figure_format),
         plot_de_dot,width = 9,height = 9)
  
  # # 小提琴图-----
  Data.all.markers %>%
    group_by(cluster) %>%
    slice_max(n = 3, order_by = avg_log2FC) -> cluster_max
  ncluster <- length(unique(cluster_max$cluster))
  for (i in 1:ncluster) {
    genes <- cluster_max$gene[cluster_max$cluster == unique(cluster_max$cluster)[i]]
    assign(paste0("plot_de_violin_",i),
           VlnPlot(Data, features = genes, slot = "counts", log = TRUE)+
             plot_annotation(title = paste0("cluster_",i))
             )
  }
  plot_de_violin <- eval(parse(text = paste("plot_de_violin_",1:ncluster,sep = "",collapse = "/")))
  ggsave(paste0(results_dir,"_differential_expression_genes.",figure_format),
         plot_de_violin,width = 9,height = 3*ncluster)
  
  # 差异基因umap图-----
  for (i in 1:ncluster) {
    genes <- cluster_max$gene[cluster_max$cluster == unique(cluster_max$cluster)[i]]
    assign(paste0("plot_de_umap_",i),
           FeaturePlot(Data, 
                       features = genes,
                       reduction = "umap",
                       ncol = 3)
    )
  }
  plot_de_umap <- eval(parse(text = paste("plot_de_umap_",
                                          1:ncluster,sep = "",collapse = "/")))
  ggsave(paste0(results_dir,"_de_genes_umap.",figure_format),
         plot_de_umap,width = 12,height = 3*ncluster,limitsize = FALSE)
  
  # 热图-----
  Data.all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  plot_heatmap <- DoHeatmap(Data, features = top10$gene) + NoLegend()
  ggsave(paste0(results_dir,"_heatmap.",figure_format),
         plot_heatmap,width = 6,height = 3,scale = 3)
  
  # 火山图-----
  # volcano figure
  cut_off_pvalue = 1e-4
  cut_off_logFC = 1
  Data.all.markers$Sig = ifelse(
    Data.all.markers$p_val < cut_off_pvalue & 
      abs(Data.all.markers$avg_log2FC) >= cut_off_logFC, 
    ifelse(Data.all.markers$avg_log2FC > cut_off_logFC ,'Up','Down'),
    'None'
  )
  clusters <- unique(Data.all.markers$cluster)
  ncluster <- length(clusters)
  for (i in 1:ncluster) {
    cluster_markers <- subset(Data.all.markers,
                              cluster == clusters[i],
                              select = c(avg_log2FC,p_val_adj,Sig,gene)
    )
    markers_index <- subset(cluster_markers,
                            cluster_markers$p_val_adj < 0.05 &
                              abs(cluster_markers$avg_log2FC) >= 1)
    markers_index <- union(order(abs(markers_index$avg_log2FC),decreasing = T)[1:5],
                           order(markers_index$p_val_adj,decreasing = F)[1:5])
    markers_top <- cluster_markers[markers_index,]
    assign(paste0("plot_volcano_",i),
           ggplot(cluster_markers, 
                  aes(x = avg_log2FC, y = -log10(p_val_adj), 
                      colour = Sig)) +
             geom_point(alpha=0.4, size=3.5) +
             scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
             # 辅助线
             geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
             geom_hline(yintercept = -log10(cut_off_pvalue),
                        lty=4,col="black",lwd=0.8) +
             # 坐标轴
             labs(x="log2(Fold Change)",
                  y="-log10 (P-value)")+
             theme_bw()+
             ggtitle(paste0("Cluster_",clusters[i]))+
             # 图例
             theme(plot.title = element_text(hjust = 0.5), 
                   legend.position="right", 
                   legend.title = element_blank()
             )+
             theme(
               plot.background = element_rect(fill = "white")
             )+
             # 基因标签
             ggrepel::geom_text_repel(
               data = cluster_markers[markers_index,],
               aes(label = gene),
               size = 3,
               box.padding = unit(0.5, "lines"),
               point.padding = unit(0.8,"lines"),
               segment.color = "black",
               show.legend = FALSE
             )
    )
  }
  plot_text <- paste0("plot_volcano <- ggarrange(",
                      paste("plot_volcano_",1:ncluster,
                            sep = "",collapse = ", "),
                      ",ncol = 3,nrow = ",ceiling(length(clusters)/3),
                      ",common.legend = TRUE, legend = \"right\")"
  )
  eval(parse(text = plot_text))
  ggsave(paste0(results_dir,"_volcano.",figure_format),
         plot_volcano,width = 18,height = 6*ceiling(length(clusters)/3),
         scale = 1)
  
  cat("----------DifferentialExpressing finished----------")
  return(Data)
}