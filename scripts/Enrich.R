Enrich <- function(data_name, results_dir, figure_format = "png"){
  dir_enrich <- paste0(results_dir,"/Enrich/")
  if (!file.exists(dir_enrich))  dir.create(dir_enrich)
  Markers <- read.csv(
    paste0(
      results_dir,
      "/DifferentialExpressing/",
      data_name,"_all_markers_by_annotation.csv"
    )
  )
  Markers %>%
    group_by(cluster) %>%
    slice_max(n = 100, order_by = avg_log2FC) -> cluster_max
  genes <- Markers$gene
  genes <- as.character(genes)
  
  #将SYMBOL格式转为ENSEMBL和ENTERZID格式 
  genes_entrezid <- bitr(genes, fromType = "SYMBOL",
                         toType = c("ENSEMBL", "ENTREZID"),
                         OrgDb = org.Hs.eg.db)
  
  # GO-----
  GO <- c("MF","BP","CC")
  for (i in 1:3) {
    go <- GO[i]
    go_dir <- paste0(dir_enrich,data_name,"_",go,".rds")
    if (file.exists(go_dir)) {
      enrich_go <- readRDS(go_dir)
    }else{
      enrich_go <- enrichGO(gene = genes_entrezid$ENTREZID,
                            OrgDb = "org.Hs.eg.db",
                            keyType = "ENTREZID",
                            ont = go,
                            readable = TRUE)
      # saveRDS(enrich_go,file = go_dir)
    }
    ## 条形图
    plot_bar <- barplot(enrich_go, showCategory=20,title = paste0("EnrichmentGO_",go,"_bar")) #条状图，按p从小到大排，绘制前20个Term
    ## 点图
    plot_dot <- dotplot(enrich_go,title=paste0("EnrichmentGO_",go,"_dot"))#点图，按富集的数从大到小的
    plot_go <- plot_bar+plot_dot
    ggsave(paste0(dir_enrich,data_name,"_",go,".",figure_format),
           plot_go,width = 6,height = 3,scale = 3)
    ## acyclic_graph
    acyclic_graph <- goplot(enrich_go)
    ggsave(paste0(dir_enrich,data_name,"_",go,"_acyclic_graph.",figure_format),
           acyclic_graph,width = 6,height = 3,scale = 3,bg = "white")
  }
  
  # KEGG-----
  # enrich_kegg <- enrichKEGG(gene = test1$ENTREZID,
  #                  organism = 'hsa', #KEGG可以用organism = 'hsa'
  #                  pvalueCutoff = 0.05)
  # plot_KEGG_bar <- barplot(enrich_kegg,title="Enrichment KEGG bar")
  # plot_KEGG_dot <- dotplot(enrich_kegg,title="Enrichment KEGG dot")
  # plot_KEGG <- plot_KEGG_bar+plot_KEGG_dot
  # ggsave(paste0(dir_enrich,data_name,"_KEGG.",figure_format),plot_KEGG,width = 6,height = 3,scale = 3)
  
  # Reactome
  enrich_reactome <- enrichPathway(gene = genes_entrezid$ENTREZID,
                                  pvalueCutoff = 0.05,
                                  readable = TRUE)
  plot_reactome_bar <- barplot(enrich_reactome,title="Enrichment Reactome bar")
  plot_reactome_dot <- dotplot(enrich_reactome,title="Enrichment Reactome dot")
  plot_reactome <- plot_reactome_bar+plot_reactome_dot
  ggsave(paste0(dir_enrich,data_name,"_reactome.",figure_format),plot_reactome,width = 6,height = 3,scale = 3)
  
  cat("--------------------Enrich finished--------------------")
}