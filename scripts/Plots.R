plot_volcano <- function(data, FC = 1.5, PValue = 0.05) {
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
    )
  return(volcano_plot)
}

plot_enrichment <- function(genes, orgdb = "org.Hs.eg.db"){
  #将SYMBOL格式转为ENSEMBL和ENTERZID格式 
  genes_entrezid <- clusterProfiler::bitr(
    genes, fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = orgdb
  )
  
  enrichment_plot <- list(
    GO_MF = NULL, GO_BP = NULL, GO_CC = NULL, kegg = NULL, reactome = NULL
  )
  # GO
  GO <- c("MF","BP","CC")
  for (i in 1:3) {
    go <- GO[i]
    enrich_go <- enrichGO(gene = genes_entrezid$ENTREZID,
                          OrgDb = orgdb,
                          keyType = "ENTREZID",
                          ont = go,
                          readable = TRUE)
    ## 条形图
    bar_plot <- barplot(enrich_go, showCategory=20,title = paste0("EnrichmentGO_",go,"_bar")) #条状图，按p从小到大排，绘制前20个Term
    ## 点图
    dot_plot <- dotplot(enrich_go,title=paste0("EnrichmentGO_",go,"_dot"))#点图，按富集的数从大到小的
    if (nrow(enrich_go) > 1) {
      ## acyclic_graph
      acyclic_plot <- goplot(enrich_go)
      enrichment_plot[[paste0("GO_", go)]] <- (bar_plot + dot_plot)/acyclic_plot
    }else{
      enrichment_plot[[paste0("GO_", go)]] <- (bar_plot + dot_plot)
    }
  }
  
  # KEGG
  enrich_kegg <- enrichKEGG(gene = genes_entrezid$ENTREZID,
                   organism = 'hsa', #KEGG可以用organism = 'hsa'
                   pvalueCutoff = 0.05)
  KEGG_bar_plot <- barplot(enrich_kegg, title="Enrichment KEGG bar")
  KEGG_dot_plot <- dotplot(enrich_kegg, title="Enrichment KEGG dot")
  enrichment_plot[["kegg"]] <- KEGG_bar_plot + KEGG_dot_plot
  
  # Reactome
  enrich_reactome <- enrichPathway(gene = genes_entrezid$ENTREZID,
                                   pvalueCutoff = 0.05,
                                   readable = TRUE)
  reactome_bar_plot <- barplot(enrich_reactome,title="Enrichment Reactome bar")
  reactome_dot_plot <- dotplot(enrich_reactome,title="Enrichment Reactome dot")
  enrichment_plot[["reactome"]] <- reactome_bar_plot + reactome_dot_plot
  
  return(enrichment_plot)
}
