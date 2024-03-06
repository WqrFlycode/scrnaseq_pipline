Enrich <- function(genes, orgdb = "org.Hs.eg.db"){
  #将SYMBOL格式转为ENSEMBL和ENTERZID格式 
  genes_entrezid <- clusterProfiler::bitr(
    genes, 
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = orgdb
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
  }
  
  # KEGG
  enrich_kegg <- enrichKEGG(gene = genes_entrezid$ENTREZID,
                            organism = 'hsa', #KEGG可以用organism = 'hsa'
                            pvalueCutoff = 0.05)
  
  # Reactome
  enrich_reactome <- enrichPathway(gene = genes_entrezid$ENTREZID,
                                   pvalueCutoff = 0.05,
                                   readable = TRUE)
  
  cat("--------------------Enrich finished--------------------")
  enrichment <- list(
    enrich_go = enrich_go, 
    enrich_kegg = enrich_kegg, 
    enrich_reactome = enrich_reactome
  )
  return(enrichment)
}
