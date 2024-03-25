# marker_list: 列表，每个类的差异基因
Enrich <- function(marker_list, Species, ncl){
  if(Species == "human") {
    orgdb <- "org.Hs.eg.db"
  } else if(Species == "mouse") {
    orgdb <- "org.Mm.eg.db"
  } else if(Species == "pig") {
    orgdb <- "org.Ss.eg.db"
  } else if(Species == "chicken") {
    orgdb <- "org.Gg.eg.db"
  }
  
  #将SYMBOL格式转为ENSEMBL和ENTERZID格式 
  entrezid_ls <- lapply(marker_list, function(x){
    clusterProfiler::bitr(
      x, 
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = orgdb
    )
  })
  
  if(Species == "human") {
    Organisms <- "hsa" # Homo sapiens (human)
  } else if(Species == "mouse") {
    Organisms <- "mmu" # Mus musculus (house mouse)
  } else if(Species == "pig") {
    Organisms <- "org.Ss.eg.db"
  }
  
  cl <- makeCluster(ncl)
  clusterExport(
    cl,
    c("entrezid_ls", "orgdb", "Organisms", "Species"),
    envir = environment()
  )
  enrichment <- list()
  
  cat("\n %%%%% enrich GO %%%%% \n")
  GO <- c("MF","BP","CC")
  for (go in GO) {
    clusterExport(cl, "go", envir = environment())
    enrichment[[go]] <- parLapply(
      cl,
      entrezid_ls,
      fun = function(x){
        enrich_go <- clusterProfiler::enrichGO(
          gene = x$ENTREZID,
          OrgDb = orgdb,
          keyType = "ENTREZID",
          ont = go,
          readable = TRUE
        )
        # dplyr::filter(enrich_go@result, p.adjust < 0.05)
        enrich_go@geneSets <- list()
        enrich_go@universe <- "NULL"
        return(enrich_go)
      }
    )
  }
  
  cat("\n %%%%% enrich KEGG %%%%% \n")
  enrichment[["kegg"]] <- parLapply(
    cl,
    entrezid_ls,
    fun = function(x){
      tryCatch({
        enrich_kegg <- clusterProfiler::enrichKEGG(
          gene = x$ENTREZID,
          organism = Organisms,
          pvalueCutoff = 0.05
        )
        # dplyr::filter(enrich_kegg@result, p.adjust < 0.05)
        enrich_kegg@geneSets <- list()
        enrich_kegg@universe <- "NULL"
        return(enrich_kegg)
      }, error = function(e){
        message("%%% Failed to download KEGG data %%%")
        list()
      })
    }
  )
  
  cat("\n %%%%% enrich Reactome %%%%% \n")
  enrichment[["Reactome"]] <- parLapply(
    cl,
    entrezid_ls,
    fun = function(x){
      tryCatch({
        enrich_reactome <- ReactomePA::enrichPathway(
          gene = x$ENTREZID,
          pvalueCutoff = 0.05,
          organism = Species,
          readable = TRUE
        )
        # dplyr::filter(enrich_reactome@result, p.adjust < 0.05)
        enrich_reactome@geneSets <- list()
        enrich_reactome@universe <- "NULL"
        return(enrich_reactome)
      }, error = function(e){
        message("%%% Failed to enrich reactome %%%")
        list()
      })
    }
  )
  stopCluster(cl)
  names(enrichment)[match(GO,names(enrichment))] <- c("GO_MF", "GO_BP", "GO_CC")
  
  cat("--------------------Enrich finished--------------------\n")
  return(enrichment)
}
