FindDEG <- function(seurat_object, metanames = NULL, ncl = 1){
  cl <- makeCluster(ncl)
  all.markers <- list()
  for (metaname in metanames) {
    cat("\n %%%%% run FindMarkers by", metaname, "%%%%% \n")
    if(!class(seurat_object@meta.data[[metaname]]) == "factor") {
      seurat_object@meta.data[[metaname]] <- as.factor(
        seurat_object@meta.data[[metaname]]
      )
    }
    Seurat::Idents(seurat_object) <- seurat_object@meta.data[[metaname]]
    clusterName <- levels(seurat_object@meta.data[[metaname]])
    names(clusterName) <- clusterName
    clusterExport(
      cl,
      c("clusterName", "seurat_object"),
      envir = environment()
    )
    deg_results <- pbapply::pblapply(
      X = clusterName,
      FUN = function(x) {
        clmarker <- Seurat::FindMarkers(
          object = seurat_object,
          ident.1 = x
        )
        clmarker <- subset(clmarker, p_val < 0.01)
        return(clmarker)
      },
      cl = cl
    )
    for(i in seq(length(deg_results))) {
      deg_results[[i]] <- data.frame(
        deg_results[[i]],
        cluster = as.factor(
          rep(
            names(deg_results)[i],
            nrow(deg_results[[i]])
          )
        )
        ,
        gene = rownames(deg_results[[i]])
      )
    }
    deg_results <- do.call(rbind,deg_results)
    rownames(deg_results) <- NULL
    all.markers[[metaname]] <- deg_results
  }
  stopCluster(cl)
  
  cat("----------DifferentialExpressing finished----------")
  return(all.markers)
}
