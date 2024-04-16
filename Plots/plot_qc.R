plot_qc <- function(Data, output_dir, status) {
  cat("\n %%%%% plot QC index %%%%% \n")
  
  results_dir <- paste0(output_dir, Data@tools$info$data_name, "_", status, "_")
  
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
  setbreak <- function(Range) {
    maxint <- signif(Range[2], digits = 1)
    intbreaks <- seq.int(0, maxint, length.out = 3)[-1]
    breaks <- c(Range[1], intbreaks, Range[2])
    breaks <- sort(breaks)
    breaks <- breaks[!duplicated(signif(breaks,digits = 1))]
    return(breaks)
  }
  plot_qc[[1]] <- plot_qc[[1]] + 
    scale_y_continuous(breaks = setbreak(range(Data$nFeature_RNA)))
  plot_qc[[2]] <- plot_qc[[2]] + 
    scale_y_continuous(breaks = setbreak(range(Data$nCount_RNA)))
  ggsave(
    paste0(results_dir,"violin.png"),
    plot_qc,
    width = 5+3*length(unique(Data$orig.ident)),
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
  ggsave(
    paste0(results_dir,"scatter.png"),
    plot_FeatureScatter,
    width = 12, height = 6
  )
  rm(plot_FeatureScatter, plot_qc, qc_index)
}
