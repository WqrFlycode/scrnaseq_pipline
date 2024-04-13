QualityControl <- function(
    Data,
    min_nFeature = 100, min_nCount = 100,
    max_percent_mito = 30, max_percent_ribo = 50
){
  info <- Data@tools$info
  # create QC dir
  info$dir$qc <- paste0(info$data_name, "_qc/")
  QC_dir <- paste0(info$dir$dir, info$dir$qc)
  if (!dir.exists(QC_dir)) dir.create(QC_dir)
  
  # origin info
  # compute the proportion of mito-genes expression
  if (info$Species == "mouse") {
    pattern_mt <- "^mt-"
  }else{
    pattern_mt <- "^MT-"
  }
  Data[["percent.mito"]] <- PercentageFeatureSet(Data, pattern = pattern_mt)
  # compute the proportion of ribo-genes expression
  Data[["percent.ribo"]] <- PercentageFeatureSet(Data, pattern = "^RP[SL]")
  # 计算红血细胞基因比例
  # if(! "percent.hb" %in% names(Data@meta.data)){
  #   Data[["percent.hb"]] <- PercentageFeatureSet(Data, pattern = "^Hb[^(p)]")
  # }
  # percentage times 100
  
  plot_qc(Data, output_dir = QC_dir, status = "raw")
  
  # QC-----
  cat("\n %%%%% run quality control %%%%% \n") # 3sigma, mad
  Data <- subset(
    Data,
    subset =
      nFeature_RNA > min_nFeature &
      nCount_RNA > min_nCount &
      percent.mito < max_percent_mito &
      percent.ribo < max_percent_ribo,
    features = row.names(Data)[rowSums(GetAssayData(Data,slot = "counts")>0)>0]
  )
  # box
  outlier_range <- function(x, n = 1.5){
    l <- sum(fivenum(x)*c(0,1+n,0,-n,0))
    u <- sum(fivenum(x)*c(0,-n,0,1+n,0))
    return(c(l,u))
  }
  max_nFeature <- outlier_range(Data$nFeature_RNA,n = 3)[2] # lower & upper
  max_nCount <- outlier_range(Data$nCount_RNA,n = 3)[2]
  Data <- subset(
    Data,
    subset =
      nFeature_RNA < max_nFeature &
      nCount_RNA < max_nCount
  )
  info$qc_threshold <- c(
    min_nFeature, max_nFeature,
    min_nCount, max_nCount,
    max_percent_mito, max_percent_ribo
  )
  names(info$qc_threshold) <- c(
    "min_nFeature", "max_nFeature",
    "min_nCount", "max_nCount",
    "max_percent_mito", "max_percent_ribo"
  )
  info$dim$filter <- dim(Data)
  names(info$dim$filter) <- c("gene", "cell")
  data_dim <- rbind(info$dim$raw, info$dim$filter, info$dim$raw-info$dim$filter)
  rownames(data_dim) <- c("raw","filter","remove")
  cat("\n");print(data_dim)
  cat("\nremove percent\n")
  print(  round(data_dim[3,]/data_dim[1,]*100, digits = 2))
  
  plot_qc(Data, output_dir = QC_dir, status = "qc")
  
  
  # save qc data
  info$filename$qc <- paste0(info$data_name, "_qc_seurat.rds")
  Data@tools$info <- info
  qc_data_path <- paste0(info$dir$dir,info$dir$seurat,info$filename$qc)
  saveRDS(Data, qc_data_path)
  cat("\nsave qc data to: \n", qc_data_path)
  
  # save info
  info_path <- paste0(info$dir$dir,info$dir$seurat,info$filename$info)
  saveRDS(info, info_path)
  cat("\nsave info to: \n", info_path, "\n")
  
  cat("\n----------Quality control finished----------\n")
  return(Data)
}
