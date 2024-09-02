QualityControl <- function(
    Data, qc_times = 3, mincells = 1,
    min_nFeature = NULL, max_nFeature = NULL,
    min_nCount = NULL, max_nCount = NULL,
    max_percent_mito = NULL, max_percent_ribo = NULL
){
  info <- Data@tools$info
  # create QC dir
  info$dir$qc <- paste0(info$data_name, "_qc/")
  QC_dir <- paste0(info$dir$dir, info$dir$qc)
  if (!dir.exists(QC_dir)) dir.create(QC_dir)
  
  # origin info
  # compute the proportion of mito-genes expression
  mtgenes <- c(
    # "^MT-"
    "MT-ATP6","MT-ATP8",
    "MT-CO1","MT-CO2","MT-CO3",
    "MT-CYB","MT-CYTB",
    "MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
    
    # "^MT"
    "MTATP6","MTATP8",
    "MTCO1","MTCO2","MTCO3",
    "MTCYB","MTCYTB",
    "MTND1","MTND2","MTND3","MTND4","MTND4L","MTND5","MTND6",
    
    # without "MT"
    "ATP6","ATP8",
    "CO1","CO2","COII","Cu","CO3","COX1","COX2","COX3",
    "CYB","CYTB","cytb",
    "ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
    
    # uppercase letters
    "mt-Atp6","mt-Atp8",
    "mt-Co1","mt-Co2","mt-Co3",
    "mt-Cyb","mt-Cytb",
    "mt-Nd1","mt-Nd2","mt-Nd3","mt-Nd4","mt-Nd4l","mt-Nd5","mt-Nd6",
    
    # losercase
    "mt-atp6","mt-atp8",
    "mt-co1","mt-co2","mt-co3",
    "mt-cyb","mt-cytb",
    "mt-nd1","mt-nd2","mt-nd3","mt-nd4","mt-nd4l","mt-nd5","mt-nd6"
  )
  included_mtgenes <- mtgenes[mtgenes %in% rownames(Data)]
  cat("mito genes:\n");print(included_mtgenes)
  Data[["percent.mito"]] <- PercentageFeatureSet(Data, features = included_mtgenes)
  # compute the proportion of ribo-genes expression
  Data[["percent.ribo"]] <- PercentageFeatureSet(Data, pattern = "^RP[SL]")
  # 计算红血细胞基因比例
  # if(! "percent.hb" %in% names(Data@meta.data)){
  #   Data[["percent.hb"]] <- PercentageFeatureSet(Data, pattern = "^Hb[^(p)]")
  # }
  # percentage times 100
  
  plot_qc(Data, output_dir = QC_dir, status = "raw")
  
  # QC-----
  if(is.null(min_nFeature)) min_nFeature = 100
  if(is.null(min_nCount)) min_nCount = 100
  if(is.null(max_percent_mito)) max_percent_mito = 30
  if(is.null(max_percent_ribo)) max_percent_ribo = 50
  
  cat("\n %%%%% run quality control %%%%% \n") # 3sigma, mad
  Data <- subset(
    Data,
    subset =
      nFeature_RNA > min_nFeature &
      nCount_RNA > min_nCount &
      percent.mito < max_percent_mito &
      percent.ribo < max_percent_ribo,
    features = row.names(Data)[
      rowSums(GetAssayData(Data, layer= "counts") > 0) >= mincells
    ]
  )
  # box
  outlier_range <- function(x, n = 1.5){
    l <- sum(fivenum(x)*c(0,1+n,0,-n,0))
    u <- sum(fivenum(x)*c(0,-n,0,1+n,0))
    return(c(l,u))
  }
  if(is.null(max_nFeature)) {
    max_nFeature <- outlier_range(Data$nFeature_RNA,n = qc_times)[2] # lower & upper
  }
  if(is.null(max_nCount)) max_nCount <- outlier_range(Data$nCount_RNA,n = qc_times)[2]
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
  
  seurat_dir <- paste0(info$dir$dir,info$dir$seurat)
  if(!dir.exists(seurat_dir)) dir.create(seurat_dir)
  # save qc data
  info$filename$qc <- paste0(info$data_name, "_qc_seurat.rds")
  Data@tools$info <- info
  qc_data_path <- paste0(seurat_dir,info$filename$qc)
  saveRDS(Data, qc_data_path)
  cat("\nsave qc data to: \n", qc_data_path)
  
  # save info
  info_path <- paste0(seurat_dir,info$filename$info)
  saveRDS(info, info_path)
  cat("\nsave info to: \n", info_path, "\n")
  
  cat("\n----------Quality control finished----------\n")
  return(Data)
}
