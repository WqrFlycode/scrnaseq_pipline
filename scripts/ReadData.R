saveinfo <- function(data_name, data_dir, Species, rawdim){
  info <- list()
  info$data_name <- data_name
  info$Species <- Species
  info$dir$dir <- paste0(dirname(data_dir), "/")
  info$dir$data <- paste0(basename(data_dir), "/")
  info$dir$seurat <- paste0(basename(dirname(data_dir)),"_seurat/")
  info$filename$info <- paste0(data_name,"_info.rds")
  info$filename$raw <- paste0(data_name, "_raw_seurat.rds")
  info$dim$raw <- rawdim
  names(info$dim$raw) <- c("gene", "cell")
  cat("\n----------create", data_name, "info----------")
  return(info)
}

saveSeuratData <- function(Data) {
  info <- Data@tools$info
  seurat_dir <- paste0(info$dir$dir, info$dir$seurat)
  if(!dir.exists(seurat_dir)) dir.create(seurat_dir)
  raw_data_path <- paste0(seurat_dir, info$filename$raw)
  saveRDS(Data, raw_data_path)
  cat("\nsave raw seurat to: \n", raw_data_path, "\n")
}

ReadData_10X <- function(data_dir, filename, data_name = "case", Species = NULL){
  files <- list.files(data_dir)
  exist_files <- files[grep(filename, files)]
  exist_data_names <- rep(NA, 3)
  names(exist_data_names) <- c("matrix", "barcodes", "features")
  if (length(grep("matrix.mtx", exist_files)) == 1) {
    exist_data_names["matrix"] <- exist_files[grep("matrix.mtx", exist_files)]
  }
  if (length(grep("barcodes.tsv", exist_files)) == 1) {
    exist_data_names["barcodes"] <- exist_files[grep("barcodes.tsv", exist_files)]
  }
  if (length(grep("features.tsv", exist_files)) == 1) {
    exist_data_names["features"] <- exist_files[grep("features.tsv", exist_files)]
  }
  if (length(grep("genes.tsv", exist_files)) == 1) {
    exist_data_names["features"] <- exist_files[grep("genes.tsv", exist_files)]
  }
  
  # 判断三个文件是否存在，并读取
  if (!any(is.na(exist_data_names))) {
    sample_files <- paste0(data_dir,exist_data_names)
    # read expression matrix
    Data <- Matrix::readMM(file = sample_files[1])
    
    # read cells
    barcode.names = read.delim(
      sample_files[2],
      header = FALSE,
      stringsAsFactors = FALSE
    )
    colnames(Data) = barcode.names$V1
    
    # read genes
    gene.names = read.delim(
      sample_files[3],
      header = FALSE, 
      stringsAsFactors = FALSE
    )
    if (ncol(gene.names) == 1) {
      rownames(Data) <-  gene.names$V1  
    } else {
      rownames(Data) <-  gene.names$V2
    }
    
    # create Seurat object
    Data <- CreateSeuratObject(counts = Data, project = data_name)
    
    # save parameters to info
    info <- saveinfo(
      data_name = data_name,
      data_dir = data_dir,
      Species = Species,
      rawdim = dim(Data)
    )
    Data@tools$info <- info
    saveSeuratData(Data)
    cat("----------Read 10X data", info$data_name, "finished----------\n")
    return(Data)
  }else{
    stop("Missing 10X files")
  }
}

ReadData_h5 <- function(data_dir, filename_prefix, data_name = "case", Species = NULL) {
  files <- list.files(data_dir)
  filename <- files[grep(filename_prefix, files)]
  Data <- Read10X_h5(paste0(data_dir, filename))
  Data <- CreateSeuratObject(Data)
  # save info
  info <- saveinfo(
    data_name = data_name,
    data_dir = data_dir,
    Species = Species,
    rawdim = dim(Data)
  )
  Data@tools$info <- info
  saveSeuratData(Data)
  cat("----------Read h5 data", info$data_name, "finished----------\n")
  return(Data)
}

ReadData_txt <- function(data_dir, filename_prefix, data_name = "case", Species = NULL, trans = FALSE) {
  files <- list.files(data_dir)
  filename <- files[grep(filename_prefix, files)]
  Data <- data.table::fread(paste0(data_dir, filename))
  print(Data[1:3,1:3]);print(dim(Data))
  
  if(trans) {
    barcodes <- Data$V1
    genes <- names(Data)[-1]
    Data <- Data[,-1]
    Data <- t(as.matrix(Data))
  } else {
    barcodes <- names(Data)[-1]
    genes <- Data$V1
    Data <- Data[,-1]
    Data <- as.matrix(Data)
  }
  
  Data <- as(Data, "CsparseMatrix")
  colnames(Data) <- barcodes
  rownames(Data) <- genes
  print(Data[1:3,1:3])
  
  Data <- CreateSeuratObject(counts = Data)
  # save info
  info <- saveinfo(
    data_name = data_name,
    data_dir = data_dir,
    Species = Species,
    rawdim = dim(Data)
  )
  Data@tools$info <- info
  saveSeuratData(Data)
  cat("----------Read count", info$data_name, "finished----------\n")
  return(Data)
}

ReadData_csv <- function(data_dir, data_name = "case", results_dir = NULL){
  if (is.null(results_dir)) {
    results_dir <- paste0(data_dir,"/",data_name,"_results/")
  }
  # create result direction
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)
  }
  cat("## results_dir: \n",results_dir, "\n")
  
  if (file.exists(paste0(results_dir,data_name,".rds"))) {
    print("Seurat object existed")
    data <- readRDS(paste0(results_dir,data_name,".rds"))
  }else{
    file_name <- list.files(data_dir,pattern = "*.csv")
    cat("## Read data: \n", file_name, "\n")
    data <- read.csv(file = paste0(data_dir, file_name))
    gene_names <- data[,1]
    data <- as.matrix(data[,-1])
    rownames(data) <- gene_names
    data <- CreateSeuratObject(counts = data)
  }
  
  # save path
  data@tools$results_dir <- results_dir
  data@tools$data_name <- data_name
  
  # save rds
  saveRDS(data,file = paste0(results_dir,"/",data_name,".rds"))
  
  print("----------Read data finished----------\n")
  return(data)
}

ReadData_merge_10X <- function(data_name = "case", data_dir, results_dir = NULL){
  if (is.null(results_dir)) {
    results_dir <- paste0(data_dir,"/",data_name,"_results/")
  }
  # create result direction
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)
  }
  if (file.exists(paste0(results_dir,data_name,".rds"))) {
    data <- readRDS(paste0(results_dir,data_name,".rds"))
  }else{
    files <- list.files(data_dir)
    # 文件夹中的所有样本
    samples <- sapply(
      files[grep("matrix.mtx",files)], 
      function(x) {
        substr(
          x = x,
          start = 1,
          stop = regexpr("matrix.mtx.gz",x)-2
        )
      },
      USE.NAMES = FALSE
    )
    # 创建Seurat对象list
    for (i in 1:length(samples)) {
      exist_files <- files[grep(samples[i], files)]
      list(
        matrix = exist_files[grep("matrix.mtx.gz", exist_files)],
        barcodes = exist_files[grep("barcodes.tsv.gz", exist_files)],
        genes = exist_files[grep("features.tsv.gz", exist_files)]
      )
      is_exist <- c(
        length(grep("matrix.mtx.gz", exist_files)) == 1,
        length(grep("barcodes.tsv.gz", exist_files)) == 1,
        length(grep("features.tsv.gz", exist_files)) == 1 | 
          length(grep("genes.tsv.gz", exist_files)) == 1
      )
      # 判断三个文件是否存在，并读取
      if (all(is_exist)) {
        sample_files <- exist_files[is_exist]
        sample_files <- paste0(data_dir,sample_files)
        # read expression matrix
        data <- Matrix::readMM(file = sample_files[1])
        # read cells
        barcode.names = read.delim(
          sample_files[2],
          header = FALSE,
          stringsAsFactors = FALSE
        )
        colnames(data) = barcode.names$V1
        # read genes
        gene.names = read.delim(
          sample_files[3],
          header = FALSE, 
          stringsAsFactors = FALSE
        )
        if (ncol(gene.names) == 1) {
          rownames(data) <-  gene.names$V1  
        } else {
          rownames(data) <-  gene.names$V2
        }
        # create Seurat object
        data <- CreateSeuratObject(counts = data, project = samples[i])
        print(paste0("Read data ", samples[i]))
      } 
    }
    
    if (length(data) == 1) {
      data <- data[[1]]
    } else {
      data_name_list <- sapply(data, function(x) x@project.name)
      # integrated all samples
      data <- merge(
        x = data[[1]], 
        y = data[-1], 
        add.cell.ids = data_name_list, 
        project = data_name
      )
    }
    
    # save path
    data@tools$results_dir <- results_dir
    data@tools$data_name <- data_name
    
    # save rds
    saveRDS(data,file = paste0(results_dir,"/",data_name,".rds"))
  }
  
  print("----------Read data finished----------")
  return(data)
}

Read_results <- function(info_dir) {
  info <- readRDS(info_dir)
  Data <- readRDS(paste0(info$results_dir, info$data_name,"_results.rds"))
  Data@tools$all.markers <- readRDS(info$all.markers_dir)
  Data@tools$cds <- readRDS(info$cds_dir)
  if("cellchat_dir" %in% names(info)){
    Data@tools$cellchat <- readRDS(info$cellchat_dir)
  }
  return(Data)
}