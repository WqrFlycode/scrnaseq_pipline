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
    barcodes <- Data[[1]]
    genes <- names(Data)[-1]
    Data <- Data[,-1]
    Data <- t(as.matrix(Data))
  } else {
    barcodes <- names(Data)[-1]
    genes <- Data[[1]]
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

# type: SYMBOL, ENSEMBL
transferGene <- function(genes, fromtype, totype, ref) {
  index <- which(genes %in% ref[[fromtype]]) # index of fromtype genes
  fromtype_genes <- genes[index]
  totype_genes <- ref[[totype]][match(fromtype_genes,ref[[fromtype]])]
  genes[index] <- totype_genes
  return(genes)
}
