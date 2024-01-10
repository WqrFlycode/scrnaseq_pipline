

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
  
  print("----------Read data finished----------")
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
