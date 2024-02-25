merge_data <- function(seurat_objects_path, results_dir, data_name) {
  data <- lapply(
    seurat_objects_path,
    function(x){
      data <- readRDS(x)
      data$orig.ident <- data@tools$data_name
      return(data)
    }
  )
  # data_name_list <- sapply(data, function(x) x@tools$data_name)
  
  # integrated all samples
  data <- merge(
    x = data[[1]], 
    y = data[-1], 
    # add.cell.ids = data_name_list, 
    project = data_name
  )
  
  # save path
  if (!dir.exists(results_dir)) dir.create(results_dir)
  results_dir <- paste0(results_dir, "/", data_name, "_results/")
  if(!dir.exists(results_dir)) dir.create(results_dir)
  data@tools$results_dir <- results_dir
  data@tools$data_name <- data_name
  
  # save rds
  saveRDS(data,file = paste0(results_dir,"/",data_name,".rds"))
  return(data)
}


#   
# dirs <- list.dirs("D:/Data/scRNA-seq/GSE227425_RAW/")[-1]
# dirlist <- strsplit(dirs,split = "/")
# data_list <- sapply(1:10, function(x) {dirlist[[x]][5]})
# 
# for (i in 1:10) {
#   Data <- ReadData(
#     data_name = data_list[i], 
#     data_dir = paste0("D:/Data/scRNA-seq/GSE227425_RAW/", data_list[i],"/"), 
#     results_dir = NULL
#   )
#   
#   alldata <- lapply(1:10, function(i) {
#     readRDS(
#       paste0(
#         "D:/Data/scRNA-seq/GSE227425_RAW/", 
#         data_list[i],"/",
#         data_list[i],"_results/",
#         data_list[i],".rds"
#       )
#     )
#   })
#   pig <- merge(
#     x = alldata[[1]], 
#     y = alldata[-1], 
#     add.cell.ids = data_list, 
#     project = "pigall"
#   )
#   
#   pig12 <- merge(pig[[1]],pig[[2]])
#   Data$orig.ident
#   
#   
#   pig12 <- lapply(X = pig12, FUN = function(x) {
#     x <- NormalizeData(x)
#     x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#   })
#   
#   features <- SelectIntegrationFeatures(object.list = pig12)
#   
#   Data <- FindIntegrationAnchors(object.list = pig12)
#   Data <- IntegrateData(anchorset = Data)
#   
# }
# 
# Datalist <- readRDS("D:/Data/scRNA-seq/GSE227425_RAW/GSE227425_RAW.rds")
# 
# Datalist <- list(Data1,Data2)
# Datalist <- Datalist[1:2]
# 
# # normalize and identify variable features for each dataset independently
# Datalist <- lapply(X = Datalist, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# 
# # select features that are repeatedly variable across datasets for integration
# features <- SelectIntegrationFeatures(object.list = Datalist)
# 
# Data.anchors <- FindIntegrationAnchors(object.list = Datalist, anchor.features = features)
# 
# # this command creates an 'integrated' data assay
# Data <- IntegrateData(anchorset = Data.anchors)
# 
# DefaultAssay(Data) <- "integrated"
# 
# Data <- ScaleData(Data, verbose = FALSE)
# Data <- RunPCA(Data, npcs = 30, verbose = FALSE)
# Data <- RunUMAP(Data, reduction = "pca", dims = 1:30)
# Data <- FindNeighbors(Data, reduction = "pca", dims = 1:30)
# Data <- FindClusters(Data, resolution = 0.5)
# 
# # Visualization
# p1 <- DimPlot(Data, reduction = "umap", group.by = "orig.ident")
# p2 <- DimPlot(Data, reduction = "umap", label = TRUE, repel = TRUE)
# p1 + p2
# 
# DimPlot(Data, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident")
# DimPlot()
