source("scripts/library_source.R")

# read data
## data_dir:    单细胞数据文件路径
## filename:    10X的三个文件名的前缀
## data_name:   数据名称，默认为"case"，名称不宜过长
## results_dir: 输出结果的路径，默认NULL, 结果保存在数据文件夹中
## 注意：
## 1.第一次读取数据后seurat对象被保存在rds文件中，data_dir和results_dir被保存在seurat对象中
## 2.所有路径中使用“/”，最后一级路径后面也要带“/”
## 3.文件名称应为*matrix.mtx.gz; *barcodes.tsv.gz; *

Data <- ReadData_10X(
  data_dir = "D:/Data/scRNA-seq/mouse/GSE226169_RAW/",
  filename = "GSM7066259",
  data_name = "GSM7066259",
  Species = "mouse", # "human"; "mouse"; "zebrafish"
  results_dir = NULL
)

# analysis
Data <- analysis_scrnaseq(
  Data,
  min_nFeature = 100, min_nCount = 100, 
  max_percent_mito = 30, max_percent_ribo = 50,
  pc_num = 50,resolution = 0.5,
  ref_singler_dir = "D:/Workspace/R/scRNA/celldex_ref_dataset/",
  ref_cell_dex = "ImmGenData",
  ref_markers = NULL
)

Plot_seurat(Data = Data, figure_format = "png")

info <- Data@tools$info
rmarkdown::render(
  input = "scripts/report.Rmd",
  params = list(info = info),
  output_file = paste0(
    info$results_dir, 
    info$data_name,
    "_report"
  )
)

# 并行
library(future)
plan(sequential)
plan(multisession, workers = 2)
plan()

# ShinyApp
scConf = createConfig(Data)
makeShinyApp(
  Data, 
  scConf, 
  gene.mapping = TRUE,
  shiny.title = "ShinyCell Quick Start",
  shiny.dir = paste0(data_dir, "shinyApp/")
) 
