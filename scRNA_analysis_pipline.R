source("scripts/library_source.R")

# read data
Data <- ReadData_10X(
  data_dir = "E:/Scripts/Rproject/scRNA_analysis/case/",
  filename = "case", # 10X的三个文件名的前缀
  data_name = "case",
  Species = "human", # "human"; "mouse"; "zebrafish"
  results_dir = NULL
)

# analysis
Data <- analysis_scrnaseq(
  Data,
  min_nFeature = 100, min_nCount = 100, 
  max_percent_mito = 30, max_percent_ribo = 50,
  pc_num = 50,resolution = 0.5,
  ref_singler_dir = "E:/Scripts/Rproject/scRNA_analysis/celldex_ref_dataset/",
  ref_cell_dex = "HumanPrimaryCellAtlasData",
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

# 读取数据-----
## data_name ：数据名称，默认为"case"，名称不宜过长
## data_dir ：单细胞数据文件路径
## data_format: "10x","rds","txt","h5","csv"
## results_dir ：输出结果的路径，默认为NULL, 将结果保存在数据文件夹data_name_results中
## 注意：
## 1.第一次读取数据后seurat对象被保存在rds文件中，data_dir和results_dir被保存在seurat对象中
## 2.所有路径中使用“/”，最后一级路径后面也要带“/”
## 3.文件名称应为*matrix.mtx.gz; *barcodes.tsv.gz; *

# data_dir <- "D:/Data/scRNA-seq/GSE227425_RAW/GSM710024_snAB1/"
# Data <- ReadData(
#   data_name = "GSM710024_snAB1", 
#   data_dir = data_dir, 
#   results_dir = NULL
# )

## read csv
data_dir <- "E:/Data/单细胞数据-王章/人的4个样本-GSE130646_RAW/GSM3746212_Muscle_1_Counts.csv/"
Data <- ReadData_csv(
  data_name = "GSM3746212_Muscle_1",
  data_dir = data_dir,
  results_dir = NULL
)

## 整合样本
paths <- c(
  "E:/Data/单细胞数据-王章/人的4个样本-GSE130646_RAW/GSM3746212_Muscle_1_Counts.csv/GSM3746212_Muscle_1_results/GSM3746212_Muscle_1.rds",
  "E:/Data/单细胞数据-王章/人的4个样本-GSE130646_RAW/GSM3746213_Muscle_2_Counts.csv/GSM3746213_Muscle_2_results/GSM3746213_Muscle_2.rds",
  "E:/Data/单细胞数据-王章/人的4个样本-GSE130646_RAW/GSM3746214_Muscle_3_Counts.csv/GSM3746214_Muscle_3_results/GSM3746214_Muscle_3.rds",
  "E:/Data/单细胞数据-王章/人的4个样本-GSE130646_RAW/GSM3746215_Muscle_4_Counts.csv/GSM3746215_Muscle_4_results/GSM3746215_Muscle_4.rds"
)
paths <- c(
  "E:/Data/单细胞数据-王章/老鼠的2个样本-GSE110878_RAW/GSM3520458_20171018_uninjured_wt_filtered.csv/GSM3520458_results/GSM3520458.rds",
  "E:/Data/单细胞数据-王章/老鼠的2个样本-GSE110878_RAW/GSM3520459_20180917_uninjured_wt_filtered.csv/GSM3520459_results/GSM3520459.rds"
)
Data <- integration(
  seurat_objects_path = paths,
  results_dir = "E:/Data/单细胞数据-王章/老鼠的2个样本-GSE110878_RAW/",
  data_name = "GSE110878"
)
Data <- readRDS("E:/Data/单细胞数据-王章/老鼠的2个样本-GSE110878_RAW/GSE110878_results/GSE110878.rds")

# 质量控制-----
## Data：seurat对象
Data <- QualityControl(Data = Data)

# 降维-----
## pc_num：降维后保留的主成分数量，默认为50
Data <- DimensionalityReduction(Data = Data, pc_num = 50)

# 聚类-----
## resolution：该参数影响分类数目，数值越大分类数目越多，最大值为1
## re：逻辑值，TRUE重新分析，FALSE使用上次分析结果作图

Data <- Clustering(Data = Data, resolution = 1, re = TRUE)

# 细胞注释-----
## ref_dir：SingleR参考数据集路径，"scRNA_analysis/celldex_ref_dataset/"
## ref_markers：参考标记基因列表
# markers_list <- list(
#   Adipocyte = c("ADIPOQ", "PPARG", "CD36"),
#   Endothelial_cells = c("MMRN1", "FLT4", "PARD6G", "JAM2", "ENG", "PLVAP", "GPNM8", "DAB2", "CDSD"),
#   Epithelial_cells = c("CSN2", "PGR", "ESR1", "PRLR", "LALBA", "CSN3", "EPCAM", "KRT18", "ELF5", "TPM2", "CLDN4"),
#   # Fibroblasts
#   Immune_cells = c("CD136", "PTPRC", "ZEB2", "BANK1", "CST3", "COX2", "CD3E", "KLRK1", "ATP6"),
#   Myoepithelial_cells = c("MYH11", "ACAT2", "MYLK", "COL3A1", "COL1A2", "COL1A1", "DCN", "POSTN", "COL6A3"),
#   Precursor_cells = c("IL1R1", "DLG2", "LARGE1", "NAV2", "ZNF710", "SERINC5", "RABEP1", "MSI2", "WWOX", "LIFR", "KIF15", "TPX2", "SMC4", "TOP2A")
# )
"celldex_ref_dataset/BlueprintEncodeData.rds"
"celldex_ref_dataset/MouseRNAseqData.rds"
Data <- Annotation(
  Data = Data, 
  ref_singler_dir = "celldex_ref_dataset/",ref_cell_dex = "MouseRNAseqData",
  # ref_markers = markers_list,
  re = FALSE
)

# 差异基因-----

Data <- DifferentialExpressing(Data = Data, re = FALSE)

# 富集分析-----
## results_dir：同上，用来读取DifferentialExpressing中的data_name_all_markers.csv
Enrich(data_name = Data@tools$data_name, results_dir = Data@tools$results_dir)

saveRDS(Data, paste0(Data@tools$results_dir,Data@tools$data_name,"_result.rds"))
# rm(Data,data_dir)

# 轨迹分析-----
## by_cluster：细胞类型或分组
Data <- TrajectoryAnalysis(Data = Data, by_cluster = Data$cluster_MouseRNAseqData)

# 细胞通讯-----
Data <- CellCommunication(Data = Data, by_cluster = names(Data@meta.data)[7], re = FALSE)

scConf = createConfig(Data)
makeShinyApp(
  Data, 
  scConf, 
  gene.mapping = TRUE,
  shiny.title = "ShinyCell Quick Start",
  shiny.dir = paste0(data_dir, "shinyApp/")
) 
