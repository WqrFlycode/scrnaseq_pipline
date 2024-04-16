rm(list = ls())
source("scripts/library_source.R")

# info <- readRDS("D:/Data/scRNA-seq/case/case_info.rds")
datadir <- "D:/Data/scRNA-seq/case/case_origin/"
Data <- ReadData_10X(
  data_dir = datadir,
  filename = "case",
  data_name = "case",
  Species = "human" # "human"; "mouse"; "zebrafish"
)
info <- Data@tools$info

datadir <- "D:/Data/scRNA-seq/G_E_00012/G_E_00012_data/"
Data <- ReadData_h5(
  data_dir = datadir,
  filename_prefix = "GSM4502481_chicken_heart_scRNAseq_D14_RV",
  data_name = "D14_RV",
  Species = "chicken"
)

Data <- readRDS(paste0(info$dir$dir, info$filename$raw))
Data <- readRDS("D:/Data/scRNA-seq/case/case_seurat/case_qc_seurat.rds")
Data <- QualityControl(
  Data,
  min_nFeature = 100, min_nCount = 100,
  max_percent_mito = 30, max_percent_ribo = 50
)
info <- Data@tools$info
table(Data$orig.ident)

# Data <- readRDS(paste0(info$dir$dir, info$dir$seurat, info$filename$qc))
# analysis
Data <- ElementaryPipline(Data,ngenes = 3000)
Data <- AdvancedAnalysis(
  Data, resolution = 0.8,
  # ref_singler_dir = "E:/Scripts/Rproject/scRNA_analysis/celldex_ref_dataset/",
  ref_singler_dir = "D:/Workspace/R/scRNA/celldex_ref_dataset/",
  ref_cell_dex = "HumanPrimaryCellAtlasData",
  ref_markers = NULL,
  ncl = 4,
  run_enrich = TRUE, run_TA = TRUE, run_CC = TRUE
)

Data <- readRDS("D:/Data/scRNA-seq/case/case_rds_files/case_results_seurat.rds")
Plot_seurat(Data = Data)

info <- Data@tools$info
rmarkdown::render(
  input = "scripts/report.Rmd",
  params = list(info = info),
  output_file = paste0(
    info$dir$dir,
    info$dir$results,
    info$data_name,
    "_report"
  )
)

# ShinyApp
scConf = createConfig(Data)
makeShinyApp(
  Data, 
  scConf, 
  gene.mapping = TRUE,
  shiny.title = "ShinyCell Quick Start",
  shiny.dir = paste0(info$dir$dir,info$dir$results, "shinyApp/")
) 

plot_heatmap(
  Data,
  genes = Data@assays$RNA@var.features[1:15],
  metaname = "singler_by_cluster_main"
)
DegTwo(Data,id1 = "T_cells",id2 = "Monocyte",metaname = "singler_by_cluster_main")
DegTwo(Data,id1 = "T_cells",metaname = "singler_by_cluster_main")

#X1.把所有的Data$seurat_clusters换成Idents(Data)
#√2.添加cellchat pig DB，若没有自动转为human
#√3.根据resolution1.5进行类注释
#X4.保存metadata.rds
#X5.enrich_table文件夹名称
#√6.var.feature 3000
