# 分析
suppressMessages(library(Seurat))
# suppressMessages(library(harmony))
suppressMessages(library(SingleR))
# suppressMessages(library(dplyr))
suppressMessages(library(monocle3))
# suppressMessages(library(SeuratWrappers))
suppressMessages(library(CellChat))
# 作图
# suppressMessages(library(ggplot2))
# suppressMessages(library(patchwork))
# suppressMessages(library(ggpubr))
# suppressMessages(library(ggrepel))
# suppressMessages(library(RColorBrewer))
# 富集
suppressMessages(library(clusterProfiler))
suppressMessages(library(ReactomePA))
# suppressMessages(library(pathview))
# suppressMessages(library(org.Hs.eg.db))

# library(ShinyCell)

# function_dir：函数的R脚本文件路径
function_dir <- "scripts/"
source(paste0(function_dir,"/ReadData.R"))
source(paste0(function_dir,"/Integration.R"))
source(paste0(function_dir,"/Annotation.R"))
source(paste0(function_dir,"/Enrich.R"))
source(paste0(function_dir,"/TrajectoryAnalysis.R"))
source(paste0(function_dir,"/CellCommunication.R"))
source(paste0(function_dir,"/analysis_scrnaseq.R"))
source(paste0(function_dir,"/Plots.R"))
# rm(function_dir)
