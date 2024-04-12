# library-----------------------------------------------------------------------
# 分析
suppressMessages(library(Seurat))
# suppressMessages(library(harmony))
suppressMessages(library(SingleR))
# suppressMessages(library(dplyr))
suppressMessages(library(monocle3))
# suppressMessages(library(SeuratWrappers))
suppressMessages(library(CellChat))
suppressMessages(library(parallel))
# suppressMessages(library(findPC))
suppressMessages(library(clustree))
# suppressMessages(library(cluster))
# suppressMessages(library(tidyr))
# suppressMessages(library(writexl))

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

suppressMessages(library(ShinyCell))

# source------------------------------------------------------------------------
function_dir <- "scripts/"
source(paste0(function_dir,"/ReadData.R"))
source(paste0(function_dir,"/QualityControl.R"))
source(paste0(function_dir,"/Integration.R"))
source(paste0(function_dir,"/Annotation.R"))
source(paste0(function_dir,"/Enrich.R"))
source(paste0(function_dir,"/TrajectoryAnalysis.R"))
source(paste0(function_dir,"/CellCommunication.R"))
source(paste0(function_dir,"/analysis_scrnaseq.R"))
source(paste0(function_dir,"/Plots.R"))
function_dir <- "Plots/"
source(paste0(function_dir,"plot_enrich.R"))
source(paste0(function_dir,"plot_qc.R"))
# rm(function_dir)


# install packages--------------------------------------------------------------
# 安装SeuratWrappers
# 适应于4.3版的Seurat
# remotes::install_github('satijalab/seurat-wrappers@community-vignette')

# 安装cellchat
# devtools::install_github("jokergoo/ComplexHeatmap")
# devtools::install_github("immunogenomics/presto")
# devtools::install_github("jinworks/CellChat")

# 安装orgdb
# https://bioconductor.org/packages/release/BiocViews.html#___OrgDb 
# BiocManager::install("org.Mm.eg.db")