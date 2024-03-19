# check required packages-------------------------------------------------------
require_packages <- c(
  "Seurat", "harmony", "SingleR", "dplyr", "monocle3", "CellChat",
  "ggplot2", "patchwork", "ggpubr", "ggrepel", "RColorBrewer",
  "clusterProfiler", "ReactomePA", "pathview", "org.Hs.eg.db"
)
required <- c(
  "4.3.0.1", "1.2.0", "2.2.0", "1.1.3", "1.3.4", "1.6.1",
  "3.4.3", "1.1.3", "0.6.0", "0.9.3", "1.1.3", "4.8.3",
  "1.44.0", "1.40.0", "3.17.0"
)
all_packages <- .packages(all.available = TRUE)
local <- sapply(require_packages, function(x) {
  if (x %in% all_packages) paste0(packageVersion(x)) else NA
})
packages_info <- data.frame(required, local)
rm(require_packages, required, all_packages, local)


# library-----------------------------------------------------------------------
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

# source------------------------------------------------------------------------
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