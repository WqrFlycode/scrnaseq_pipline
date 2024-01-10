```
ps <- c("Matrix","Seurat","ggplot2","patchwork","patchwork","ggpubr","ggrepel","SingleR",
  "dplyr","org.Hs.eg.db","clusterProfiler","pathview","ReactomePA","monocle3","SeuratWrappers",
  "CellChat")

n=length(ps)
vid=rep(NA,n)

for (i in 1:n) {
  vid[i]=paste0("r-",ps[i],"=",packageVersion(ps[i]))
}
paste0(vid,collapse = " ")
```
# 安装SeuratWrappers
适应于4.3版的Seurat
```
remotes::install_github('satijalab/seurat-wrappers@community-vignette')
```

# 安装cellchat
```
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("immunogenomics/presto")
devtools::install_github("jinworks/CellChat")
```
