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

# 安装orgdb
# https://bioconductor.org/packages/release/BiocViews.html#___OrgDb 
```
BiocManager::install("org.Mm.eg.db")
```

# 各部分程序物种支持
Annotation:
enrich:
  GO:
  Reactome: "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"
CellChat:

# Issues

## cellchat::netVisual_circle

Error in i_set_edge_attr(x, attr(value, "name"), index = value, value = attr(value,  : 
  Length of new attribute value must be 1 or 79, the number of target edges, not 75

solve: 

The package igaph 1.4.0 is not compatible with cellchat. When I changed igraph 1.4.0 to 1.3.5, the problem was solved.