# support for species
Annotation: human, mouse

enrich:

  GO:
  
  Reactome: "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"

CellChat: human, mouse, zebrafish


<!--

# Issues

## cellchat::netVisual_circle

Error in i_set_edge_attr(x, attr(value, "name"), index = value, value = attr(value,  : 
  Length of new attribute value must be 1 or 79, the number of target edges, not 75

solve: 

The package igaph 1.4.0 is not compatible with cellchat. When I changed igraph 1.4.0 to 1.3.5, the problem was solved.

-->


# future
1.seurat v5
2.analyze and plot step by step
3.使用相对路径，便于移动目录
4.save meta data

<!--
# read data
## data_dir:    单细胞数据文件路径
## filename:    10X的三个文件名的前缀
## data_name:   数据名称，默认为"case"，名称不宜过长
## results_dir: 输出结果的路径，默认NULL, 结果保存在数据文件夹中
## 注意：
## 1.第一次读取数据后seurat对象被保存在rds文件中，data_dir和results_dir被保存在seurat对象中
## 2.所有路径中使用“/”，最后一级路径后面也要带“/”
## 3.文件名称应为*matrix.mtx.gz; *barcodes.tsv.gz; *
-->
