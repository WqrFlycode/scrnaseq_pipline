CellCommunication <- function(Data, by_cluster = "singler_by_cluster"){
  # if(!dir.exists(paste0(Data@tools$results_dir,"/CellCommunication/")))
  #   dir.create(paste0(Data@tools$results_dir,"/CellCommunication/"))
  # results_dir <- paste0(
  #   Data@tools$results_dir,
  #   "/CellCommunication/",
  #   Data@tools$data_name
  # )
  # 
  # net_dir <- paste0(results_dir,"_cellchat_",
  #                   c("net","netP","group_size"),".rds")
  Species <- Data@tools$info$Species
  
  cat("\n %%%%% create cellchat object %%%%% \n")
  data.input <- GetAssayData(Data, assay = "RNA", slot = "data")
  meta <- subset(Data@meta.data, select = by_cluster)
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = by_cluster)
  # cellchat <- createCellChat(object = Data, meta = Data@meta.data, group.by = by_cluster)
  
  cat("\n %%%%% input ligand-receptor database %%%%% \n")
  if(Species == "human"){
    cellchat@DB <- CellChatDB.human
  } else if(Species == "mouse"){
    cellchat@DB <- CellChatDB.mouse
  } else if(Species == "zebrafish"){
    cellchat@DB <- CellChatDB.zebrafish
  } else {
    stop("指定species")
  }
  
  cat("\n %%%%% subsetData %%%%% \n")
  cellchat <- subsetData(cellchat)
  cat("\n %%%%% identifyOverExpressedGenes %%%%% \n")
  cellchat <- identifyOverExpressedGenes(cellchat)
  cat("\n %%%%% identifyOverExpressedInteractions %%%%% \n")
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # project gene expression data onto PPI
  # (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  # cellchat <- projectData(cellchat,PPI.human)
  
  cat("\n %%%%% infer ligand-receptor level cell communication network %%%%% \n")
  cat("\n   %%%%% computeCommunProb %%%%% \n")
  cellchat <- computeCommunProb(cellchat, population.size = TRUE)
  cat("\n   %%%%% filterCommunication %%%%% \n")
  cellchat <- filterCommunication(cellchat,min.cells = 10)
  
  cat("\n %%%%% infer signal pathway level cell communication network %%%%% \n")
  cellchat <- computeCommunProbPathway(cellchat)
  
  cellchat <- aggregateNet(cellchat)
  
  # cellchat@data <- matrix()
  # cellchat@data.signaling <- matrix()
  # cellchat@meta <- data.frame()
  
  # # 细胞通讯网络圈图
  # png( 
  #   filename = paste0(results_dir,"_interactions.",figure_format),
  #   width = 1080*2,
  #   height = 1080,
  #   units = "px",
  #   bg = "white",
  #   res = 100)
  # par(mfrow = c(1,2),xpd = TRUE, mar = c(1,1,1,1))
  # netVisual_circle(cellchat_net$count,vertex.weight = group_size,
  #                  weight.scale = T,label.edge = F,title.name = "number of interactions")
  # netVisual_circle(cellchat_net$weight,vertex.weight = group_size,
  #                  weight.scale = T,label.edge = F,title.name = "Interaction weights/strength")
  # dev.off() # close
  # 
  # # signal from a type of cell
  # mat <- cellchat_net$count
  # png( 
  #   filename = paste0(results_dir,"_each_cell_interaction.",figure_format),
  #   width = 1080*3,
  #   height = 1080*3,
  #   units = "px",
  #   bg = "white",
  #   res = 450)
  # par(mfrow = c(ceiling(nrow(mat)/3),3), xpd=TRUE, mar = c(1,1,1,1))
  # for (i in 1:nrow(mat)) {
  #   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  #   mat2[i,] <- mat[i,]
  #   netVisual_circle(mat2, vertex.weight = group_size, weight.scale = TRUE, arrow.width = 0.2,
  #                    arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  # }
  # dev.off()
  # 
  # # hierarchy figure
  # pathways.show <- cellchat_netP$pathways[1]
  # png(
  #   filename = paste0(results_dir,"_hierarchy.",figure_format),
  #   width = 1920*2,
  #   height = 1920*(3/4)*2,
  #   units = "px",
  #   bg = "white",
  #   res = 450)
  # netVisual_aggregate(
  #   cellchat, 
  #   signaling = pathways.show, 
  #   vertex.receiver = 1:length(levels(Data$ref_markers_clusters)),
  #   layout = "hierarchy"
  # )
  # dev.off()
  
  cat("\n ----------CellCommunication finished---------- \n")
  return(cellchat)
}

# cellchat分析小鼠数据时的bug：
# Please check `object@data.signaling` and ensure that you have run `subsetData` and that the data matrix `object@data.signaling` looks