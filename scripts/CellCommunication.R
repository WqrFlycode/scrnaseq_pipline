CellCommunication <- function(Data, by_cluster = "seurat_clusters", 
                              re = FALSE, figure_format = "png"){
  if(!dir.exists(paste0(Data@tools$results_dir,"/CellCommunication/")))
    dir.create(paste0(Data@tools$results_dir,"/CellCommunication/"))
  results_dir <- paste0(
    Data@tools$results_dir,
    "/CellCommunication/",
    Data@tools$data_name
  )
  
  net_dir <- paste0(results_dir,"_cellchat_",
                    c("net","netP","group_size"),".rds")
  if (FALSE) {
    cellchat_net <- readRDS(net_dir[1])
    cellchat_netP <- readRDS(net_dir[2])
    group_size <- readRDS(net_dir[3])
  }else{
    # create cellchat object
    cellchat <- createCellChat(object = Data, group.by = by_cluster)
    
    # input 配体-受体 database
    data("CellChatDB.human")
    
    # select a subset
    CellChatDB.use <- subsetDB(
      CellChatDB.human, 
      search = unique(CellChatDB.human$interaction$annotation)[1]
    )
    cellchat@DB <- CellChatDB.use
    
    # preprocess
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat,PPI.human)
    
    # infer 配体-受体 level cell communication network
    cellchat <- computeCommunProb(cellchat,raw.use = FALSE,population.size = TRUE)
    cellchat <- filterCommunication(cellchat,min.cells = 10)
    
    # infer signal pathway level cell communication network
    cellchat <- computeCommunProbPathway(cellchat)
    
    cellchat <- aggregateNet(cellchat)
    group_size <- as.numeric(table(cellchat@idents))
    cellchat_net <- cellchat@net
    cellchat_netP <- cellchat@netP
    saveRDS(cellchat_net, net_dir[1])
    saveRDS(cellchat_netP, net_dir[2])
    saveRDS(group_size, net_dir[3])
  }
  
  # 细胞通讯网络圈图
  png( 
    filename = paste0(results_dir,"_interactions.",figure_format),
    width = 1080*2,
    height = 1080,
    units = "px",
    bg = "white",
    res = 100)
  par(mfrow = c(1,2),xpd = TRUE, mar = c(1,1,1,1))
  netVisual_circle(cellchat_net$count,vertex.weight = group_size,
                   weight.scale = T,label.edge = F,title.name = "number of interactions")
  netVisual_circle(cellchat_net$weight,vertex.weight = group_size,
                   weight.scale = T,label.edge = F,title.name = "Interaction weights/strength")
  dev.off() # close
  
  # signal from a type of cell
  mat <- cellchat_net$count
  png( 
    filename = paste0(results_dir,"_each_cell_interaction.",figure_format),
    width = 1080*3,
    height = 1080*3,
    units = "px",
    bg = "white",
    res = 450)
  par(mfrow = c(ceiling(nrow(mat)/3),3), xpd=TRUE, mar = c(1,1,1,1))
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i,] <- mat[i,]
    netVisual_circle(mat2, vertex.weight = group_size, weight.scale = TRUE, arrow.width = 0.2,
                     arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()
  
  # hierarchy figure
  pathways.show <- cellchat_netP$pathways[1]
  png(
    filename = paste0(results_dir,"_hierarchy.",figure_format),
    width = 1920*2,
    height = 1920*(3/4)*2,
    units = "px",
    bg = "white",
    res = 450)
  netVisual_aggregate(
    cellchat, 
    signaling = pathways.show, 
    vertex.receiver = 1:length(levels(Data$ref_markers_clusters)),
    layout = "hierarchy"
  )
  dev.off()
  
  print("----------CellCommunication finished----------")
  return(Data)
}