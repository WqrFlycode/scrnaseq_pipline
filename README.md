# support for species
Annotation:

  "human", "rat"

enrich:

  GO:
  
  Reactome: "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"

CellChat:


<!--

# Issues

## cellchat::netVisual_circle

Error in i_set_edge_attr(x, attr(value, "name"), index = value, value = attr(value,  : 
  Length of new attribute value must be 1 or 79, the number of target edges, not 75

solve: 

The package igaph 1.4.0 is not compatible with cellchat. When I changed igraph 1.4.0 to 1.3.5, the problem was solved.

-->