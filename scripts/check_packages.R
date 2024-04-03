# check required packages-------------------------------------------------------
require_packages <- c(
  "Seurat", "harmony", "SingleR", "dplyr", "monocle3", "SeuratWrappers",
  "CellChat", "parallel", "findPC", "clustree", "cluster", "tidyr", "writexl",
  "ggplot2", "patchwork", "ggpubr", "ggrepel", "RColorBrewer",
  "clusterProfiler", "ReactomePA", "pathview", "org.Hs.eg.db",
  "ShinyCell"
)
required <- c(
  "4.3.0.1", "1.2.0", "2.4.0", "1.1.4", "1.3.4", "0.3.1",
  "1.1.0", "4.3.1", "1.0", "0.5.1", "2.1.4", "1.3.0", "1.5.0",
  "3.4.4", "1.1.3", "0.6.0", "0.9.4", "1.1.3",
  "4.10.0", "1.46.0", "1.42.0", "3.18.0",
  "2.1.0"
)
all_packages <- .packages(all.available = TRUE)
local <- sapply(require_packages, function(x) {
  if (x %in% all_packages) paste0(packageVersion(x)) else NA
})
packages_info <- data.frame(required, local)
rm(require_packages, required, all_packages, local)
print(packages_info)
