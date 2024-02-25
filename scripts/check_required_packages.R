require_packages <- c(
  "Seurat", "harmony", "SingleR", "dplyr", "monocle3", "CellChat",
  "ggplot2", "patchwork", "ggpubr", "ggrepel", "RColorBrewer",
  "clusterProfiler", "ReactomePA", "pathview", "org.Hs.eg.db"
)
required <- c("4.3.0.1", "1.2.0", "2.2.0", "1.1.3", "1.3.4", "1.6.1",
                      "3.4.3", "1.1.3", "0.6.0", "0.9.3", "1.1.3", "4.8.3",
                      "1.44.0", "1.40.0", "3.17.0")
all_packages <- .packages(all.available = TRUE)
local <- sapply(require_packages, function(x) {
  if (x %in% all_packages) paste0(packageVersion(x)) else NA
})
data.frame(required, local)
