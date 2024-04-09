plot_enrich <- function(enrichment_dir, results_dir) {
  enrichment <- readRDS(enrichment_dir)
  for(dbn in names(enrichment)) {
    db_results <- enrichment[[dbn]]
    for(cl in names(db_results)) {
      cl_results <- db_results[[cl]]
      nr <- cl_results@result %>%
        filter(pvalue < 0.05, p.adjust < 0.05) %>%
        nrow()
      if(nr > 0) {
        ## 条形图，按p从小到大排，绘制前20个Term
        bar_plot <- barplot(
          cl_results, showCategory = 20, title = paste0(dbn,"_bar")
        )
        ## 点图，按富集的数从大到小的
        dot_plot <- dotplot(
          cl_results, showCategory = 20, title = paste0(dbn,"_dot")
        )
        if(dbn %in% names(enrichment)[1:3] & nrow(cl_results) > 1) {
          acyclic_plot <- goplot(cl_results)
          enrich_plot <- (bar_plot + dot_plot)/acyclic_plot
          enrich_plot_h <- 18
        } else {
          enrich_plot <- bar_plot + dot_plot
          enrich_plot_h <- 10
        }
        ggsave(
          paste0(results_dir, "_", cl, "_", dbn, ".svg"),
          enrich_plot,
          width = 18, height = enrich_plot_h, # scale = 3,
          bg = "white"
        )
      }
    }
  }
  return(enrich_plot)
}  
