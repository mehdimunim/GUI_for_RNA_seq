library(tidyverse)

clustering <- function(expMatrix, n.clusters) {
  kmeans.obj <- kmeans(expMatrix, n.clusters)
  return(kmeans.obj)
  
}
get_expression_matrix <- function(dds,
                                  expTable,
                                  p.cutoff,
                                  threshold) {
  m <- apply(expTable,
             1,
             function(comp) {
               res <-
                 get_diff_analysis(dds, p.cutOff, threshold, comp[["treated"]], comp[["ref"]])
               return (res$log2FoldChange)
             })
  m <- data.frame(m)
  m$gene_id <- rownames(m)
  return(m)
}

filter_exp <- function(expMatrix) {
  perc <- 100*sum(is.na(expMatrix))/nrow(expMatrix)
  print(paste0(round(perc,2)," % of NA in expression matrix"))
  print("Filtering NA")
  cleaned_expMatrix <-na.omit(expMatrix)
  return (cleaned_expMatrix)
}



plot_cluster_sizes <- function(expMatrix, kmeans.object) {
  exp_cluster_obj <- cbind(expMatrix, cluster = kmeans.object$cluster)
  ggplot(exp_cluster_obj, 
         aes(x = cluster, 
             fill = cluster))+
    geom_bar()
}

plot_exp_profile <- function(expMatrix) {
  # Inspired from: https://stackoverflow.com/a/63273303
  as.data.frame(expMatrix) %>%
    pivot_longer(cols=1:6, names_to="comparison", values_to="logFC") %>%
    ggplot(aes(x=comparison, y=logFC)) +
    geom_line(aes(group=gene_id), size=0.5, alpha=0.3, color="black")
    }

plot_mean_exp_profile <- function(expMatrix, kmeans.object) {
  
}