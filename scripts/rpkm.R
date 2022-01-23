if(!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)

rpkm <- function(condition_table) {
  # Get the average RPKM distributions 
  # counts: number of reads per gene
  # lengths: lengths of the genes in counts
  # RPKM = numberOfReads / ( geneLength/1000 * totalNumReads/1,000,000 )
  counts <- as.numeric(condition_table$count)
  lengths <- as.numeric(condition_table$gene_length)
  rpk <- counts/(lengths/1000)
  coef <- sum(rpk) / 1e6
  return(rpk/coef)
}  

plot_rpkm <- function(rpkm, title, save = FALSE, output_folder ="", c = "grey") {
  # Plot an RPKM distribution
  # if save is true, save graph in the output_folder
  hist(log(rpkm), 
       main = title,
       col = c,
       xlab = "Log(RPKM)",
       ylab = "Density of log(RPKM)",
       breaks = length(rpkm)/100,
       freq = FALSE
       )
  lines(density(log(rpkm)), col = "red", lwd = 3)
}

plot_rpkm_comparison <- function(rpkm1, rpkm2, title, save = FALSE, output_folder ="") {
  # Plot and compare RPKM distributions
  # if save is true, save graph in the output_folder
  df1 <- data.frame(rpkm1)
  df2 <- data.frame(rpkm2)
  df1$id <- 1
  df2$id <- 2
  data <- rbind(df1, df2)
  ggplot(data) +
    aes(1)
    geom_histogram(bins = 100, col ="blue") +
    geom_histogram(data = df2, bins = 100, col ="red") +
    labs(title = title)
  if (save) { 
    filename <- "rpkm_output_" + title + ".pdf"
    ggsave(filename, path = output_folder, dpi=500)
  }
}