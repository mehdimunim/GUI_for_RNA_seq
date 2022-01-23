if(!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)

rpkm <- function(condition_table) {
  # Get the average RPKM distributions 
  # counts: number of reads per gene
  # lengths: lengths of the genes in counts
  # RPKM = numberOfReads / ( geneLength/1000 * totalNumReads/1,000,000 )
  counts <- condition_table$count
  lengths <- condition_table$length
  rpk <- counts/(lengths/1000)
  coef <- sum(rpk) / 1e6
  return(rpk/coef)
}  

plot_rpkm <- function(rpkm, title) {
  # Plot an rpkm distribution
  ggplot(rpkm) +
    aes(x) +
    geom_point() +
    labs(title = title)
}

plot_rpkm_comparison <- function(rpkm, rpkm2) {
  
  
}