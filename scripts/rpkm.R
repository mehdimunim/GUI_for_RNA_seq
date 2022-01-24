if(!require(tidyverse))
  install.packages("tidyverse")
library(tidyverse)

rpkm <- function(condition_table) {
  # Get the average RPKM distributions
  # counts: number of reads per gene
  # lengths: lengths of the genes in counts
  # RPKM = numberOfReads / ( geneLength/1000 * totalNumReads/1,000,000 )
  counts <- as.numeric(condition_table$count)
  lengths <- as.numeric(condition_table$gene_length)
  rpk <- counts / (lengths / 1000)
  coef <- sum(rpk) / 1e6
  return(rpk / coef)
}

plot_lrpkm <- function(rpkm,
                       title,
                       c1 = "lightgrey",
                       c2 = "red",
                       add = FALSE) {
  # Plot a Log RPKM distribution
  
  # Create histogram with log(rpkm)
  lrpkm <- log(rpkm)
  h <- hist(lrpkm,
            breaks = length(rpkm) / 100,
            plot = FALSE)
  
  # Create curve from the counts
  multiplier <- h$counts / h$density
  d <- density(lrpkm)
  d$y <- d$y * multiplier[1]
  
  # plot with graphical options
  plot(
    h,
    main = title,
    col = c1,
    xlab = "Log(RPKM)",
    ylab = "Density of log(RPKM)",
    add = add
  )
  lines(d, col = c2, lwd = 2)
}

plot_lrpkm_comparison <- function(rpkm1,
                                 rpkm2,
                                 title,
                                 save = FALSE,
                                 output_folder = "") {
  # Plot and compare RPKM distributions
  # if save is true, save graph in the output_folder
  plot_lrpkm(rpkm1, title, c1 = rgb(1,215/255,0,0.2), c2 = "red")
  plot_lrpkm(rpkm2,
            "",
            c1 = rgb(175/255,238/255,238/255,0.2),
            c2 = "blue",
            add = TRUE)
}