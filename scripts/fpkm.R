fpkm <- function(condition_table) {
  # Get the FPKM calculation
  # RPKM seems deprecated (https://bioinformatics.stackexchange.com/a/69)
  counts <- as.numeric(condition_table$count)
  effective_lengths <- as.numeric(condition_table$gene_length)
  res <- exp(log(counts) - log(effective_lengths) - log(sum(counts) + log(1E9)))
  return(res)
}

get_fpkm <- function(file_path, length_table, mapping_table){ 
  ### read counts and group in two tables
  count_table <- read_count_table(file_path)
  table <-
    get_condition_table(count_table, length_table, mapping_table)
  l <- list(fpkm(table))
  return(l)
  }

plot_lfpkm <- function(fpkm,
                       title,
                       c1 = "lightgrey",
                       c2 = "red",
                       add = FALSE) {
  # Plot a Log fpkm distribution
  
  # Create histogram with log(fpkm)
  lfpkm <- log(fpkm)
  h <- hist(lfpkm,
            breaks = length(fpkm) / 100,
            plot = FALSE)
  
  # Create curve from the counts
  multiplier <- h$counts / h$density
  d <- density(lfpkm)
  d$y <- d$y * multiplier[1]
  
  # plot with graphical options
  plot(
    h,
    main = title,
    col = c1,
    xlab = "Log(fpkm)",
    ylab = "Count of log(fpkm)",
    add = add
  )
  lines(d, col = c2, lwd = 2)
}

plot_lfpkm_comparison <- function(fpkm1,
                                 fpkm2,
                                 title,
                                 save = FALSE,
                                 output_folder = "") {
  # Plot and compare fpkm distributions
  # if save is true, save graph in the output_folder
  plot_lfpkm(fpkm1, title, c1 = rgb(1,215/255,0,0.2), c2 = "red")
  plot_lfpkm(fpkm2,
            "",
            c1 = rgb(175/255,238/255,238/255,0.2),
            c2 = "blue",
            add = TRUE)
}