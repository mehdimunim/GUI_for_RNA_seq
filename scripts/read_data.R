if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)
read_count_table <- function(file_path, mapping_table) {
  counts <- read.table(file_path, sep ="\t")
  colnames(counts) <- c("gene_id", "count")
  # ID JGI (Joint Genome Institute)
  return(counts)
}

read_mapping_table <- function(file_path) {
  mapping_table <- read.xlsx(file_path)
  return(mapping_table)
}

read_length_table <- function(annot_path) {
  # get the lengths of the genes in the annot_path
  # annot_path: path of the GFF annotation file
  # The rows of the GFF correspond to exons
  data <- read.table(annot_path, sep ="\t")
  # Reduce ID field (column 9) to the gene ID
  gene_id <- apply(data[9], 1, function(x) unlist(strsplit(x, "\\Q|\\E|="))[2])
  gene_length <- data[5] - data[4] + 1
  lengths <- data.frame(gene_id, gene_length)
  # Now, we have several lengths for the same gene_id because of several exons
  # Therefore, we sum up the exon length to get the gene length
  df <- aggregate(gene_length, by=list(gene_id), FUN=sum)
  colnames(df) <- c("gene_id", "gene_length")
  return(df)
}

get_condition_table <- function(count_table, length_table, mapping_table) {
  count_extended <- merge(count_table, mapping_table, by.x = "gene_id", by.y = "ID-JGI")
  full_table <- merge(length_table, count_extended, by.x = "gene_id", by.y = "ID-IFPEN")
  return(full_table)
}