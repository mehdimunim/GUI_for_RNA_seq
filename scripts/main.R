setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# Main Script
source("read_data.R")
source("rpkm.R")
source("group.R")

## Genome files
annot_path <- "../data/genome/QM6aAnnotationIFPEN2021strict.gff"
mapping_table_path <- "../data/genome/MappingTable_geneID.xlsx"

## Load sample files paths and info
directory <- "../data/expression"
paths <- list.files(directory, pattern = "*.tsv", full.names = TRUE, recursive =TRUE)
replics <- c(rep(seq(1,6),4), 
             rep(seq(1,4),2), 
             rep(c(1,2), 2))
condition <- c(rep("glucose_24",6), 
               rep("glucose_48",6), 
               rep("lactose_24",6),
               rep("lactose_48",6), 
               rep("mixed_g75_l25_24",4), 
               rep("mixed_g75_l25_48",4), 
               rep("mixed_g90_l10_24",2),
               rep("mixed_g90_l10_48",2))

sampleFiles <- data.frame(paths, replics, condition)
names(sampleFiles) <- c("paths", "replicates", "conditions")


## Calculate and plot log(RPKM)
### Load reference genome info
length_table <- read_length_table(annot_path)
mapping_table <- read_mapping_table(mapping_table_path)

### Get RPKM vector for each sample files
sampleFiles$rpkm <- sapply(sampleFiles$paths, function(x) get_rpkm(x, length_table, mapping_table))

### group by experience
analysis <- group_by_analysis(sampleFiles)

### Plot RPKM comparison
apply(analysis, 1, function(x) plot_lrpkm_comparison(x$rpkm1, x$rpkm2, x$analysis))