setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# Main Script
source("read_data.R")
source("fpkm.R")
source("group.R")
source("test.R")

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


## Calculate and plot log(fpkm)
### Load reference genome info
length_table <- read_length_table(annot_path)
mapping_table <- read_mapping_table(mapping_table_path)

### Get fpkm vector for each sample files
sampleFiles$fpkm <- sapply(sampleFiles$paths, function(x) get_fpkm(x, length_table, mapping_table))

### group by experience
analysis <- group_by_analysis(sampleFiles)

### Plot fpkm comparison
apply(analysis, 1, function(x) plot_lfpkm_comparison(x$fpkm1, x$fpkm2, x$analysis))

## Student test
alpha <- 0.01
apply(analysis, 1, function(x) summarize_test(x$fpkm1, x$fpkm2, x$analysis, alpha))
