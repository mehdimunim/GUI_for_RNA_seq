# Main Script 
rm(list = ls())
source("read_data.R")
source("rpkm.R")

## Input directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory <- "../data"
annot_path <- "../data/genome/QM6aAnnotationIFPEN2021strict.gff"
mapping_table_path <- "../data/genome/MappingTable_geneID.xlsx"

## Loading count file paths
#sampleFiles <- list.files(directory, pattern = "(*_1[0-2])|(*_[1-9]).tsv")

## Calculate RPKM
### Reference genome
length_table <- read_length_table(annot_path)
mapping_table <- read_mapping_table(mapping_table_path)

### RNA-seq reads
path1 <- paste(directory,"expression/glucose/24/GSM2188043_expression_1.tsv",sep ="/")
path2 <- paste(directory,"expression/glucose/48/GSM2188049_expression_19.tsv",sep ="/")

### read counts and group in two tables
count_table1 <- read_count_table(path1)
count_table2 <- read_count_table(path2)
table_24 <- get_condition_table(count_table1, length_table, mapping_table)
table_48 <- get_condition_table(count_table2, length_table, mapping_table)

### RPKM
rpkm_24 <- rpkm(table_24)
rpkm_48 <- rpkm(table_48)

#plot_rpkm(rpkm_24, "RPKM test glucose 24 48 rep1", c = "grey")
plot_rpkm_comparison(rpkm_24, rpkm_48, "RPKM test glucose 24 48 rep1")


