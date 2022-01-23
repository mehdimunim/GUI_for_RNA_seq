# Analyse différentielle du champignon Trichoderma en utilisant DESeq2
# sources : https://bioconductor.riken.jp/packages/3.6/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## Clean memory
rm(list = ls())

library(DESeq2)

## Input directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory <- "../data"
annot_path <- "../data/genome/QM6aAnnotationIFPEN2021strict.gff"
mapping_table_path <- "../data/genome/MappingTable_geneID.xlsx"

## Loading htseq-count files
sampleFiles <- list.files(directory, pattern = "(*_1[0-2])|(*_[1-9]).tsv")
sampleCondition <- c(rep("glucose", 6), rep("lactose", 6))
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition, levels=c("glucose", "lactose"))

# Calcuate RPKM

rpkm(read_counts(paste0("../data/",sampleFiles[1])), read_lengths(annot_path))