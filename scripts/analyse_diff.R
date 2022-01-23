# Analyse différentielle du champignon Trichoderma en utilisant DESeq2
# sources : https://bioconductor.riken.jp/packages/3.6/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## Clean memory
rm(list = ls())

library(DESeq2)

## Input directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory <- "../data"
annot_path <- "../QM6aAnnotationIFPEN2021strict.gff"
  
## Loading htseq-count files
sampleFiles <- list.files(directory, pattern = "(*_1[0-2])|(*_[1-9]).tsv")
sampleCondition <- c(rep("glucose", 6), rep("lactose", 6))
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition, levels=c("glucose", "lactose"))

# Calcuate RPKM
# RPKM = numberOfReads / ( geneLength/1000 * totalNumReads/1,000,000 )

rpkm <- function(counts, lengths) {
  rpk <- counts/(lengths/1000)
  coef <- sum(rpk) / 1e6
  rpk/coef
}  

read_counts <- function(file_path) {
  file <- read.table(file_path, sep ="\t")
  file
}

read_lengths <- function(annot_path) {
  file <- read.table(annot_path, sep ="\t")
  df <- data.frame(apply(file[9], 1, function(x) unlist(strsplit(x, "_|="))[1]), file[5] - file[4] + 1)
  aggregate(df[2], by=df[1], FUN=sum) 
}

rpkm(read_counts(paste0("../data/",sampleFiles[1])), read_lengths(annot_path))