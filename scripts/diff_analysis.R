# Differential analysis
library(DESeq2)
library(ggplot2)

get_dds <-
  function(sampleFiles,
           sampleCondition) {
    sampleTable <- data.frame(sampleName = sampleFiles,
                              fileName = sampleFiles,
                              condition = sampleCondition)
    
    sampleTable$condition <- factor(sampleTable$condition)
    
    l <-
      lapply(sampleFiles, function(fn)
        read.table(fn,
                   fill = TRUE,
                   header = TRUE))
    tbl <- sapply(l, function(a)
      a[, ncol(a)])
    colnames(tbl) <- sampleFiles
    rownames(tbl) <- l[[1]]$V1
    rownames(sampleTable) <- sampleTable[, 1]
    ddsHTSeq <-
      DESeqDataSetFromMatrix(
        countData = tbl,
        colData = sampleTable[, -(1:2), drop = FALSE],
        design = ~ condition
      )
    ddsHTSeq <-
      DESeq(
        ddsHTSeq,
        test = "Wald",
        fitType = "mean",
        parallel = FALSE,
        BPPARAM = bpparam()
      )
    
    return (ddsHTSeq)
    
  }

get_diff_analysis <- function(dds,
                              pCutoff,
                              threshold,
                              cond.changed,
                              conda.ref)
{
  res <- results(
    dds,
    alpha = pCutoff,
    lfcThreshold = threshold,
    contrast = c("condition", cond.changed, conda.ref)
  )
  
  return (res)
}

extract_de_genes <- function(res, up.or.down) {
  # get the metadata from results
  p.cutoff <- metadata(res)$alpha
  threshold <- metadata(res)$lfcThreshold
  
  # subset of differentially expressed genes
  de <- subset(res, padj < p.cutoff)
  
  if (up.or.down == "up") {
    return(subset(de, log2FoldChange > threshold))
  }
  else {
    return(subset(de, log2FoldChange < -1 * threshold))
  }
  
}

plot_diff_results <- function(analysis_label, up, down) {
  ggplot(NULL,
         aes(x = c("low", "up"),
             y = c(length(up$padj), length(down$padj)))) +
    geom_col() +
    labs(title = analysis_label,
         x = "Differentially expressed genes",
         y = "Count of genes")
}
