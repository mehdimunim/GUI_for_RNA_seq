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
        colData = sampleTable[,-(1:2), drop = FALSE],
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

# plot_diff_results <- function(analysis_label, up, down) {
#   ggplot(NULL,
#          aes(x = c("low", "up"),
#              y = c(length(up$padj), length(down$padj)))) +
#     geom_col() +
#     labs(title = analysis_label,
#          x = "Differentially expressed genes",
#          y = "Count of genes")
# }

get_up_down_table <- function(dds,
                              expTable,
                              p.cutoff,
                              threshold) {
  m <- apply(expTable,
             1,
             function(comp) {
               # extract number of up and down genes
               res <-
                 get_diff_analysis(dds, p.cutOff, threshold, comp[["treated"]], comp[["ref"]])
               up <- extract_de_genes(res, "up")
               down <- extract_de_genes(res, "down")
               return (c(length(up$padj), length(down$padj)))
             })
  m <- t(m)
  m <- data.frame(m)
  names(m) <- c("up", "down")
  m$label <- rownames(expTable)
  
  # Reformat data frame
  up_down_table <- data.frame(
    count =
      c(m$up, m$down),
    type = c(rep("up", length(m$up)), rep("down", length(m$down))),
    label = rep(m$label,2)
             )
  
  
  return(up_down_table)
}

plot_diff_results <-
  function(up_down_table, fontsize) {
    
    ggplot(up_down_table,
           aes(
             x = type,
             y = count)
           ) +
      geom_col() +
      facet_wrap(~label) +
      theme(strip.text = element_text(size = fontsize))
  }
