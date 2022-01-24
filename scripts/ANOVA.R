# Analysis of variance


test_difference <- function(rpkm1, rpkm2, analysis_label) {
  # Use ANOVA to test difference between two conditions
  df <- data.frame(as.data.frame(rpkm1), as.data.frame(rpkm2))
  names(df) <- c("rpkm1", "rpkm2")
  one.way <- aov(rpkm1 ~ rpkm2, data = df) 
  return(one.way)
}