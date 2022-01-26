test_difference <- function(fpkm1, fpkm2) {
  # Hypothesis: we suppose lfpkm1 and lfpkm2 follow binomial laws
  # The variance are unknown but equal
  # We want to determine if the two mean are significantly different
  r1 <- log(unlist(unname(fpkm1)))
  r2 <- log(unlist(unname(fpkm2)))
  # avoid -Inf values 
  r1 <- r1[!is.infinite(r1)]
  r2 <- r2[!is.infinite(r2)]
  test <- t.test(r1, r2,var.equal = TRUE, null.value = -Inf)
  return(test)
}

is_test_significant <- function(test, alpha) {
  return (test$p.value < alpha)
}

print_test_summary <- function(test, conclusion, analysis_label) {
  print(analysis_label)
  print("Two Sample t-test")
  print("H0: true difference in means is equal to 0")
  print(paste("p-value:",round(test$p.value,8)))
  if(conclusion == FALSE) {
    print("Unsignificant average expression")
  } 
  else {
    print("Significant average expression")
  }
  
}

summarize_test <- function(fpkm1, fpkm2, analysis_label, alpha) {
  test <- test_difference(fpkm1, fpkm2)
  conclusion <- is_test_significant(test, alpha)
  print_test_summary(test, conclusion, analysis_label)
}

