if(!require(dplyr)) install.packages(dplyr)
library(dplyr) 
group_by_analysis <- function(fpkm_table) {
  names <- c(
    "Lactose vs Glucose 24h",
    "Lactose vs Glucose 48h",
    "90% Glucose, 10% Lactose vs Glucose 24h",
    "90% Glucose, 10% Lactose vs Glucose 48h",
    "75% Glucose, 25% Lactose vs Glucose 24h",
    "75% Glucose, 25% Lactose vs Glucose 48h"
  )
  res <- data.frame(names,
                    cbind(sapply(
                      c(
                        "lactose_24",
                        "lactose_48",
                        "mixed_g90_l10_24",
                        "mixed_g90_l10_48",
                        "mixed_g75_l25_24",
                        "mixed_g75_l25_48"
                      ),
                      FUN = function(x) avg(fpkm_table, x)
                    )),
                    cbind(sapply(
                      c(
                        "glucose_24",
                        "glucose_48",
                        "glucose_24",
                        "glucose_48",
                        "glucose_24",
                        "glucose_48"
                      ),
                      FUN = function(x) avg(fpkm_table, x)
                    )))
  names(res) <- c("analysis", "fpkm1", "fpkm2")
  return(res)
}

avg <- function(fpkm_table, condition_label) {
  fpkms <-
    filter(fpkm_table,  condition == condition_label)$fpkm # select fpkm by condition
  res <- rowMeans(data.frame(fpkms)) # average condition
  res <- list(res)
  return(res)
}

get_condition_comparisons <- function() {
  
  res <- data.frame(treated = c(
                        "lactose_24",
                        "lactose_48",
                        "mixed_g90_l10_24",
                        "mixed_g90_l10_48",
                        "mixed_g75_l25_24",
                        "mixed_g75_l25_48"
                      ),
                    ref = c(
                        "glucose_24",
                        "glucose_48",
                        "glucose_24",
                        "glucose_48",
                        "glucose_24",
                        "glucose_48"
                      ))
  rownames(res) <- c(
    "Lactose vs Glucose 24h",
    "Lactose vs Glucose 48h",
    "90% Glucose, 10% Lactose vs Glucose 24h",
    "90% Glucose, 10% Lactose vs Glucose 48h",
    "75% Glucose, 25% Lactose vs Glucose 24h",
    "75% Glucose, 25% Lactose vs Glucose 48h"
  )
  
  return(res)
}


