if(!require(dplyr)) install.packages(dplyr)
library(dplyr) 
group_by_analysis <- function(RPKM_table) {
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
                        "glucose_24",
                        "glucose_48",
                        "mixed_g90_l10_24",
                        "mixed_g90_l10_48",
                        "mixed_g75_l25_24",
                        "mixed_g75_l25_48"
                      ),
                      FUN = function(x) avg(RPKM_table, x)
                    )),
                    cbind(sapply(
                      c(
                        "lactose_48",
                        "lactose_48",
                        "glucose_24",
                        "glucose_48",
                        "glucose_24",
                        "glucose_48"
                      ),
                      FUN = function(x) avg(RPKM_table, x)
                    )))
  names(res) <- c("analysis", "rpkm1", "rpkm2")
  return(res)
}

avg <- function(RPKM_table, condition_label) {
  rpkms <-
    filter(RPKM_table,  condition == condition_label)$rpkm # select rpkm by condition
  res <- rowMeans(data.frame(rpkms)) # average condition
  res <- list(res)
  return(res)
}



