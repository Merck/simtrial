# Helper functions used by test-double_programming_mb_weight.R

test_mb_weight <- function(x, delay = 4) {
  out <- NULL
  for (i in unique(x$stratum)) {
    outi <- x[x$stratum == i, ]
    outi.sort <- outi[order(outi$tte), ]
    locmaxt <- length(outi.sort$tte[outi.sort$tte <= delay]) # location of the maximum timepoint (tte) that is less or equal to the input 'delay'
    outi$mb_weight <- NA
    outi$mb_weight[1:locmaxt] <- 1 / outi$s[1:locmaxt]
    outi$mb_weight[(locmaxt + 1):nrow(outi)] <- outi$mb_weight[locmaxt]
    out <- rbind(out, outi)
  }

  return(out)
}
